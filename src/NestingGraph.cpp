#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>
#include "common/nfp_lib.hpp"
#include "common/geometry_util.hpp"
#include "common/graph.hpp"
#include <utility>

#if defined(_WIN32)
  #include <windows.h>
#elif defined(__APPLE__)
  #include <mach-o/dyld.h>
#else //Linux and Unix-like systems
  #include <unistd.h>
  #include <limits.h>
#endif

using json = nlohmann::json;
using namespace std;
namespace fs = std::filesystem;



std::filesystem::path executable_path() {
#if defined(_WIN32)
    std::wstring buf;
    buf.resize(32768); // max path for many Win32 APIs
    DWORD len = GetModuleFileNameW(nullptr, buf.data(), (DWORD)buf.size());
    if (len == 0) return {};
    buf.resize(len);
    return std::filesystem::path(buf);

#elif defined(__APPLE__)
    uint32_t size = 0;
    _NSGetExecutablePath(nullptr, &size);          // get required size
    std::vector<char> buf(size);
    if (_NSGetExecutablePath(buf.data(), &size) != 0) return {};
    // _NSGetExecutablePath may return a path with symlinks; canonicalize if you want
    return std::filesystem::weakly_canonical(std::filesystem::path(buf.data()));

#else // Linux and many Unix-like systems with /proc
    std::vector<char> buf(4096);
    for (;;) {
        ssize_t n = readlink("/proc/self/exe", buf.data(), buf.size());
        if (n < 0) return {};
        if ((size_t)n < buf.size()) {
            return std::filesystem::path(std::string(buf.data(), (size_t)n));
        }
        buf.resize(buf.size() * 2);
    }
#endif
}


bool is_existing_directory(const std::string& s) {
    if (s.find('\0') != std::string::npos) return false;
    std::error_code ec;
    return std::filesystem::is_directory(std::filesystem::path(s), ec) && !ec;
}


void writeNfpInfpJson(const std::string& outputdir,
                      const nlohmann::json& polygons)  // polygons can be any JSON value (array/object)
{
    // if not os.path.exists(outputdir): os.mkdir(outputdir)
    if (!fs::exists(outputdir)) {
        fs::create_directories(outputdir); // like mkdir, but also creates parents if needed
    }

    // nfp_infp = json.dumps(polygons, indent=2)
    const std::string nfp_infp = polygons.dump(2);

    // with open(outputdir+'/nfp-infp.json', 'w') as file: file.write(...)
    const fs::path outPath = fs::path(outputdir) / "nfp-infp.json";
    std::ofstream out(outPath, std::ios::out | std::ios::trunc);
    if (!out) {
        throw std::runtime_error("Failed to open file for writing: " + outPath.string());
    }
    out << nfp_infp;
}

struct BBox {
    double xMin, xMax, yMin, yMax;
};

using PtIJ16   = std::array<int16_t, 2>;   // (i, j) grid indices
using PatternIJ = std::vector<PtIJ16>;

using PtD      = std::array<double, 2>;    // (x, y) world coords (double)
using PatternPos = std::vector<PtD>;

using NfpIndexTableT = std::unordered_map<std::string,
                        std::unordered_map<std::string, PatternIJ>>;

using NfpPosTableT   = std::unordered_map<std::string,
                        std::unordered_map<std::string, PatternPos>>;

static inline void buildNfpPatternsForPoly(
    const std::string& key,
    nlohmann::json& poly,                 // must be non-const because we mutate nfp vertices by pivot shifting
    const nlohmann::json& pivot,          // {x,y}
    double gx, double gy,
    NfpIndexTableT& nfpIndexTable,
    NfpPosTableT& nfpPosTable
) {
    if (!poly.contains("nfps")) return;

    const double px = pivot.at("x").get<double>();
    const double py = pivot.at("y").get<double>();

    for (auto& nfp : poly.at("nfps")) {

        // --- Compute bbox while pivot-zeroing nfp vertices ---
        double vxMax = std::numeric_limits<double>::lowest();
        double vxMin = std::numeric_limits<double>::max();
        double vyMax = std::numeric_limits<double>::lowest();
        double vyMin = std::numeric_limits<double>::max();

        auto& verts = nfp.at("VERTICES");
        for (auto& v : verts) {
            v["x"] = v.at("x").get<double>() - px;
            v["y"] = v.at("y").get<double>() - py;

            const double x = v.at("x").get<double>();
            const double y = v.at("y").get<double>();

            vxMax = (x > vxMax) ? x : vxMax;
            vxMin = (x < vxMin) ? x : vxMin;
            vyMax = (y > vyMax) ? y : vyMax;
            vyMin = (y < vyMin) ? y : vyMin;
        }

        const auto polyVerts = NFPTool::parseVertices(verts);

        const int iStart = static_cast<int>(std::floor(vxMin / gx)) - 1;
        const int iEnd   = static_cast<int>(std::ceil (vxMax / gx)) + 3;
        const int jStart = static_cast<int>(std::floor(vyMin / gy)) - 1;
        const int jEnd   = static_cast<int>(std::ceil (vyMax / gy)) + 3;

        // Local buffers (freed each iteration after move)
        PatternIJ  pattern;
        PatternPos pattern_pos;

        // reserve rough capacity to reduce reallocations
        // (bbox cells) = (iEnd-iStart)*(jEnd-jStart)
        const std::size_t cells = static_cast<std::size_t>(std::max(0, iEnd - iStart))
                                * static_cast<std::size_t>(std::max(0, jEnd - jStart));

        pattern.reserve(cells);
        pattern_pos.reserve(cells);

        for (int i = iStart; i < iEnd; ++i) {
            for (int j = jStart; j < jEnd; ++j) {
                const double tx = static_cast<double>(i) * gx;
                const double ty = static_cast<double>(j) * gy;

                if (GeometryUtil::pointInPolygon({tx, ty}, polyVerts) == 1) {
                    pattern.push_back({static_cast<int16_t>(i), static_cast<int16_t>(j)});
                    pattern_pos.push_back({tx, ty});
                }
            }
        }

        const std::string other = nfp.at("POLYGON").get<std::string>();
        nfpIndexTable[key][other] = std::move(pattern);
        nfpPosTable[key][other]   = std::move(pattern_pos);
    }
}

struct GenPoint {
    double x;
    double y;
    uint64_t id;
};

// board: JSON array of {"x":..., "y":...}
static inline std::vector<GenPoint>
generatePoints(double length,double width, double freq)
{
    if (freq == 0.0) {
        throw std::invalid_argument("error: freq must be non-zero");
    }

    // Python: stepsx = int((maxx-minx)/freq)
    // Python int() truncates toward 0, not floor.
    const int stepsx = static_cast<int>(length / freq);
    const int stepsy = static_cast<int>(width / freq);

    std::vector<GenPoint> points;
    points.reserve(static_cast<std::size_t>(stepsx + 1) * static_cast<std::size_t>(stepsy + 1));

    for (int i = 0; i <= stepsx; ++i) {
        for (int j = 0; j <= stepsy; ++j) {
            points.push_back(GenPoint{
               (static_cast<double>(i) * freq),
               (static_cast<double>(j) * freq),
               0
            });
        }
    }

    return points;
}




struct LayerPoints {
    string polygon;
    unsigned int layer = 0;
    pair<uint64_t, uint64_t> indexRange; // (start_index, end_index)
    vector<GenPoint> innerFitPoints; // {x, y, accId}
};

struct LayerMatrix {
    string polygon;
    int layer = 0;
    int nx = 0; // int(length/gx)
    int ny = 0; // int(width/gy)
    vector<uint64_t> data; // (nx+1)*(ny+1), filled with UINT64_MAX

    uint64_t& at(int ix, int iy) {
        return data[static_cast<size_t>(ix) * static_cast<size_t>(ny + 1)
                  + static_cast<size_t>(iy)];
    }
};

using LayerPolyT = vector<pair<int, string>>;

struct LayersResult {
    LayerPolyT layerPoly;
    vector<LayerPoints> layerOfPoint;
    vector<LayerMatrix> layerMatrix;
    uint64_t valIndex = 0;
};

static inline LayersResult generateLayers(
    const json& polygons,
    int freq,
    double gx,
    double gy,
    double length,
    double width,
    bool type_oriented
) {
    LayersResult out;

    int layer = 0;
    uint64_t accIndex = 0;

    const int nx = static_cast<int>(length / gx);
    const int ny = static_cast<int>(width  / gy);

    for (auto& [key, val] : polygons.items()) {
        const std::string polyKey = key;
        const json& mainPiece = val;

        unsigned int quantity = mainPiece.value("QUANTITY", 1);
        if (type_oriented && quantity > 1) quantity = 1;

        const GeometryUtil:: Polygon innerfit = NFPTool::parseVertices(mainPiece.at("innerfit"));
        const double innerFitArea = std::abs(GeometryUtil::polygonArea(innerfit));
        const double boardArea = length * width;
        for (int q = 0; q < quantity; ++q) {
            auto boardPoints = generatePoints(length,width, freq);
            
            
            LayerPoints lp;
            lp.polygon = polyKey;
            lp.layer = layer;
            lp.innerFitPoints.reserve(boardArea / innerFitArea + 16); //TODO: estimate better using area ratio
            
            LayerMatrix lm;
            lm.polygon = polyKey;
            lm.layer = layer;
            lm.nx = nx;
            lm.ny = ny;
            lm.data.assign(static_cast<std::size_t>(nx + 1) * static_cast<std::size_t>(ny + 1),
                           std::numeric_limits<uint64_t>::max()); // == -1 in uint64

            lp.indexRange.first = accIndex;

            for (const auto& p : boardPoints) {
                const double x = p.x;
                const double y = p.y;

                if (GeometryUtil::pointInPolygon({static_cast <double>(x), static_cast <double>(y)}, innerfit) != 0) {
                    lp.innerFitPoints.push_back(GenPoint{x, y, accIndex});

                    // Python: innerMatrix[int(x/gx)][int(y/gy)] = AccIndex
                    // Here x,y are uint64; gx,gy are double. If gx/gy are integers, cast once.
                    const int ix = static_cast<int>(static_cast<double>(x) / gx);
                    const int iy = static_cast<int>(static_cast<double>(y) / gy);

                    if (ix >= 0 && ix <= nx && iy >= 0 && iy <= ny) {
                        lm.at(ix, iy) = accIndex;
                    }
                    ++accIndex;
                }
            }
            lp.indexRange.second = accIndex-1;

            out.layerOfPoint.push_back(std::move(lp));
            out.layerMatrix.push_back(std::move(lm));
            out.layerPoly.push_back({layer, polyKey});
            ++layer;
        }
        out.valIndex = accIndex;
    }

    return out;
}

void processGroupPair(
    std::vector<GenPoint> PointsA,
    PatternIJ& NFPA_B,
    LayerMatrix& innerMatrixB,
    BBox BBoxB,
    double gx,
    double gy,
    Graph& graph
) {
    for (const auto& pA : PointsA) {
        const unsigned int idxOffset = static_cast<int>(pA.x / gx);
        const unsigned int idyOffset = static_cast<int>(pA.y / gy);
        for (const auto& nfp_pt : NFPA_B) {
            const int16_t i_offset = nfp_pt[0]+ idxOffset;
            const int16_t j_offset = nfp_pt[1]+ idyOffset;

            const double xB = (static_cast<double>(i_offset) * gx);
            const double yB = (static_cast<double>(j_offset) * gy);

            // check if (xB,yB) is within BBoxB
            if (xB < BBoxB.xMin || xB > BBoxB.xMax || yB < BBoxB.yMin || yB > BBoxB.yMax) {
                continue;
            }

            // get grid indices for innerMatrixB
            const int ixB = static_cast<int>(xB / gx);
            const int iyB = static_cast<int>(yB / gy);

            if (ixB < 0 || ixB > innerMatrixB.nx || iyB < 0 || iyB > innerMatrixB.ny) {
                continue;
            }

            const uint64_t idB = innerMatrixB.at(ixB, iyB);
            if (idB != std::numeric_limits<uint64_t>::max()) {
                // valid point found in B
                graph.addEdge(pA.id, idB);
            }
        }
    }
}

string findFileinDirectory(const string& directory, const string& filename){
    error_code ec;
    for (const fs::directory_entry& e : fs::directory_iterator(directory,ec)) {
        if (e.path().filename() == filename) {
            return e.path().string();
        }
    }
    return "";
}

//TODO: change datasetPath to datasetPaths
vector<json> loadJson(string jsonPath,unordered_map<string,string> datasetPathMap,bool width_set=false,double width=0,bool length_set=false,double length=0){
    vector<json> datasets;
    std::ifstream file(jsonPath);
    if (!file) {
        std::cerr << "Failed to open " << jsonPath << "\n";
        exit(1);
    }

    json dataset;
    try{
        
        file >> dataset;
    }
    catch (const nlohmann::json::parse_error& e) {
        cerr << "Error: Failed to parse JSON file " << jsonPath << ": " << e.what() << "\n";
        exit(1);
    }
    // case1: json is an object containing a single dataset/instance
    if (dataset.is_object()){

        //check if is an instance
        if (dataset.contains("dataset") && dataset.contains("outputName") && dataset.contains("width") && dataset.contains("length")){
            string datasetName = dataset["dataset"].get<string>();
            try{
                dataset["dataset"] = datasetPathMap.at(datasetName);
            }
            catch(const out_of_range&){
                 if(fs::exists(dataset["dataset"].get<string>()) && fs::is_regular_file(dataset["dataset"].get<string>())){
                    //do nothing, the path is correct
                 }else{
                    cerr << "Error: dataset name " << datasetName << " not found in provided dataset directories.\n";
                    exit(1);
                }
            }
           
            datasets.push_back(dataset);
            return datasets;
        }
        

        //Make a wrapper json instance
        json wrapper;
        wrapper["dataset"] = jsonPath;
        wrapper["outputName"] = filesystem::path(jsonPath).stem().string();
        //set width and length

        if(dataset.contains("width")){
            wrapper["width"] = dataset["width"].get<double>();
            dataset.erase("width");
        }
        else{
            if (!width_set){
                cerr << "Error: width not specified in json or as argument.\n";
                exit(1);
            }
            wrapper["width"] = width;
        }
        
        if(dataset.contains("length")){
            wrapper["length"] = dataset["length"].get<double>();
            dataset.erase("length");
        }
        else{
            if (!length_set){
                cerr << "Error: length not specified in json or as argument.\n";
                exit(1);
            }
            wrapper["length"] = length;
        }
        wrapper["quantity"] = json::object();
        for (auto& [key, val] : dataset.items()){
            if(!val.contains("VERTICES") || !val["VERTICES"].is_array()){
                cerr << "Error: polygon " << key << " missing VERTICES array.\n";
                exit(1);
            }
            wrapper["quantity"][key] = dataset[key].value("QUANTITY",1);
        }

        datasets.push_back(wrapper);
        return datasets;
    }
    else if (dataset.is_array()){
        //case2: json is an array of instances
        for (auto& set : dataset){
                //check if is an instance
                if (set.contains("dataset") && set.contains("outputName") && set.contains("width") && set.contains("length")){
                    string datasetName = set["dataset"].get<string>();
                try{
                     set["dataset"] = datasetPathMap.at(datasetName);
                }
                catch(const out_of_range&){
                    //try to see if the dataset is the path instead of filename
                    if(fs::exists(set["dataset"].get<string>()) && fs::is_regular_file(set["dataset"].get<string>())){
                    //do nothing, the path is correct
                    }else{
                        cerr << "Error: dataset name " << datasetName << " not found in provided dataset directories.\n";
                        exit(1);
                    }
                }
                datasets.push_back(set);
            }
            else{
                cerr << "Error: one of the dataset in the array is missing required fields.\n";
                exit(1);
            }
        }
        return datasets;
    }
    else{
        cerr << "Error: JSON root must be an object or an array.\n";
        exit(1);
    }
}

static bool parse_bool(const char* s) {
    std::string v(s);
    if (v == "1" || v == "true" || v == "TRUE" || v == "True") return true;
    if (v == "0" || v == "false" || v == "FALSE" || v == "False") return false;
    throw std::invalid_argument("Expected boolean (0/1/true/false), got: " + v);
}

int main(int argc, char** argv) {
    auto dir_path = executable_path().parent_path();
    const string Usagephrase = "--instances <file|dir> Loads one instance JSON file, or a directory of instance JSON files. repeat for multiple file/directory nested directory will be ingnored.\n"
               "--datasets <dir> directory of dataset repeat for multiple directory\n"
               "--outputdir <dir>[optiional default ./result] \n"
               "--typeOriented <0|1> [optional, default false] \n"
               "--cliqueCovering <0|1> [optional, default false] \n"
               "--width <number> [optional, mandatory if dataset doenst contain it] ignored if specified in json \n"
               "--length <number>[optional, mandatory if dataset doenst contain it]  ignored if specified in json\n";
    if (argc < 2) {
        cerr
            << "Usage: " << argv[0]
            << Usagephrase;
        return 1;
    }

    
    bool type_oriented = false; //one of the arguments
    vector<string> inputPaths; //one of the arguments
    vector<string> datasetPaths; //one of the arguments
    string outputdir = ""; //one of the arguments
    bool cliqueCovering = false; //one of the arguments
    json dataset;
    double width = 0.0; //one of the arguments
    double length = 0.0; //one of the arguments
    

    bool type_oriented_param = false;
    bool cliqueCovering_param = false;
    bool length_set = false;
    bool width_set = false;
    
    for (int argi = 1; argi < argc; argi += 2) {
        if (argi + 1 >= argc) {
            std::cerr << "Missing value for argument: " << argv[argi] << "\n";
            return 1;
        }

        std::string key = argv[argi];

        if (key == "--instances") {
            inputPaths.push_back(argv[argi + 1]);
        } else if (key == "--outputdir") {
            outputdir = argv[argi + 1];
        } else if (key == "--type_oriented") {
            type_oriented = parse_bool(argv[argi + 1]);
            type_oriented_param = true;
        } else if (key == "--cliqueCovering") {
            cliqueCovering = parse_bool(argv[argi + 1]);
            cliqueCovering_param = true;
        } else if (key == "--width") {
            width = std::stod(argv[argi + 1]);
            width_set = true;
        } else if (key == "--length") {
            length = std::stod(argv[argi + 1]);
            length_set = true;
        } else if (key == "--datasets") {
            datasetPaths.push_back(argv[argi + 1]);
        } else {
            std::cerr << "Unknown argument: " << key << "\n";
            std::cerr << "Usage:" << argv[0] << " \n" << Usagephrase;
            return 1;
        }
    }
    
    if (datasetPaths.empty())
    {
        cerr << "ERROR: Dataset path is required. Use -datasetPath <dataset path>\n";
        return 1;
    }
    

    if (inputPaths.empty()){
        cerr << "ERROR: Input path is required. Use -input <instance json path> <datasest json path> <dataset/instance folder>\n";
        return 1;
    }

    if (outputdir.empty()){
        cerr << "Warning: Output directory not specified, using default ./results \n";
        outputdir = dir_path.string() + "/results/";
    }

    if (!type_oriented_param){
        cerr << "Warning: type_oriented parameter not set, defaulting to false.\n";
    }

    if (!cliqueCovering_param){
        cerr << "Warning: cliqueCovering parameter not set, defaulting to false.\n";
    }

    

    unordered_map<string, string> datasetPathMap; //checking duplicates dataset names and fast insert dataset paths
    for (const auto& datasetPath : datasetPaths){
        if (!is_existing_directory(datasetPath)){
            cerr << "Error: dataset path is not a valid folder directory: " << datasetPath << "\n";
            return 1;
        }
        for (const auto& entry : fs::directory_iterator(datasetPath)){
            if (entry.path().extension() == ".json"){
                string datasetJsonPath = entry.path().string();
                string datasetName = entry.path().stem().string();
                
                if (datasetPathMap.find(datasetName) != datasetPathMap.end()){
                    cerr << "Error: Duplicate dataset name found: " << datasetName << " at " 
                         << datasetJsonPath << "\n" 
                         << datasetPathMap.at(datasetName) << "\n";
                    return 1;
                }
                datasetPathMap[datasetName] = datasetJsonPath;
                datasetPathMap[entry.path().filename().string()] = datasetJsonPath;
            }
        }
    }

    vector <json> datasets;
    // Load dataset json
    //case1: inputPath is a file path
    error_code ec;
    for (const auto& inputPath : inputPaths){
        if (filesystem::is_regular_file(inputPath, ec)) {
            ifstream in(inputPath);
            if (!in) {
                throw runtime_error("Failed to open input file: " + inputPath);
            }
            // Load dataset json
            try {
                vector<json> tmpJson = loadJson(inputPath,datasetPathMap,width_set,width,length_set,length);
                datasets.insert(datasets.end(), tmpJson.begin(), tmpJson.end());

            } catch (const exception& e) {
                cerr << "Error parsing JSON file(613): " << e.what() << "\n";
                cerr << "File: " << inputPath << "\n";
                return 1;
            }
        }
        // case2: inputPath is a directory path
        else if (filesystem::is_directory(inputPath, ec)) { 
            
            try {
                //iterate through all json files in the directory
                unordered_map<string, int> outputNameMap; //for checking duplicate dataset names
                for (const fs::directory_entry& e : fs::directory_iterator(inputPath)) {      
                    vector<json> tmpJson = loadJson(e.path().string(),datasetPathMap,width_set,width,length_set,length);
                    for (auto &set : tmpJson){
                        string name = set["outputName"].get<string>();
                        if (outputNameMap.find(name) != outputNameMap.end()){
                            cerr << "Error: Duplicate dataset name found: " << name << "\n";
                            return 1;
                        }
                        outputNameMap[name] = 1;
                    }
                    datasets.insert(datasets.end(), tmpJson.begin(), tmpJson.end());
                }
            } catch (const exception& e) {
                cerr << "Error parsing JSON file: " << e.what() << "\n";
                cerr << "File: " << inputPath << "\n";
                return 1;
            }
        }
        else {
            cerr << "Error: input path is not a valid file or directory: " << inputPath << "\n";
            return 1;
        }
    }
    

    for (auto &set : datasets){

        string selected =  set["dataset"].get<string>();
        //if (!is_existing_directory(selected)){
        //    cerr << "dataset path is not a valid folder directory: " << selected << "\n";
        //    return 1;
        //}

        string outputname = set["outputName"].get<string>();
        double width = set["width"].get<double>();
        double length = set["length"].get<double>();
        const string outputDataset = outputdir + outputname;
        unsigned int step = 1;
        unsigned int gx = step;
        unsigned int gy = step;

        error_code ec;
        std::ifstream file(selected);
        if (!file) throw std::runtime_error("Cannot open: " + selected);
        json dataset = json::parse(file);

        //remove board info from dataset to prevent interference with NFP processing and polygon iteration
        dataset.erase("rect");
        dataset.erase("length");
        dataset.erase("width");


        json polygons;
        try {
            polygons = NFPTool::processNFP(dataset, length, width);
        } catch (const exception& e) {
            cerr << "Error during NFP processing: " << e.what() << "\n";
            return 1;
        }

        // write the dataset with NFP into JSON file

        writeNfpInfpJson(outputDataset, polygons);

        nlohmann::json board = polygons.at("rect");
        polygons.erase("rect"); 
        size_t num_polygon = polygons.size();

        //quantity can either be a single unsigned int for all polygons or a dictionary of per-polygon quantities

        if (set.contains("quantity")){
            if(set["quantity"].is_number_unsigned()) {
                for (auto& [key, val] : polygons.items()){
                    polygons[key]["QUANTITY"] = set["quantity"].get<unsigned int>();
                }
            }
            else{
                for (auto& [key, val] : set["quantity"].items()){
                    polygons[key]["QUANTITY"] = set["quantity"][key].get<unsigned int>();
                }
            }
        }

        int total_polygon = 0;
        double total_area = 0.0;

        for (auto& [key, val] : polygons.items()){
            total_polygon += polygons[key]["QUANTITY"].get<unsigned int>();
            total_area += polygons[key]["QUANTITY"].get<double>() * std::abs(GeometryUtil::polygonArea(NFPTool::parseVertices(polygons[key]["VERTICES"])));
        }

        NfpPosTableT NfpPosTable;
        NfpIndexTableT NfpIndexTable;

        for (auto& [key, poly] : polygons.items()){

            json& pivot = poly.at("VERTICES").at(0);
            json& innerfit = poly.at("innerfit");

            //pivot the polygon
            for (auto& v : polygons.at(key).at("VERTICES")) {
                v["x"] = v.at("x").get<double>() - pivot.at("x").get<double>();
                v["y"] = v.at("y").get<double>() - pivot.at("y").get<double>();
            } 

            //Compute innerfit BBox
            BBox Innerbbox{ numeric_limits<double>::max(),  //xMin
                            numeric_limits<double>::lowest(), //xMax
                             numeric_limits<double>::max(),  //yMin
                             numeric_limits<double>::lowest() //yMax
                    };
            
            for (auto& v : innerfit) {
                const double x = v.at("x").get<double>();
                const double y = v.at("y").get<double>();
                if (x > Innerbbox.xMax) Innerbbox.xMax = x;
                if (x < Innerbbox.xMin) Innerbbox.xMin = x;
                if (y > Innerbbox.yMax) Innerbbox.yMax = y;
                if (y < Innerbbox.yMin) Innerbbox.yMin = y;
            }

            poly["innerfit_BoundingBox"]["xMin"] = floor(Innerbbox.xMin / gx);
            poly["innerfit_BoundingBox"]["xMax"] = ceil(Innerbbox.xMax / gx);
            poly["innerfit_BoundingBox"]["yMin"] = floor(Innerbbox.yMin / gy);
            poly["innerfit_BoundingBox"]["yMax"] = ceil(Innerbbox.yMax / gy);

            //Generate NFP Pattern 

            
            
            buildNfpPatternsForPoly(key,poly,pivot,gx,gy,NfpIndexTable,NfpPosTable);
            
            /*
            for (auto& nfp : poly.at("nfps")) {
                std::cout << "  nfp original: " << nfp.at("VERTICES").dump() << "\n";

                double vxMax = std::numeric_limits<double>::lowest();
                double vxMin =  std::numeric_limits<double>::max();
                double vyMax = std::numeric_limits<double>::lowest();
                double vyMin =  std::numeric_limits<double>::max();

                // for v in nfp['VERTICES']: v -= Pivot; update min/max
                for (auto& v : nfp.at("VERTICES")) {
                    v["x"] = v.at("x").get<double>() - pivot.at("x").get<double>();
                    v["y"] = v.at("y").get<double>() - pivot.at("y").get<double>();

                    const double x = v.at("x").get<double>();
                    const double y = v.at("y").get<double>();

                    if (x > vxMax) vxMax = x;
                    if (x < vxMin) vxMin = x;
                    if (y > vyMax) vyMax = y;
                    if (y < vyMin) vyMin = y;
                }

                
                const int iStart = floor(vxMin / gx) - 1;
                const int iEnd   = ceil(vxMax / gx) + 3;
                const int jStart = floor(vyMin / gy) - 1;
                const int jEnd   = ceil(vyMax / gy) + 3;

                for (int i = iStart; i < iEnd; ++i) {
                    for (int j = jStart; j < jEnd; ++j) {
                        const double tx = static_cast<double>(i) * gx;
                        const double ty = static_cast<double>(j) * gy;

                        // if key == 'PIECE 4' and nfp['POLYGON'] == 'PIECE 5': print testpoint
                        if (key == "PIECE 4" && nfp.value("POLYGON", "") == "PIECE 5") {
                            cout << " testing point: {\"x\":" << tx << ",\"y\":" << ty << "}\n";
                        }
                        
                        // res = point_in_polygon(testpoint, nfpVertices)
                        if (GeometryUtil::pointInPolygon(GeometryUtil::Point{tx,ty}, NFPTool::parseVertices(nfp.at("VERTICES"))) == 1) {
                            // nfp_pattern.append([i*gx,j*gy,i,j]) then cast to int16
                            // NOTE: i*gx and j*gy are doubles; Python stores them, then np.int16 casts.
                            // Here we mimic by truncating toward 0 before int16.
                            const int16_t ii = static_cast<int16_t>(i);
                            const int16_t jj = static_cast<int16_t>(j);
                            
                            pattern.push_back({ii, jj});
                            pattern_pos.push_back({tx, ty});
                        }
                    }
                }
                
                NfpIndexTable[key][nfp.at("POLYGON").get<std::string>()] = move(pattern);
                NfpPosTable[key][nfp.at("POLYGON").get<std::string>()] = move(pattern_pos);

            }
            */
        }
        //generate Layer points (vertex)

        LayersResult layers = generateLayers(
            polygons,
            step,
            static_cast<double>(gx),
            static_cast<double>(gy),
            static_cast<double>(length),
            static_cast<double>(width),
            type_oriented
        );
        // Build graph

        Graph graph(static_cast<int>(layers.valIndex)); // valIndex is 1 + max vertex id
        for (const auto& [layerA, polyA] : layers.layerPoly) {
            const auto& pointsA = layers.layerOfPoint[layerA].innerFitPoints;

            for (const auto& [layerB, polyB] : layers.layerPoly) {
                if (layerA == layerB)
                {
                    if(!type_oriented)
                    {
                        continue;
                    }
                }

                auto& nfpPatternB = NfpIndexTable.at(polyA).at(polyB);
                if (nfpPatternB.empty()) continue;

                // find innerMatrixB
                auto& innerMatrixB = layers.layerMatrix[layerB];

                // find BBoxB
                const auto& bboxBJson = polygons.at(polyB).at("innerfit_BoundingBox");
                BBox BBoxB{
                    bboxBJson.at("xMin").get<double>(),
                    bboxBJson.at("xMax").get<double>(),
                    bboxBJson.at("yMin").get<double>(),
                    bboxBJson.at("yMax").get<double>()
                };

                processGroupPair(
                    pointsA,
                    nfpPatternB,
                    innerMatrixB,
                    BBoxB,
                    static_cast<double>(gx),
                    static_cast<double>(gy),
                    graph
                );
            }
        }

        // Optionally compute clique covering
        EdgeCliqueCover max1Cover;
        EdgeCliqueCover max1MinEKCover;
        uint64_t NFPEdges = graph.getNumEdges();
        uint32_t NumberVertices = graph.getNumVertices();
        

        if (cliqueCovering){

            if(!type_oriented){
                //Add clique edges for each layer if not type oriented.
                for (const auto& [layer, poly] : layers.layerPoly) {
                    const uint64_t start = layers.layerOfPoint[layer].indexRange.first;
                    const uint64_t end = layers.layerOfPoint[layer].indexRange.second;
                    for (uint64_t i = start; i < end; ++i) {
                        for (uint64_t j = i + 1; j <= end; ++j) {
                            graph.addEdge(i, j);
                        }
                    }
                }
            }
            // Implement clique covering algorithm
            cout << "Starting clique covering: Max1\n";
            max1Cover = maximum1Heuristic(graph);
            //cout << "Start "
            max1MinEKCover = expandKouHeuristic(graph, max1Cover);
        }

        // Compute graph statistics
        uint64_t cliqueEdges = 0;
        if(!type_oriented){
            for (const auto& [layer, poly] : layers.layerPoly) {
                const uint64_t start = layers.layerOfPoint[layer].indexRange.first;
                const uint64_t end = layers.layerOfPoint[layer].indexRange.second;
                const uint64_t numPoints = end - start + 1;
                cliqueEdges += (numPoints * (numPoints - 1)) / 2;
            }
        }
        double density = (2.0 * (NFPEdges+cliqueEdges)) / (NumberVertices * (NumberVertices - 1));

        cout << "Graph statistics:\n";
        cout << "  Dataset name:     " << outputname << "\n";
        cout << "  Number of vertices: " << graph.getNumVertices() << "\n";
        cout << "  Number of NFP edges: " << NFPEdges << "\n";
        cout << "  Number of clique edges:   " << cliqueEdges << "\n";
        cout << "  Number of edges:    " << NFPEdges+cliqueEdges << "\n";
        cout << "  Density:            " << density << "\n";
        if(cliqueCovering){
            cout << "  Number of cliques (clique covering): " << max1MinEKCover.size() << "\n";
        }
        cout << "Finished processing dataset: " << outputname << "\n";
        // Output graph to file
        if(cliqueCovering){
            const string cliqueOutputdir = outputDataset + "/edgecover.txt";
            max1MinEKCover.writeIntoFile(cliqueOutputdir);
        }
        else{
            const string graphOutputPath = outputDataset + "/graph.csv";
            graph.writeEdgesToFile(graphOutputPath);
        }

        
        

        // Write pointsCoordinate
        const string pointsOutputPath = outputDataset + "/pointCoordinate.txt";
        ofstream pointsOut(pointsOutputPath, ios::out | ios::trunc);
        if (!pointsOut) {
            throw runtime_error("Failed to open file for writing: " + pointsOutputPath);
        }
        pointsOut << "##format: (Layer,x,y,id)\n";
        for (const auto& lp : layers.layerOfPoint) {
            for (const auto& p : lp.innerFitPoints) {
                pointsOut << lp.layer << "," << p.x << "," << p.y << "," << p.id << "\n";
            }
        }
        
        //Write layerpolygon
        const string LayerPoly = outputDataset + "/LayerPoly.txt";
        ofstream LayerPolyOut(LayerPoly, ios::out | ios::trunc);
        if (!LayerPolyOut) {
            throw runtime_error("Failed to open file for writing: " + LayerPoly);
        }
        LayerPolyOut << "##format:  Layer polygon\n";
        for (const auto& lp : layers.layerOfPoint) {
            LayerPolyOut << lp.layer << "\t" << lp.polygon << "\n";
        }

        //Write Cliques
        const string cliquesOutputPath = outputDataset + "/cliques.csv";
        ofstream cliquesOut(cliquesOutputPath, ios::out | ios::trunc);
        if (!cliquesOut) {
            throw runtime_error("Failed to open file for writing: " + cliquesOutputPath);
        }
        cliquesOut << "##format: (start_index, end_index)\n";
        for (const auto& lp : layers.layerOfPoint) {
            cliquesOut << lp.indexRange.first << "," << lp.indexRange.second << "\n";
        }

        //Write metadata
        const string metadataOutputPath = outputDataset + "/metadata.csv";
        ofstream metadataOut(metadataOutputPath, ios::out | ios::trunc);
        if (!metadataOut) {
            throw runtime_error("Failed to open file for writing: " + metadataOutputPath);
        }
        

        metadataOut <<"Name :\t" <<outputname << "\n";
        metadataOut <<"Total Pieces :\t" << total_polygon << "\n";
        metadataOut <<"Type of Pieces :\t" << num_polygon << "\n";
        metadataOut <<"Board Width :\t" << width << "\n";
        metadataOut <<"Board Length :\t" << length << "\n";
        metadataOut <<"Number of Nodes :\t" << graph.getNumVertices() << "\n";
        metadataOut <<"Number of Edges :\t" << graph.getNumEdges() << "\n";
        if(!type_oriented){  
            metadataOut <<"Number of Clique Edges :\t" << cliqueEdges << "\n";
        }
        else{
            metadataOut <<"Number of Clique Edges :\t" << -1 << "\n";
        }

        if(cliqueCovering){
            metadataOut <<"Clique Cover :\t" << max1MinEKCover.size() << "\n";
        }
        
        metadataOut <<"Total Polygon Area :\t" << total_area << "\n";
        metadataOut.close();

    }
    return 0;
}