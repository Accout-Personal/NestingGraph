#include "common/nfp_lib.hpp"
#include "common/geometry_util.hpp"
#include "common/graph.hpp"
#include "common/cancel_exit.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <charconv>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>
#include <utility>
#include <limits>
#include <algorithm>
#include <filesystem>
#include <system_error>

using namespace std;

namespace fs = std::filesystem;

bool REMOVE_OLD_FOLDER = true;
bool CLIQUE_COVERING_ADD_LAYER_CLIQUE = true;
bool OUTPUT_ADD_LAYER_CLIQUE = true;
bool OUT_OLD_METADATA = false;


std::string to_string_fixed(double x, int decimals) {
    std::array<char, 128> buf{};
    auto [ptr, ec] = std::to_chars(buf.data(), buf.data() + buf.size(),
                                   x, std::chars_format::fixed, decimals);
    if (ec != std::errc{}) throw std::runtime_error("to_chars failed");
    return std::string(buf.data(), ptr);
}

// argv0 = argv[0] dal main
fs::path executable_path_from_argv0(const char* argv0) {
    if (!argv0 || !*argv0) return {};
    std::error_code ec;

    fs::path p = fs::path(argv0);

    // Se è relativo, rendilo assoluto rispetto alla working dir
    if (p.is_relative())
        p = fs::absolute(p, ec);

    // Canonicalizzazione "soft" (non fallire se non esiste)
    fs::path c = fs::weakly_canonical(p, ec);
    return ec ? p : c;
}


bool is_existing_directory(const std::string& s) {
    if (s.find('\0') != std::string::npos) return false;
    std::error_code ec;
    return std::filesystem::is_directory(std::filesystem::path(s), ec) && !ec;
}

 // polygons can be any JSON value (array/object)
void writeNfpInfpJson(const std::string& outputdir,
                      const nlohmann::json& polygons){

    // nfp_infp = json.dumps(polygons, indent=2)
    const std::string nfp_infp = polygons.dump(2);

    // with open(outputdir+'/nfp-infp.json', 'w') as file: file.write(...)
    const fs::path outPath = fs::path(outputdir);
    std::ofstream out(outPath, std::ios::out | std::ios::trunc);
    if (!out) {
        throw std::runtime_error("Failed to open file for writing: " + outPath.string());
    }
    out << nfp_infp;
}

GeometryUtil::BBox computeBoundingBox(json poly){
    GeometryUtil::Polygon geoPoly;
    geoPoly.reserve(poly.size());
    for (auto const & v: poly){
        geoPoly.push_back(GeometryUtil::Point{v["x"].get<double>(), v["y"].get<double>()});
    }
    return GeometryUtil::computeBoundingBox(geoPoly);

}



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

        // --- pivot-zeroing nfp vertices ---
        auto& verts = nfp.at("VERTICES");
        for (auto& v : verts) {
            v["x"] = v.at("x").get<double>() - px;
            v["y"] = v.at("y").get<double>() - py;
        }

        //Compute Bounding Box of NFP
        GeometryUtil::BBox nfpBBox = computeBoundingBox(verts);

        const auto polyVerts = GeometryUtil::parseVertices(verts);

        const int iStart = static_cast<int>(std::floor(nfpBBox.xMin / gx)) - 1;
        const int iEnd   = static_cast<int>(std::ceil (nfpBBox.xMax / gx)) + 3;
        const int jStart = static_cast<int>(std::floor(nfpBBox.yMin / gy)) - 1;
        const int jEnd   = static_cast<int>(std::ceil (nfpBBox.yMax / gy)) + 3;

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

struct PointCoord {
    unsigned int layer;
    double x;
    double y;
    uint64_t id;
};

// board: JSON array of {"x":..., "y":...}
static inline std::vector<GenPoint>
generatePoints(double length,double width, double freq){
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
    vector<vector<unsigned int>> Ymap;
    double MaxX = 0;
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
    bool type_oriented){

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

        const GeometryUtil:: Polygon innerfit = GeometryUtil::parseVertices(mainPiece.at("innerfit"));
        const double innerFitArea = std::abs(GeometryUtil::polygonArea(innerfit));
        
        const double boardArea = length * width;
        for (int q = 0; q < quantity; ++q) {
            auto boardPoints = generatePoints(length,width, freq);
            
            
            LayerPoints lp;
            lp.polygon = polyKey;
            lp.layer = layer;
            //std::cout << "Board Area: " << boardArea << " InnerfitArea " << innerFitArea << endl;
            lp.innerFitPoints.reserve(static_cast<unsigned int>(innerFitArea + 16)); //TODO: estimate better using area ratio
            
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
                    lp.MaxX = x > lp.MaxX ? x : lp.MaxX; 
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

            //Make Y MAP
            //lp.Ymap = vector<vector <unsigned int>>();
            int MaxY = 0;
            vector<GenPoint> LayerCoords = lp.innerFitPoints;
            sort(LayerCoords.begin(), LayerCoords.end(),
            [](const GenPoint& a, const GenPoint& b) {
                if (a.y != b.y) return a.y < b.y;  // Sort by y value increasing order
                return a.x < b.x; // then by x increasing  order.
            });
            //now the max should be the last row
            //unsigned int MaxY = LayerCoords.back().y;
            lp.Ymap.resize(LayerCoords.back().y - LayerCoords.front().y + 1);
            MaxY = 0;
               
            /*for (auto& row : LayerCoords) {
                cout << row[0] << " " << row[1] << " " << row[2] << " "<< row[3] << " \n";
            }*/
               
            for (unsigned int j = 0 ; j< LayerCoords.size(); j++)
            {
                 if (LayerCoords[j].y > MaxY)
                 {
                     MaxY = LayerCoords[j].y;
                 }
                 //cout << "pushing back Y: "<< MaxY << " X:" << LayerCoords[j][1] << "\n";
                 lp.Ymap[MaxY - LayerCoords.front().y].push_back(LayerCoords[j].id);

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
    GeometryUtil::BBox BBoxB,
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
            if(!set.contains("dataset"))
            {
                cerr << "Error: dataset : "<<set.dump(0) <<"is missing \"dataset\" field" << endl;
                exit(1);
            }

            if(!set.contains("outputName"))
            {
                cerr << "Error: dataset : "<<set.dump(0) <<"is missing \"outputName\" field" << endl;
                exit(1);
            }

            if(!set.contains("width"))
            {
                cerr << "Error: dataset : "<<set.dump(0) <<"is missing \"width\" field" << endl;
                exit(1);
            }

            if(!set.contains("length"))
            {
                cerr << "Error: dataset : "<<set.dump(0) <<"is missing \"length\" field" << endl;
                exit(1);
            }

            
            //check if is an instance
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

unordered_map<string,string> MapDuplicateDataset(json &polydataset){
    vector<string> keylist;
    unordered_map<string,GeometryUtil::Polygon> polygons;
    for (auto [k,v]:polydataset.items()){
        polygons[k] = GeometryUtil::parseVertices(v["VERTICES"]);
        keylist.push_back(k);
    }
    unordered_map<string,string> finalMap;
    vector<string> removeList;
    for (auto polyA:keylist){
        for (auto polyB:keylist){
            if (polyA == polyB) continue;
            if (GeometryUtil::samePolygonVertices(polygons[polyA],polygons[polyB]) && abs(GeometryUtil::polygonArea(polygons[polyA])) == abs(GeometryUtil::polygonArea(polygons[polyB]))){
                if(finalMap.find(polyB)!=finalMap.end()){
                    //B hasnt been processed, so we points polyA to polyB and add polyA at remove list
                    finalMap[polyA] = polyB;
                    removeList.push_back(polyA);
                } // Otherwise the B is already been processed, so only points back to self.
            }
        }
        if(finalMap.find(polyA)==finalMap.end()){
            finalMap[polyA] = polyA;
        }
    }

    return finalMap;
}

std::vector<double> split_by_token(const std::string& s, const std::string& token) {
    std::vector<double> out;
    if (token.empty()) { // avoid infinite loop
        out.push_back(std::stoi(s));
        return out;
    }

    std::size_t pos = 0;
    while (true) {
        std::size_t next = s.find(token, pos);
        if (next == std::string::npos) {
            out.push_back(std::stod(s.substr(pos)));
            break;
        }
        out.push_back(std::stod(s.substr(pos, next - pos)));
        pos = next + token.size();
    }
    return out;
}

std::vector<double> sanitizeCuts(std::vector<double> cuts,double boardLength){
    sort(cuts.begin(),cuts.end());
    vector<double> cutsSanitized;
    for(auto cut:cuts){
        if (cut<0){
            cerr << "WARNING: negative cut length, ignore.\n";
            continue;
        }
        if(cut >= boardLength){
            cerr << "WARNING: cut is at least longer than board length, nothing will be cut\n";
        }
        cutsSanitized.push_back(cut);
            
    }
    return cutsSanitized;
}

unsigned int computeTotalVertices(vector<LayerPoints> layerOfPoint){
    unsigned int totalVertice = 0;
    for (auto layer:layerOfPoint){
        totalVertice += layer.indexRange.second-layer.indexRange.first+1;
    }
    return totalVertice;
}


vector<uint32_t> GetRemoveNodeList(vector<LayerPoints> layerOfPoint, double cut,double boardLength){

    //Truncates the cut and length since we're operating on the dotted board.
        
    int globalOffset = int(boardLength)-int(cut);
    vector<uint32_t> removeNodes;
    for (auto layer:layerOfPoint){
        int endIndex = layer.Ymap[0].size()-globalOffset;
        endIndex = endIndex>=0 ? endIndex : 0;
        for (auto xarr:layer.Ymap){
            
            for (int i = xarr.size()-1 ;i>=endIndex;i--){
                removeNodes.push_back(xarr[i]);
            }
        }
    }

    return removeNodes;

}

pair <vector<uint32_t>,vector<bool>> MakeMap(vector<uint32_t>  removeNodes,uint32_t totalVertice){
    vector<uint32_t> Map(totalVertice,1);
    vector<bool> Mask(totalVertice, true);
    for (uint32_t rem:removeNodes){
        Map[rem] = UINT32_MAX;
        Mask[rem] = false;
    }

    uint32_t increment=0;
    for (uint32_t&s:Map){
        if(s==1) s=increment++;
    }

    return { Map,Mask };
}


// Replaces all occurrences of from[i] with to[i] in the file at `path`.
// - from/to must have the same size
// - from[i] must not be empty (to avoid infinite loops)
// Returns true if the file content changed.
// Throws std::runtime_error on I/O errors or invalid arguments.
inline bool SearchReplaceInFile(const std::string& path,
                               const std::vector<std::string>& from,
                               const std::vector<std::string>& to)
{
    if (from.size() != to.size()) {
        throw std::runtime_error("SearchReplaceInFile: 'from' and 'to' sizes differ.");
    }
    for (std::size_t i = 0; i < from.size(); ++i) {
        if (from[i].empty()) {
            throw std::runtime_error("SearchReplaceInFile: from[" + std::to_string(i) + "] is empty.");
        }
    }

    // Read whole file
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("SearchReplaceInFile: failed to open for reading: " + path);
    }
    std::ostringstream ss;
    ss << in.rdbuf();
    std::string content = ss.str();

    bool changed = false;

    // Apply replacements in order: from[i] -> to[i]
    for (std::size_t i = 0; i < from.size(); ++i) {
        const std::string& f = from[i];
        const std::string& t = to[i];

        std::size_t pos = 0;
        while ((pos = content.find(f, pos)) != std::string::npos) {
            content.replace(pos, f.size(), t);
            pos += t.size();
            changed = true;
        }
    }

    if (!changed) return false;

    // Write back (overwrite)
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) {
        throw std::runtime_error("SearchReplaceInFile: failed to open for writing: " + path);
    }
    out.write(content.data(), static_cast<std::streamsize>(content.size()));
    if (!out) {
        throw std::runtime_error("SearchReplaceInFile: failed while writing: " + path);
    }

    return true;
}

bool CopyFile(const std::string& source,
              const std::string& destination,
              bool overwrite = false){
    std::error_code ec;

    fs::copy_options options = overwrite
        ? fs::copy_options::overwrite_existing
        : fs::copy_options::none;

    bool result = fs::copy_file(source, destination, options, ec);

    if (ec) {
        std::cerr << "Copy failed: " << ec.message() << "\n";
        return false;
    }

    return result;
}

void addLayerPolyClique(Graph &graph, LayersResult & layers){
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



int main(int argc, char** argv) {
    install_cancel_exit_handlers(130); // For non-0 exit code 
    auto dir_path = executable_path_from_argv0(argv[0]).parent_path();
    const string Usagephrase = "--instances <file|dir> Loads one instance JSON file, or a directory of instance JSON files. repeat for multiple file/directory nested directory will be ingnored.\n"
               "--datasets <dir> directory of dataset repeat for multiple directory\n"
               "--outputdir <dir>[optional default ./result] \n"
               "--typeOriented <0|1> [optional, default false] \n"
               "--cliqueCovering <0|1> [optional, default false] \n"
               "--width <number> [optional, mandatory if dataset doenst contain it] ignored if specified in json \n"
               "--length <number>[optional, mandatory if dataset doenst contain it]  ignored if specified in json\n"
               "--format <graph|lp|all> [optional, graph by default]\n"
               "--cuts <number,...,number>[optional, generate cut graph for single instance, ignored for multiple instances]\n"
               "--apply_cuts <0|1> [Optional (default: false). Generates the cut graph using the cuts defined in each instance. If none are defined, no cuts are applied.]";
    if (argc < 2) {
        cerr
            << "Usage: " << argv[0]
            << Usagephrase;
        return 1;
    }

    
    bool type_oriented = false; //one of the arguments
    vector<string> inputPaths; //one of the arguments
    vector<string> datasetPaths; //one of the arguments
    vector<double> cuts; //one of the arguments
    string outputdir = ""; //one of the arguments
    bool cliqueCovering = false; //one of the arguments
    json dataset;
    double width = 0.0; //one of the arguments
    double length = 0.0; //one of the arguments
    

    bool type_oriented_param = false;
    bool cliqueCovering_param = false;
    bool length_set = false;
    bool width_set = false;
    bool cuts_set = false;
    bool enable_cut = false;
    bool enable_cut_set = false;
    bool singleInstace = false;
    
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
        } else if (key == "--typeOriented") {
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
        } else if (key == "--cuts") {
            cuts = split_by_token(argv[argi + 1],",");
            cuts_set = true;
            enable_cut = true;
        } else if (key == "--apply_cuts") {
            enable_cut = parse_bool(argv[argi + 1]);
            enable_cut_set = true;
        } else {
            std::cerr << "Unknown argument: " << key << "\n";
            std::cerr << "Usage:" << argv[0] << " \n" << Usagephrase;
            return 1;
        }
    }
    
    if (datasetPaths.empty()){
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

    if (!enable_cut_set){
        cerr << "Warning: apply_cuts parameter not set, defaulting to false.\n";
    }

    
    if(REMOVE_OLD_FOLDER){
        cout << "removing existing output directory\n";
        try{
            std::filesystem::remove_all(outputdir);
        }catch(exception & e){
            cerr << "WARNING: fail to remove old directory\n";
        }
        
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

    if (datasets.size() == 1){
        singleInstace = true;
    }

    if(datasets.size() == 1 && cuts_set){
        datasets[0]["cuts"] = sanitizeCuts(cuts,dataset[0]["length"]);
    }else if(enable_cut){
        for (auto &set : datasets){
            //initialize the length as last cut (longest)
            if (set.contains("cuts")){
                set["lengthOriginal"] = set["length"].get<double>();
                set["cuts"] = sanitizeCuts(set["cuts"].get<vector<double>>(),set["length"].get<double>());
                set["length"] = set["cuts"].back();
                if (set["cuts"].is_array() && !set["cuts"].empty()) {
                    set["cuts"].erase(std::prev(set["cuts"].end()));  // remove last elemen
                }
            }
        }
        
    }else{ // no cut, make cut array empty
        for (auto &set : datasets){
            
            if (set.contains("cuts")){
                set["lengthOriginal"] = set["length"].get<double>();
                set["cuts"] = {};
            }
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
        
        string outputDataset = outputdir + outputname;
        if(singleInstace){outputDataset = outputdir;}

        if (enable_cut && set.contains("cuts") && set["cuts"].size()>0){
            if(!set.contains("cuts") || set["cuts"].size()==0) {
                cout << "WARNING: dataset: "<< outputname << " doesnt have cuts field, skipped.\n";
                continue;
            }

            outputDataset += "\\cut_"+to_string_fixed(length,1)+"\\";
        }
        
        
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

        //In case of Dataset with duplicated polygons we merge them by mapping one of them into another.
        std::cout << "Duplicate map processing..\n";
        //std::cout << dataset.dump(1) << "\n";
        unordered_map<string, string> DuplicatePolyMap = MapDuplicateDataset(dataset);

        //for (const auto& [key, value] : DuplicatePolyMap) {   // structured bindings (C++17)
        //    std::cout << key << " -> " << value << "\n";
        //}
        
        if (set.contains("quantity")){
            //Override the quantity field of dataset to 0
            for (auto& [key, val] : dataset.items()){
                dataset[key]["QUANTITY"] = 0;
            }

            if(set["quantity"].is_number_unsigned()) {
                for (auto& [key, val] : dataset.items()){
                    dataset[DuplicatePolyMap.at(key)]["QUANTITY"] = dataset[DuplicatePolyMap.at(key)]["QUANTITY"].get<unsigned int>() + set["quantity"].get<unsigned int>();
                }
            }
            else{
                for (auto& [key, val] : set["quantity"].items()){
                    dataset[DuplicatePolyMap.at(key)]["QUANTITY"] = dataset[DuplicatePolyMap.at(key)]["QUANTITY"].get<unsigned int>() + set["quantity"][key].get<unsigned int>();
                }
            }


            //Remove the polygon from dataset with 0 quantity
            vector<string> removeList;
            for (auto& [key, val] : dataset.items()){
                if(val["QUANTITY"].get<unsigned int>() == 0) removeList.push_back(key);
            }

            for (auto key:removeList){
                dataset.erase(key);
            }
            
        }

        //for (auto& [key, val] : dataset.items()){
        //    std::cout << key << " " << (dataset[key]["QUANTITY"]) << "\n";
        //}
        
        cout << "Generating NFP..\n";
        json polygons;
        try {
            //std::cout << "Processing NFP..\n";
            polygons = NFPTool::processNFP(dataset, length, width);
        } catch (const exception& e) {
            cerr << "Error during NFP processing: " << e.what() << "\n";
            return 1;
        }
        json polygonsWrite = polygons;
        //std::cout << "NFP Generated:" << polygons.dump(1) << "\n";
        // write the dataset with NFP into JSON file
        

        nlohmann::json board = polygons.at("rect");
        polygons.erase("rect"); 
        size_t num_polygon = polygons.size();

        //quantity can either be a single unsigned int for all polygons or a dictionary of per-polygon quantities

        int total_polygon = 0;
        double total_area = 0.0;

        for (auto& [key, val] : polygons.items()){
            total_polygon += polygons[key]["QUANTITY"].get<unsigned int>();
            total_area += polygons[key]["QUANTITY"].get<double>() * std::abs(GeometryUtil::polygonArea(GeometryUtil::parseVertices(polygons[key]["VERTICES"])));
        }

        NfpPosTableT NfpPosTable;
        NfpIndexTableT NfpIndexTable;

        std::cout << "Computing Innerfit Polygon BBox \n";
        for (auto& [key, poly] : polygons.items()){

            json& pivot = poly.at("VERTICES").at(0);
            json& innerfit = poly.at("innerfit");

            //pivot the polygon
            for (auto& v : polygons.at(key).at("VERTICES")) {
                v["x"] = v.at("x").get<double>() - pivot.at("x").get<double>();
                v["y"] = v.at("y").get<double>() - pivot.at("y").get<double>();
            } 

            //Compute innerfit BBox     
            //std::cout << "Computing BBox for innerfit\n" << key << " \n";       
            GeometryUtil::BBox Innerbbox = computeBoundingBox(innerfit);

            poly["innerfit_BoundingBox"]["xMin"] = floor(Innerbbox.xMin / gx);
            poly["innerfit_BoundingBox"]["xMax"] = ceil(Innerbbox.xMax / gx);
            poly["innerfit_BoundingBox"]["yMin"] = floor(Innerbbox.yMin / gy);
            poly["innerfit_BoundingBox"]["yMax"] = ceil(Innerbbox.yMax / gy);

            //Generate NFP Pattern             
            //std::cout << "Genering NFP pattern for " << key << " \n";
            buildNfpPatternsForPoly(key,poly,pivot,gx,gy,NfpIndexTable,NfpPosTable);
            
        }
        //generate Layer points (vertex)
        //std::cout << "Generating layer vertices\n";
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

        std::cout << "Making NFP Graph\n";
        Graph graph(static_cast<int>(layers.valIndex)); // valIndex is 1 + max vertex id
        for (const auto& [layerA, polyA] : layers.layerPoly) {
            const auto& pointsA = layers.layerOfPoint[layerA].innerFitPoints;

            for (const auto& [layerB, polyB] : layers.layerPoly) {
                if (layerA == layerB){
                    if(!type_oriented){
                        continue;
                    }
                }

                auto& nfpPatternB = NfpIndexTable.at(polyA).at(polyB);
                if (nfpPatternB.empty()) continue;

                // find innerMatrixB
                auto& innerMatrixB = layers.layerMatrix[layerB];

                // find BBoxB
                const auto& bboxBJson = polygons.at(polyB).at("innerfit_BoundingBox");
                GeometryUtil::BBox BBoxB{
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

        

        uint64_t NFPEdges = graph.getNumEdges();
        uint32_t NumberVertices = graph.getNumVertices();
        Graph graphCopy;

        // Optionally compute clique covering
        EdgeCliqueCover max1Cover;
        EdgeCliqueCover max1MinEKCover;
        
        if (cliqueCovering){
            //add layer regardless the difference is earlier and later.
            if(!type_oriented){
                //graphCopy useful if needed to compute NFP edges
                if (CLIQUE_COVERING_ADD_LAYER_CLIQUE){
                    cout << "Add layer enabled: add layer before clique covering\n";
                    graphCopy.loadGraph(graph); //make a copybackup before adding new edges;
                    //Add clique edges for each layer if not type oriented. And flag is set true
                    addLayerPolyClique(graph, layers);
                    //cout << "Graph edges before cover"<< graph.getNumEdges() << "\n";
                    cout << "Starting clique covering: Max1-EK\n";
                    max1Cover = maximum1Heuristic(graph);
                    //cout << "Start "
                    max1MinEKCover = expandKouHeuristic(graph, max1Cover);
                    //cout << "Graph edges after cover"<< graph.getNumEdges() << "\n";
                }else{
                    cout << "Add layer disabled: add layer after clique covering\n";
                    cout << "Starting clique covering: Max1-EK\n";
                    max1Cover = maximum1Heuristic(graph);
                    //cout << "Start "
                    max1MinEKCover = expandKouHeuristic(graph, max1Cover);
                    addLayerPolyClique(graph, layers);
                }
                
            }
            else{
                cout << "Starting clique covering: Max1-EK\n";
                max1Cover = maximum1Heuristic(graph);
                //cout << "Start "
                max1MinEKCover = expandKouHeuristic(graph, max1Cover);
            }

            
            
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

        cout << " graph edges from Graph: " << graph.getNumEdges() << "\n";

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




        //IO:print result to disk
  
        //TODO: write metadata on top of graph
        // Output graph to file
        cout << "saving result..\n";
        if(enable_cut){
            if(singleInstace){
                outputDataset = outputdir + "\\cut_" + to_string_fixed(length,1);
            }
            else{
                outputDataset = outputdir + "/" + outputname + "\\cut_" + to_string_fixed(length,1);
            }
            
        } 
        fs::create_directories(outputDataset);

        if(cliqueCovering){
            //const string cliqueOutputdir = outputDataset + "/BLAZEWICZ1_ECC_z.txt";
            const string cliqueOutputdir = outputDataset + "/" + outputname+"_ECC_"+to_string_fixed(length,1)+".txt";
            auto cliqueCount = max1MinEKCover.countRowMap();
            unsigned TotalCliqueCount = cliqueCount;
            if(OUTPUT_ADD_LAYER_CLIQUE){
                TotalCliqueCount += total_polygon;
            }
            //# number of cliques in the edge-clique cover (# of lines below): |C| = …
            std::ofstream out(cliqueOutputdir, std::ios::out | std::ios::trunc);
            
            out <<"# G_z number of vertices: |V| = " << NumberVertices << endl;
            out <<"# G_z number of edges: |E| = " << NFPEdges+cliqueEdges << endl;
            out <<"# number of cliques in the edge-clique cover (# of lines below): |C| = " << TotalCliqueCount << endl;
            out << fixed << setprecision(1) <<"# strip length z = " << length << endl;
            out <<"# number of polygons: |P| = " << total_polygon << endl;
            //out <<"# number of cliques in the edge-clique cover (# of lines below): |C| =" << cliqueCount+ layerCliques.size() << endl;
            //Piece oriented + clique covering add layer clique.
            if(!type_oriented && OUTPUT_ADD_LAYER_CLIQUE){
                cout << "adding layer cliques..\n";                 
                cout << "clique count: " << cliqueCount << endl;
                
                for (const auto& [layer, poly] : layers.layerPoly) {
                    const uint64_t start = layers.layerOfPoint[layer].indexRange.first;
                    const uint64_t end = layers.layerOfPoint[layer].indexRange.second;
                    //cout << "Printing layer clique\n";
                    for (uint64_t i = start; i <= end;i++) {
                        out << i << " ";
                    }
                    out << "\n";
                }   
            }
            out.close();
            max1MinEKCover.writeIntoFile(cliqueOutputdir);
        }
        else{
            const string graphOutputPath = outputDataset + "/" + outputname+"_graph_"+to_string_fixed(length,1)+".csv";
            graph.writeEdgesToFile(graphOutputPath);
        }
        
        // Write pointsCoordinate
        const string pointsOutputPath = outputDataset + "/" + outputname+"_graph2Strip_"+to_string_fixed(length,1)+".txt";
        ofstream pointsOut(pointsOutputPath, ios::out | ios::trunc);
        if (!pointsOut) {
            throw runtime_error("Failed to open file for writing: " + pointsOutputPath);
        }
        pointsOut << "# ID polygon, strip coord h (height), strip coord w (width), ID of vertex\n";
        for (const auto& lp : layers.layerOfPoint) {
            for (const auto& p : lp.innerFitPoints) {
                pointsOut << lp.layer << "," << p.x << "," << p.y << "," << p.id << "\n";
            }
        }
        
        
        std::cout << "Writting NFP into JSON\n";
        const string json_output = outputDataset + "/" + outputname+"_polygons_"+to_string_fixed(length,1)+".json"; 
        writeNfpInfpJson(json_output, polygonsWrite);
        //search and substitue polygon string after write.
        vector<string> from;
        vector<string> to;
        for (const auto& lp : layers.layerOfPoint) {
           from.push_back("\""+lp.polygon+"\"");
           to.push_back("\""+to_string(lp.layer)+"\"");
        }
        SearchReplaceInFile(json_output, from, to);

       
        ////Write metadata
        if(OUT_OLD_METADATA){
            const string metadataOutputPath = outputDataset + "/metadata.csv";
            ofstream metadataOut(metadataOutputPath, ios::out | ios::trunc);
            if (!metadataOut) {
                throw runtime_error("Failed to open file for writing: " + metadataOutputPath);
            }

            metadataOut <<"Name :\t" <<outputname << "\n";
            metadataOut <<"Total polygons :\t" << total_polygon << "\n";
            metadataOut <<"Type of polygons :\t" << num_polygon << "\n";
            metadataOut <<"Board Width :\t" << width << "\n";
            metadataOut <<"Board Length :\t" << int(length) << "\n";
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
        


        //TODO: cut and write for the rest.
        //if is clique covering, cut the already computed cliques and graph for metadata. DO NOT RECOMPUTE CLIQUE.
        
        //vector<uint32_t> cutMap;
        //vector<bool>Mask;
        for (auto cut:set["cuts"]){

            if(singleInstace){
                outputDataset = outputdir + "\\cut_"+to_string(cut)+"\\";
            }
            else{
                outputDataset = outputdir + "/" + outputname + "\\cut_"+to_string(cut)+"\\";
            }
            
            fs::create_directories(outputDataset);
            vector<uint32_t> removeList = GetRemoveNodeList(layers.layerOfPoint,cut,length);
            uint32_t TotalVertices = computeTotalVertices(layers.layerOfPoint);
            auto [cutMap,Mask] = MakeMap(removeList,TotalVertices);
            //prepare layer clique for the cut graph, which is the same as the original graph but with some vertices removed using the map. 
            //In case of one layer completely removed, also decrease the number of polygons by 1. This is for the metadata of the cut graph, which is different from the original graph. 
            //The cut graph will have less vertices and edges, and possibly less polygons if some layers are completely removed. 
            vector<vector<uint32_t>> layerCliques;
            vector<PointCoord> newIndex;
            uint64_t cliqueEdges = 0;

            for (auto layer:layers.layerOfPoint){
                
                vector<uint32_t> layerclique;
                
                for (const auto& p : layer.innerFitPoints) {
                    if(Mask[p.id]){
                        newIndex.push_back(PointCoord{layer.layer,p.x,p.y,cutMap[p.id]});
                        layerclique.push_back(cutMap[p.id]);
                    }
                }
                if(layerclique.size() > 1){
                    layerCliques.push_back(layerclique);
                    auto n = layerclique.size();
                    cliqueEdges += n*(n-1)/2;
                }
                
            }
            
            //copy json
            CopyFile(json_output, outputDataset + outputname+"_polygons_"+to_string(cut)+".json", true);

            //add pointcoord
            const string pointsOutputPath = outputDataset + "/"+outputname+"_graph2Strip_"+to_string(cut)+".txt";
            ofstream pointsOut(pointsOutputPath, ios::out | ios::trunc);
            if (!pointsOut) {
                throw runtime_error("Failed to open file for writing: " + pointsOutputPath);
            }
            pointsOut << "# ID polygon, strip coord h (height), strip coord w (width), ID of vertex\n";
            for (const auto& lp : newIndex) {
                pointsOut << lp.layer << "," << lp.x << "," << lp.y << "," << lp.id << "\n";
            }
            cout << "removing nodes in cut \n";
            graph.removeNodes(removeList);
            //TODO: fix for the case of CLIQUE_COVERING_ADD_LAYER_CLIQUE inside cliquecovering case.
            uint32_t NumberEdges = graph.getNumEdges();
            uint32_t NumberVertices = graph.getNumVertices();
            //TODO: make fater computation by using the copy graph before adding layer clique, and only remove nodes on the copy graph for clique covering, and keep the original graph for edge output.
            //Since clique covering doesnt change the graph structure.
            density = (2.0*NumberEdges)/(NumberVertices*NumberVertices-1);
            cout << "Graph statistics:\n";
            cout << "  Dataset name:     " << outputname << "\n";
            cout << "  Number of vertices: " << NumberVertices << "\n";
            cout << "  Number of clique edges:   " << cliqueEdges << "\n";
            cout << "  Number of edges:    " << NumberEdges << "\n";
            cout << "  Density:            " << density << "\n";

            
            if (cliqueCovering){
                const string cliqueOutputdir = outputDataset + "/"+outputname+"_ECC_"+to_string(cut)+".txt";
                auto cliqueCount = max1MinEKCover.countRowMap(cutMap);
                unsigned TotalCliqueCount = cliqueCount;
                if(OUTPUT_ADD_LAYER_CLIQUE){
                    TotalCliqueCount += layerCliques.size();
                }
                cout << "  Number of cliques (clique covering): " << TotalCliqueCount << "\n";
                //TODO: add latest version of metadata format
                //# number of cliques in the edge-clique cover (# of lines below): |C| = …
                std::ofstream out(cliqueOutputdir, std::ios::out | std::ios::trunc);
                
                out <<"# G_z number of vertices: |V| = " << NumberVertices << endl;
                out <<"# G_z number of edges: |E| = " << NumberEdges << endl;
                out <<"# number of cliques in the edge-clique cover (# of lines below): |C| = " << TotalCliqueCount << endl;
                out <<"# number of polygons: |P| = " << layerCliques.size() << endl;
                out <<"# strip length z = " << cut << endl;

                //Piece oriented + clique covering add layer clique.
                if(!type_oriented && OUTPUT_ADD_LAYER_CLIQUE){                   
                    
                    for (size_t i = 0; i < layerCliques.size(); ++i){
                        const auto& layerClique = layerCliques[i];
                        if (layerCliques[i].size() <= 1) continue;
                        for (const auto& v : layerClique) {
                            out << v << " ";
                        }
                        
                        //if (i + 1 != layerCliques.size()) out << "\n";
                        out << "\n";
                    }  
                    
                }
                out.close();
                max1MinEKCover.writeIntoFile(cliqueOutputdir,cutMap, Mask);
                
            }else{

                const string graphOutputPath = outputDataset + "/"+outputname+"_graph_"+to_string(length)+".txt";
                 std::ofstream out(graphOutputPath, std::ios::out | std::ios::trunc);
                    out <<"# G_z number of vertices: |V| = " << NumberVertices << endl;
                    out <<"# G_z number of edges: |E| = " << NumberEdges << endl;
                    out <<"# number of polygons: |P| = " << layerCliques.size() << endl;
                    out <<"# strip length z = " << cut << endl;
                graph.writeEdgesToFile(graphOutputPath,cutMap);
            }
            

            if(OUT_OLD_METADATA){
                const string metadataOutputPath = outputDataset + "/metadata.csv";
                ofstream metadataOut(metadataOutputPath, ios::out | ios::trunc);
                if (!metadataOut) {
                    throw runtime_error("Failed to open file for writing: " + metadataOutputPath);
                }

                metadataOut <<"Name :\t" <<outputname << "\n";
                metadataOut <<"Total polygons :\t" << total_polygon << "\n";
                metadataOut <<"Type of polygons :\t" << num_polygon << "\n";
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

            cout << "Finished processing cut : " << cut << "\n";
        }


    }
    return 0;
}