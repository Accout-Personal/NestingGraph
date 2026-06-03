#include "nfp_lib.hpp"
#include <string>
#include <vector>
#include <iostream>
#include "nfp/NFPGenerator.hpp"

using json = nlohmann::json;
namespace NFPTool{

static json dumpVertices(const GeometryUtil::Polygon& p) {
    json arr = json::array();
    for (const auto& pt : p) {
        arr.push_back({ {"x", pt.x}, {"y", pt.y} });
    }
    return arr;
}

static bool hasVerticesObject(const json& obj) {
    return obj.is_object() && obj.contains("VERTICES") && obj["VERTICES"].is_array();
}


json processNFP(json &dataset,double height, double width){
    //height = lenght
    
    GeometryUtil::Polygon rect = {
        {0, 0},
        {height, 0},
        {height, width},
        {0, width}
    };

    // Collect keys of polygon entries to avoid iterating over "STRIP" or non-polygons later
    std::vector<std::string> polyKeys;
    polyKeys.reserve(dataset.size());

    for (auto it = dataset.begin(); it != dataset.end(); ++it) {
        const std::string key = it.key();
        const json& val = it.value();
        if (hasVerticesObject(val)) {
            polyKeys.push_back(key);
        }
    }

    // Replace original dataset with the copy that has polygon entries replaced with numbers.
    json datasetCopy;
    unsigned int total_polygon = 0;
    std::vector<std::string> polyKeysNew;
    polyKeysNew.reserve(dataset.size());
    for (const auto& key : polyKeys) {
        datasetCopy[std::to_string(total_polygon)] = dataset[key];
        polyKeysNew.push_back(std::to_string(total_polygon));
        total_polygon++;
    }

    dataset = datasetCopy; 
    
    // Main loop: for each fixed polygon
    for (const auto& fixedKey : polyKeysNew) {
        auto& fixedObj = dataset[fixedKey];

        GeometryUtil::Polygon fixedPoly = GeometryUtil::parseVertices(fixedObj["VERTICES"]);
        if (fixedPoly.size() < 3) {
            // Keep structure but empty results
            fixedObj["INNER-FIT"] = json::array();
            fixedObj["NO-FIT"] = json::array();
            continue;
        }

        // ---- Inner-fit for rectangle ----
        GeometryUtil::BBox polyBBox = GeometryUtil::computeBoundingBox(fixedPoly);

        //Consider first point of the polygon as pivot.
        double FixPolyxPivot = fixedObj["VERTICES"][0]["x"].get<double>();
        double FixPolyyPivot = fixedObj["VERTICES"][0]["y"].get<double>();

        fixedObj["INNER-FIT"] = json::array();
        //the pivot point must resides inside or on the boundary of the bounding box.
        // thus always: >= than xMin, <= xMax >=yMin <= yMax
        fixedObj["INNER-FIT"].push_back({{"x", FixPolyxPivot - polyBBox.xMin},{"y", FixPolyyPivot - polyBBox.yMin}}); //Bottom left 
        fixedObj["INNER-FIT"].push_back({{"x", height - polyBBox.xMax + FixPolyxPivot},{"y", FixPolyyPivot - polyBBox.yMin}}); //Bottom right
        fixedObj["INNER-FIT"].push_back({{"x", height- polyBBox.xMax + FixPolyxPivot},{"y", width - polyBBox.yMax + FixPolyyPivot}}); //Top right
        fixedObj["INNER-FIT"].push_back({{"x", FixPolyxPivot - polyBBox.xMin},{"y",  width - polyBBox.yMax + FixPolyyPivot}}); // Top left


        // ---- Outer NFPs against all polygons ----
        fixedObj["NO-FIT"] = json::array();

        for (const auto& rotKey : polyKeysNew) {
            auto& rotObj = dataset[rotKey];

            GeometryUtil::Polygon rotPoly = GeometryUtil::parseVertices(rotObj["VERTICES"]);
            if (rotPoly.size() < 3) {
                // Still push entry to match original style? Optional.
                json nfpEntry;
                nfpEntry["POLYGON"] = rotKey;
                nfpEntry["VERTICES"] = json::array();
                fixedObj["NO-FIT"].push_back(nfpEntry);
                continue;
            }
            std::vector<std::pair<double,double>> polygonPairsFixed;
            polygonPairsFixed.reserve(fixedObj["VERTICES"].size());
            std::vector<std::pair<double,double>> polygonPairsRot;
            polygonPairsRot.reserve(rotObj["VERTICES"].size());
            double RotPolyxPivot = rotObj["VERTICES"][0]["x"].get<double>();
            double RotPolyyPivot = rotObj["VERTICES"][0]["y"].get<double>();

            for (auto const &fv : fixedObj["VERTICES"])
            {
                polygonPairsFixed.emplace_back(fv.at("x").get<double>()-FixPolyxPivot,fv.at("y").get<double>()-FixPolyyPivot);

            }

            for (auto const &ov : rotObj["VERTICES"])
            {
                polygonPairsRot.emplace_back(ov.at("x").get<double>()-RotPolyxPivot,ov.at("y").get<double>()-RotPolyyPivot);

            }
            //auto nfpGen = NFPGeneratorCGAL();
            auto outNfps = nfp::processNFP(
                polygonPairsFixed, polygonPairsRot
            );

            json nfpEntry;
            nfpEntry["POLYGON"] = rotKey;
            
            nfpEntry["NFP_HOLES"] =  json::array();
            nfpEntry["NFP_SLITS"] =  json::array();
            nfpEntry["NFP_POINTS"] =  json::array();
            if (!outNfps.cycles.empty()) {
                
                for (auto cycle : outNfps.cycles){
                    json JoutCycles = json::array();
                    for (auto outV :cycle.points) JoutCycles.push_back({{"x",outV.x},{"y",outV.y}});

                    std::string entry;
                    if (cycle.role == nfp::NFPCycleRole::OuterBoundary){
                        nfpEntry["VERTICES"] = JoutCycles;
                    }
                    else if(cycle.role == nfp::NFPCycleRole::HoleBoundary){
                        nfpEntry["NFP_HOLES"].push_back(JoutCycles);
                    }

                }
                

                for (auto slit: outNfps.slits){
                    //std::cout << "slits detected: x1:" << slit.a.x << " y1:" << slit.a.y << " x2:" << slit.b.x << " y2:" << slit.b.y << std::endl;
                    json JoutSlits = json::array();
                    JoutSlits.push_back({ {"x",slit.a.x},{"y",slit.a.y} });
                    JoutSlits.push_back({ {"x",slit.b.x},{"y",slit.b.y} });
                    nfpEntry["NFP_SLITS"].push_back(JoutSlits);
                }

                json JoutPoint = json::array();
                for (auto point: outNfps.isolated_points){
                    JoutPoint.push_back({ {"x",point.x},{"y",point.y} });
                }
                
                if (! JoutPoint.empty()){
                    nfpEntry["NFP_POINTS"] = JoutPoint;
                }                
                
            } else {
                
                std::cerr << "ERROR: Empty NFP. Terminating\n";
                std::cerr << "Fixed Polygon:" << fixedKey << " RotPolygon: " << rotKey << std::endl; 
                for (auto SingleNfps: outNfps.cycles)
                {   
                    std::cerr <<"printing nfps..\n";
                    for (auto outV: SingleNfps.points) std::cerr << "x: "<< outV.x <<" y:" << outV.y <<std::endl;
                }
                
                //std::cerr << "datset: " << dataset << "\n";
                exit(1);
                nfpEntry["VERTICES"] = json::array();
            }

            fixedObj["NO-FIT"].push_back(nfpEntry);
        }
    }
    // Store rect in dataset (same top-level field as JS)
    dataset["STRIP"] = dumpVertices(rect);
    return dataset;
}
}


