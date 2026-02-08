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


json processNFP(json &dataset,double height, double width)
{
    //height = lenght
    
    GeometryUtil::Polygon rect = {
        {0, 0},
        {height, 0},
        {height, width},
        {0, width}
    };

    // Collect keys of polygon entries to avoid iterating over "rect" or non-polygons later
    std::vector<std::string> polyKeys;
    polyKeys.reserve(dataset.size());

    for (auto it = dataset.begin(); it != dataset.end(); ++it) {
        const std::string key = it.key();
        const json& val = it.value();
        if (hasVerticesObject(val)) {
            polyKeys.push_back(key);
        }
    }

    // Main loop: for each fixed polygon
    for (const auto& fixedKey : polyKeys) {
        auto& fixedObj = dataset[fixedKey];

        GeometryUtil::Polygon fixedPoly = GeometryUtil::parseVertices(fixedObj["VERTICES"]);
        if (fixedPoly.size() < 3) {
            // Keep structure but empty results
            fixedObj["innerfit"] = json::array();
            fixedObj["nfps"] = json::array();
            continue;
        }

        // ---- Inner-fit for rectangle ----
        //TODO:implement bounding box innerfit polygon
        GeometryUtil::BBox polyBBox = GeometryUtil::computeBoundingBox(fixedPoly);

        //Consider first point of the polygon as pivot.
        double FixPolyxPivot = fixedObj["VERTICES"][0]["x"].get<double>();
        double FixPolyyPivot = fixedObj["VERTICES"][0]["y"].get<double>();

        fixedObj["innerfit"] = json::array();
        //the pivot point must resides inside or on the boundary of the bounding box.
        // thus always: >= than xMin, <= xMax >=yMin <= yMax
        fixedObj["innerfit"].push_back({{"x", FixPolyxPivot - polyBBox.xMin},{"y", FixPolyyPivot - polyBBox.yMin}}); //Bottom left 
        fixedObj["innerfit"].push_back({{"x", height - polyBBox.xMax + FixPolyxPivot},{"y", FixPolyyPivot - polyBBox.yMin}}); //Bottom right
        fixedObj["innerfit"].push_back({{"x", height- polyBBox.xMax + FixPolyxPivot},{"y", width - polyBBox.yMax + FixPolyyPivot}}); //Top right
        fixedObj["innerfit"].push_back({{"x", FixPolyxPivot - polyBBox.xMin},{"y",  width - polyBBox.yMax + FixPolyyPivot}}); // Top left


        // ---- Outer NFPs against all polygons ----
        fixedObj["nfps"] = json::array();

        for (const auto& rotKey : polyKeys) {
            auto& rotObj = dataset[rotKey];

            GeometryUtil::Polygon rotPoly = GeometryUtil::parseVertices(rotObj["VERTICES"]);
            if (rotPoly.size() < 3) {
                // Still push entry to match original style? Optional.
                json nfpEntry;
                nfpEntry["POLYGON"] = rotKey;
                nfpEntry["VERTICES"] = json::array();
                fixedObj["nfps"].push_back(nfpEntry);
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

            if (!outNfps.empty() && !outNfps[0].empty()) {
                json Jout = json::array();
                for (auto outV :outNfps[0]) Jout.push_back({{"x",outV.first},{"y",outV.second}});
                nfpEntry["VERTICES"] = Jout;
                
            } else {
                
                std::cerr << "ERROR: Empty NFP. Terminating\n";
                std::cerr << "Fixed Polygon:" << fixedKey << " RotPolygon: " << rotKey << std::endl; 
                for (auto SingleNfps: outNfps)
                {   
                    std::cerr <<"printing nfps..\n";
                    for (auto outV: SingleNfps) std::cerr << "x: "<< outV.first <<" y:" << outV.second <<std::endl;
                }
                
                //std::cerr << "datset: " << dataset << "\n";
                exit(1);
                nfpEntry["VERTICES"] = json::array();
            }

            fixedObj["nfps"].push_back(nfpEntry);
        }
    }
    // Store rect in dataset (same top-level field as JS)
    dataset["rect"] = dumpVertices(rect);
    return dataset;
}
}


