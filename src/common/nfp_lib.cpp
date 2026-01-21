#include "nfp_lib.hpp"
#include <string>
#include <vector>
#include <iostream>

using json = nlohmann::json;
namespace NFPTool{
// ---------- JSON helpers ----------
GeometryUtil::Polygon parseVertices(const json& j) {
    GeometryUtil::Polygon poly;
    if (!j.is_array()) return poly;

    poly.reserve(j.size());
    for (const auto& pt : j) {
        if (!pt.is_object()) continue;
        GeometryUtil::Point p;
        p.x = pt.value("x", 0.0);
        p.y = pt.value("y", 0.0);
        poly.push_back(p);
    }
    return poly;
}

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

        GeometryUtil::Polygon fixedPoly = parseVertices(fixedObj["VERTICES"]);
        if (fixedPoly.size() < 3) {
            // Keep structure but empty results
            fixedObj["innerfit"] = json::array();
            fixedObj["nfps"] = json::array();
            continue;
        }

        // ---- Inner-fit for rectangle ----
        {
            auto inner = GeometryUtil::noFitPolygonLite(
                rect, fixedPoly,
                /*inside=*/true,
                /*searchEdges=*/false
            );

            if (!inner.empty() && !inner[0].empty()) {
                fixedObj["innerfit"] = dumpVertices(inner[0]);
            } else {
                fixedObj["innerfit"] = json::array();
            }
        }

        // ---- Outer NFPs against all polygons ----
        fixedObj["nfps"] = json::array();

        for (const auto& rotKey : polyKeys) {
            auto& rotObj = dataset[rotKey];

            GeometryUtil::Polygon rotPoly = parseVertices(rotObj["VERTICES"]);
            if (rotPoly.size() < 3) {
                // Still push entry to match original style? Optional.
                json nfpEntry;
                nfpEntry["POLYGON"] = rotKey;
                nfpEntry["VERTICES"] = json::array();
                fixedObj["nfps"].push_back(nfpEntry);
                continue;
            }

            auto out = GeometryUtil::noFitPolygonLite(
                fixedPoly, rotPoly,
                /*inside=*/false,
                /*searchEdges=*/false
            );

            json nfpEntry;
            nfpEntry["POLYGON"] = rotKey;

            if (!out.empty() && !out[0].empty()) {
                nfpEntry["VERTICES"] = dumpVertices(out[0]);
            } else {
                std::cerr << "ERROR: Empty NFP. Terminating\n";
                std::cerr << "datset: " << dataset << "\n";
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


