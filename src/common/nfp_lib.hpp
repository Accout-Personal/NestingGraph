#pragma once
#include <nlohmann/json.hpp>
#include "common/geometry_util.hpp"

using json = nlohmann::json;
namespace NFPTool{

    json processNFP(json &dataset,double height, double width);
    GeometryUtil::Polygon parseVertices(const json& j);

}