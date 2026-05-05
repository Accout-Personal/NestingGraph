#pragma once
#include <vector>
#include <string>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace GeometryUtil {
    inline constexpr double TOL = 1e-9;

    struct Point {
        double x{};
        double y{};
    
        bool operator==(const Point& other) const noexcept {
            return abs(x - other.x) < TOL && abs(y -other.y) < TOL;
        }
    };
    
    struct BBox {
        double xMin, xMax, yMin, yMax;
    };
    
    using Polygon  = std::vector<Point>;
    using Polygons = std::vector<Polygon>;
    
    
    
    
    // Public helpers you use elsewhere:
    double polygonArea(const Polygon& polygon);
    int pointInPolygon(const Point& point, const Polygon& polygon);
    bool pointInSegment(const Point& point, const Polygon& slit);
    BBox computeBoundingBox(const Polygon& poly);
    
    struct Bounds {
        double x{};
        double y{};
        double width{};
        double height{};
    };
    
    bool isRectangle(const Polygon& polygon, double tolerance = TOL);
    Polygon parseVertices(const json& j);
    bool samePolygonVertices(Polygon A, Polygon B,bool allowReverse = true);

} // namespace GeometryUtil
