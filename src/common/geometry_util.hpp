#pragma once
#include <vector>
#include <string>

namespace GeometryUtil {

struct Point {
    double x{};
    double y{};
};

struct BBox {
    double xMin, xMax, yMin, yMax;
};

using Polygon  = std::vector<Point>;
using Polygons = std::vector<Polygon>;
 

inline constexpr double TOL = 1e-9;

// Public helpers you use elsewhere:
double polygonArea(const Polygon& polygon);
int pointInPolygon(const Point& point, const Polygon& polygon);
BBox computeBoundingBox(const Polygon& poly);

struct Bounds {
    double x{};
    double y{};
    double width{};
    double height{};
};

Bounds getPolygonBounds(const Polygon& polygon);
bool isRectangle(const Polygon& polygon, double tolerance = TOL);
Polygon parseVertices(const json& j);

} // namespace GeometryUtil
