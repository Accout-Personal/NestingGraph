#pragma once
#include <vector>

namespace GeometryUtil {

struct Point {
    double x{};
    double y{};
};

using Polygon  = std::vector<Point>;
using Polygons = std::vector<Polygon>;

inline constexpr double TOL = 1e-9;

// “Public” helpers you use elsewhere:
double polygonArea(const Polygon& polygon);
int pointInPolygon(const Point& point, const Polygon& polygon);

struct Bounds {
    double x{};
    double y{};
    double width{};
    double height{};
};

Bounds getPolygonBounds(const Polygon& polygon);
bool isRectangle(const Polygon& polygon, double tolerance = TOL);

// (your existing declarations)
Polygons minkowskiDifference(const Polygon& A, const Polygon& B, double scale = 1e7);
Polygons noFitPolygonLite(const Polygon& A, const Polygon& B, bool inside, bool searchEdges);

} // namespace GeometryUtil
