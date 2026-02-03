#include "geometry_util.hpp"
#include <cmath>
#include <algorithm>

namespace GeometryUtil {

static bool almostEqual(double a, double b, double tol = TOL) {
    return std::abs(a - b) < tol;
}

static bool onSegment(const Point& A, const Point& B, const Point& p) {
    // vertical line
    if (almostEqual(A.x, B.x) && almostEqual(p.x, A.x)) {
        if (!almostEqual(p.y, B.y) &&
            !almostEqual(p.y, A.y) &&
            p.y < std::max(B.y, A.y) &&
            p.y > std::min(B.y, A.y)) {
            return true;
        } else {
            return false;
        }
    }

    // horizontal line
    if (almostEqual(A.y, B.y) && almostEqual(p.y, A.y)) {
        if (!almostEqual(p.x, B.x) &&
            !almostEqual(p.x, A.x) &&
            p.x < std::max(B.x, A.x) &&
            p.x > std::min(B.x, A.x)) {
            return true;
        } else {
            return false;
        }
    }

    // range check
    if ((p.x < A.x && p.x < B.x) ||
        (p.x > A.x && p.x > B.x) ||
        (p.y < A.y && p.y < B.y) ||
        (p.y > A.y && p.y > B.y)) {
        return false;
    }

    // exclude endpoints
    if ((almostEqual(p.x, A.x) && almostEqual(p.y, A.y)) ||
        (almostEqual(p.x, B.x) && almostEqual(p.y, B.y))) {
        return false;
    }

    // collinearity via cross product
    double cross = (p.y - A.y) * (B.x - A.x) - (p.x - A.x) * (B.y - A.y);
    if (std::abs(cross) > TOL) {
        return false;
    }

    // projection via dot product
    double dxAB = B.x - A.x;
    double dyAB = B.y - A.y;
    double dxAP = p.x - A.x;
    double dyAP = p.y - A.y;

    double dot = dxAP * dxAB + dyAP * dyAB;
    if (dot < 0.0 || almostEqual(dot, 0.0)) {
        return false;
    }

    double len2 = dxAB * dxAB + dyAB * dyAB;
    if (dot > len2 || almostEqual(dot, len2)) {
        return false;
    }

    return true;
}

double polygonArea(const Polygon& p) {
    if (p.size() < 3) return 0.0;
    long double a = 0.0;
    for (size_t i = 0, j = p.size() - 1; i < p.size(); j = i++) {
        a += (long double)p[j].x * p[i].y - (long double)p[i].x * p[j].y;
    }
    return static_cast<double>(a * 0.5);
}

int pointInPolygon(const Point& point, const Polygon& polygon) {
    if (polygon.size() < 3) {
        // JS: return null
        return -1;
    }

    bool inside = false;

    // JS supports polygon.offsetx/offsety; we assume 0 here
    const double offsetx = 0.0;
    const double offsety = 0.0;

    const size_t n = polygon.size();
    for (size_t i = 0, j = n - 1; i < n; j = i++) {
        double xi = polygon[i].x + offsetx;
        double yi = polygon[i].y + offsety;
        double xj = polygon[j].x + offsetx;
        double yj = polygon[j].y + offsety;

        // exactly on a vertex 
        if (almostEqual(xi, point.x) && almostEqual(yi, point.y)) {
            return -1;
        }

        // exactly on a segment 
        if (onSegment(Point{xi, yi}, Point{xj, yj}, point)) {
            return -1;
        }

        // ignore degenerate edges
        if (almostEqual(xi, xj) && almostEqual(yi, yj)) {
            continue;
        }

        // even-odd rule
        bool intersect = ((yi > point.y) != (yj > point.y)) &&
                         (point.x < (xj - xi) * (point.y - yi) / (yj - yi) + xi);
        if (intersect) inside = !inside;
    }

    return inside ? 1 : 0;
}

bool isRectangle(const Polygon& p, double eps) {
    if (p.size() != 4) return false;
    double minx = p[0].x, maxx = p[0].x, miny = p[0].y, maxy = p[0].y;
    for (auto& pt : p) {
        minx = std::min(minx, pt.x); maxx = std::max(maxx, pt.x);
        miny = std::min(miny, pt.y); maxy = std::max(maxy, pt.y);
    }
    auto near = [&](double a, double b){ return std::abs(a-b) <= eps; };
    int cornerHits = 0;
    for (auto& pt : p) {
        bool onCorner =
            (near(pt.x, minx) || near(pt.x, maxx)) &&
            (near(pt.y, miny) || near(pt.y, maxy));
        if (onCorner) cornerHits++;
    }
    return cornerHits == 4;
}

BBox computeBoundingBox(const Polygon& poly){
    double vxMax = std::numeric_limits<double>::lowest();
    double vxMin = std::numeric_limits<double>::max();
    double vyMax = std::numeric_limits<double>::lowest();
    double vyMin = std::numeric_limits<double>::max();
    for (const auto &p: poly)
    {
        if (p.x > vxMax) vxMax = p.x;
        if (p.x < vxMin) vxMin = p.x;
        if (p.y > vyMax) vyMax = p.y;
        if (p.y < vyMin) vyMin = p.y;
    }
    return BBox{vxMin,vxMax,vyMin,vyMax};

}


static bool samePoint(const GeometryUtil::Point& a,
                      const GeometryUtil::Point& b,
                      double eps = 1e-9)
{
    return std::abs(a.x - b.x) <= eps && std::abs(a.y - b.y) <= eps;
}



} // namespace GeometryUtil
