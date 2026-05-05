#include "geometry_util.hpp"
#include <cmath>
#include <algorithm>
#include <nlohmann/json.hpp>
#include <iostream>
using json = nlohmann::json;
namespace GeometryUtil {

    
// ---------- JSON helpers ----------
 GeometryUtil::Polygon parseVertices(const json& j) {
    GeometryUtil::Polygon poly;
    if (!j.is_array()) return poly;

    poly.reserve(j.size());
    for (const auto& pt : j) {
        if (!pt.is_object()) continue;
        GeometryUtil::Point p;
        try{
            p.x = pt.at("x").get<double>();
            p.y = pt.at("y").get<double>();
            poly.push_back(p);
        }catch(const std::exception& e){
            std::cerr << "ERROR: parse vertices failed! (make sure the vertices \"x\" and \"y\" are valid numerical values)" << std::endl;
            exit(-1);
        }
    }
    return poly;
}

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

bool pointInSegment(const Point& point, const Polygon& slit) {

    if (slit.size() != 2) return false;
    // exactly on vertices
    if (almostEqual(slit[0].x, point.x) && almostEqual(slit[0].y, point.y)) {
        return true;
    }

    if (almostEqual(slit[1].x, point.x) && almostEqual(slit[1].y, point.y)) {
        return true;
    }

    // exactly on a segment 
    if (onSegment(Point{ slit[0].x, slit[0].y }, Point{ slit[1].x, slit[1].y }, point)) {
        return true;
    }

    return false;
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

static bool pointLess(const Point& a, const Point& b) {
    if (!almostEqual(a.x ,b.x)) return a.x < b.x;
    return a.y < b.y;
}

static void stripClosingPoint(Polygon& p) {
    if (p.size() >= 2 && samePoint(p.front(),p.back())) p.pop_back();
}

static void removeConsecutiveDuplicates(Polygon& p) {
    if (p.empty()) return;
    Polygon out;
    out.reserve(p.size());
    out.push_back(p[0]);
    for (size_t i = 1; i < p.size(); ++i) {
        if (!samePoint(p[i],out.back())) out.push_back(p[i]);
    }
    p.swap(out);
}

static Polygon rotateToSmallest(const Polygon& p) {
    if (p.empty()) return p;
    auto it = std::min_element(p.begin(), p.end(), pointLess);
    Polygon out;
    out.reserve(p.size());
    out.insert(out.end(), it, p.end());
    out.insert(out.end(), p.begin(), it);
    return out;
}

static Polygon canonical(Polygon p) {
    stripClosingPoint(p);
    removeConsecutiveDuplicates(p);
    return rotateToSmallest(p);
}

bool samePolygonVertices(Polygon A,
                         Polygon B,
                         bool allowReverse){
    auto A0 = canonical(std::move(A));
    auto B0 = canonical(std::move(B));

    if (A0.size() != B0.size()) return false;
    if (A0 == B0) return true;

    if (!allowReverse) return false;

    // Canonicalize reversed B
    std::reverse(B0.begin(), B0.end());
    B0 = rotateToSmallest(B0);
    return A0 == B0;
}


} // namespace GeometryUtil
