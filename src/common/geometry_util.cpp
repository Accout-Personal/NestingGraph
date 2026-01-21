#include "geometry_util.hpp"
#include <cmath>
#include <algorithm>

// Clipper2
#include "clipper2/clipper.h"


using namespace Clipper2Lib;

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

static Point64 toP64(const Point& p, double scale) {
    return Point64(
        static_cast<int64_t>(std::llround(p.x * scale)),
        static_cast<int64_t>(std::llround(p.y * scale))
    );
}

static Point toP(const Point64& p, double scale) {
    return Point{ static_cast<double>(p.x) / scale,
                  static_cast<double>(p.y) / scale };
}

static Path64 toPath64(const Polygon& poly, double scale) {
    Path64 out;
    out.reserve(poly.size());
    for (const auto& pt : poly) out.push_back(toP64(pt, scale));
    return out;
}

static Polygon toPolygon(const Path64& path, double scale) {
    Polygon out;
    out.reserve(path.size());
    for (const auto& pt : path) out.push_back(toP(pt, scale));
    return out;
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

static bool samePoint(const GeometryUtil::Point& a,
                      const GeometryUtil::Point& b,
                      double eps = 1e-9)
{
    return std::abs(a.x - b.x) <= eps && std::abs(a.y - b.y) <= eps;
}

static GeometryUtil::Polygon minimalClean(GeometryUtil::Polygon p, double eps = 1e-9)
{
    if (p.size() < 2) return p;

    // remove closing duplicate if present
    if (samePoint(p.front(), p.back(), eps)) {
        p.pop_back();
    }

    // remove consecutive duplicates
    GeometryUtil::Polygon out;
    out.reserve(p.size());
    for (const auto& pt : p) {
        if (out.empty() || !samePoint(out.back(), pt, eps)) {
            out.push_back(pt);
        }
    }

    // if still closed after compaction
    if (out.size() > 1 && samePoint(out.front(), out.back(), eps)) {
        out.pop_back();
    }

    return out;
}

static void CollectOuterContours(const Clipper2Lib::PolyPath64& node,
                                 Clipper2Lib::Paths64& outers)
{
    for (const auto& childPtr : node)  // childPtr is usually std::unique_ptr<PolyPath64>
    {
        const auto& child = *childPtr;
        if (!child.IsHole() && !child.Polygon().empty())
            outers.push_back(child.Polygon());

        CollectOuterContours(child, outers);
    }
}


Paths64 GetOuterContours(const PolyTree64& tree)
{
    Paths64 outers;
    CollectOuterContours(tree, outers);
    return outers;
}


Polygons minkowskiDifference(const Polygon& A, const Polygon& B, double scale)
{
    if (A.size() < 3 || B.size() < 3) return {};

    Path64 Ac = toPath64(A, scale);
    Path64 Bc = toPath64(B, scale);

    for (auto& pt : Bc) { pt.x = -pt.x; pt.y = -pt.y; }

    Paths64 solution = MinkowskiSum(Ac, Bc, true);
    //Clipper64 c;
    //c.PreserveCollinear(true);
    //c.AddSubject(solution);
//
    //PolyTree64 tree;
    //c.Execute(ClipType::Union, FillRule::NonZero, tree);
    //Paths64 outers;
    //CollectOuterContours(tree, outers); 
    //TODO:extract countour from the tree.

    if (solution.empty()) return {};

    // JS-like: pick smallest signed area
    int bestIdx = 0;
    double bestArea = 0.0;
    bool bestSet = false;

    for (int i = 0; i < (int)solution.size(); ++i) {
        Polygon cand = toPolygon(solution[i], scale);
        double a = polygonArea(cand);

        if (!bestSet || a < bestArea) {
            bestArea = a;
            bestIdx = i;
            bestSet = true;
        }
    }

    Polygon nfp = toPolygon(solution[bestIdx], scale);

    // translate by B[0] like JS
    const double dx = B[0].x;
    const double dy = B[0].y;
    for (auto& pt : nfp) {
        pt.x += dx;
        pt.y += dy;
    }

    // ONLY minimal clean
    nfp = minimalClean(std::move(nfp),0.5/scale);

    return nfp.size() >= 3 ? Polygons{ nfp } : Polygons{};
}

// "as original" logic shape, assuming searchEdges is typically false
Polygons noFitPolygonLite(const Polygon& A,
                          const Polygon& B,
                          bool inside,
                          bool searchEdges)
{
    Polygons nfp;

    if (inside) {
        if (isRectangle(A, 0.001)) {
            // your working rectangle inner-fit path
            nfp = minkowskiDifference(A, B);
        } else {
            // leave as original intent placeholder
            nfp = {};
        }

        // ensure winding consistency like JS
        for (auto& p : nfp) {
            if (polygonArea(p) > 0) {
                std::reverse(p.begin(), p.end());
            }
        }
        return nfp;
    }

    // outside
    if (searchEdges) {
        // original intent placeholder
        nfp = {};
    } else {
        nfp = minkowskiDifference(A, B);
    }

    if (nfp.empty())
    {
        std::cerr << "NFP ERROR: Minkowswki results empty set.\n";
        return {};
    }

    if(false)
    {

    
    // sanity check-ish
    for (size_t i = 0; i < nfp.size(); ++i) {
        
        if (!searchEdges || i == 0) {
            if (std::abs(polygonArea(nfp[i])) < std::abs(polygonArea(A))) {
                
                std::cerr << "NFP ERROR: NFP Area smaller than original Polygon.\n"
                <<"Area NFP: " << std::abs(polygonArea(nfp[i])) << " "
                <<"Area Polygon A: " << std::abs(polygonArea(A)) << "\n"
                << "nfp size: " << nfp.size() << "\n";

                std::cerr << "Vertices NFP\n";
                for (auto v : nfp[i])
                {
                    std::cerr << v.x << " " << v.y << "\n";
                }

                std::cerr << "Vertices Polygon A\n";
                for (auto v : A)
                {
                    std::cerr << v.x << " " << v.y << "\n";
                }
                return {};
            }
        }
    }
    }

    if (nfp.empty()) return {};

    // winding + hole relationship
    for (size_t i = 0; i < nfp.size(); ++i) {
        if (polygonArea(nfp[i]) > 0) {
            std::reverse(nfp[i].begin(), nfp[i].end());
        }

        if (i > 0) {
            if (pointInPolygon(nfp[i][0], nfp[0]) == 1) {
                if (polygonArea(nfp[i]) < 0) {
                    std::reverse(nfp[i].begin(), nfp[i].end());
                }
            }
        }
    }

    return nfp;
}

} // namespace GeometryUtil
