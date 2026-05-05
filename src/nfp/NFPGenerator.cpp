#ifndef NFP_GENERATOR_H
#define NFP_GENERATOR_H


#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Nef_polyhedron_2.h>
#include <CGAL/minkowski_sum_2.h>
#include <CGAL/partition_2.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/Aff_transformation_2.h>
#include <CGAL/Extended_homogeneous.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#include <CGAL/number_utils.h>
#include "NFPGenerator.hpp"
#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include <set>

namespace nfp {
// Kernel typedefs
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2                                   Point_2;
typedef CGAL::Polygon_2<Kernel>                           Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel>                Polygon_with_holes_2;

// Partition traits - this ensures output polygons match our Polygon_2 type
typedef CGAL::Partition_traits_2<Kernel>                  Partition_traits;
typedef Partition_traits::Polygon_2                       Partition_polygon;  // Uses std::list internally

// Nef kernel typedefs - MUST use Gmpz (integer type), not Exact_rational (field type)
// Extended_homogeneous requires an integer ring for polynomial GCD operations
typedef CGAL::Extended_homogeneous<CGAL::Gmpz>            Nef_kernel;
typedef CGAL::Nef_polyhedron_2<Nef_kernel>                Nef_polyhedron;

// Standard kernel for input points to Nef constructor
typedef CGAL::Homogeneous<CGAL::Gmpz>                     Standard_kernel;
typedef Standard_kernel::Point_2                          Standard_point;

typedef int ShapeID; // Simple alias for shape identifiers


    // Helper: Decompose a concave polygon into convex pieces
     std::vector<Polygon_2> decompose(const Polygon_2& p) {
        // Use partition traits polygons for output (they use std::list internally)
        if (p.is_convex()) {
            return { p };
        }

        std::vector<Partition_polygon> partition_parts;

        CGAL::approx_convex_partition_2(
            p.vertices_begin(), 
            p.vertices_end(), 
            std::back_inserter(partition_parts),
            Partition_traits());

        
        //OPTIMAL DECOMPOSITION IS BUGGED! DO NOT USE.
        //CGAL::optimal_convex_partition_2(
        //    p.vertices_begin(), 
        //    p.vertices_end(), 
        //    std::back_inserter(partition_parts),
        //    Partition_traits());  // Pass the traits!

        
    
        // Convert partition polygons (list-based) to our Polygon_2 (vector-based)
        std::vector<Polygon_2> parts;
        parts.reserve(partition_parts.size());
        
        for (const auto& pp : partition_parts) {
            Polygon_2 converted;
            for (auto vit = pp.vertices_begin(); vit != pp.vertices_end(); ++vit) {
                converted.push_back(*vit);
            }
            parts.push_back(std::move(converted));
        }
        
        return parts;
    }

    // Helper: Reflect polygon through origin (-B)
     Polygon_2 reflect_polygon(const Polygon_2& p) {
        Polygon_2 p_reflected;
        for (auto vit = p.vertices_begin(); vit != p.vertices_end(); ++vit) {
            p_reflected.push_back(Point_2(-vit->x(), -vit->y()));
        }
        
        // Reflection reverses orientation, fix it to CCW
        if (p_reflected.is_clockwise_oriented()) {
            p_reflected.reverse_orientation();
        }
        return p_reflected;
    }

    // Helper: Convert EPECK point to Homogeneous<Gmpz> point
     Standard_point convert_to_standard_point(const Point_2& pt) {
        // Get the exact coordinates - EPECK::FT is a lazy wrapper around mpq_class
        // We use CGAL's exact() to get the underlying exact type, then convert via strings
        // to avoid direct mpq_class <-> Gmpq conversion issues
        
        auto x_exact = pt.x().exact();
        auto y_exact = pt.y().exact();
        
        // Convert through string representation to avoid gmpxx vs CGAL::Gmp type issues
        std::ostringstream oss_x, oss_y;
        oss_x << x_exact;
        oss_y << y_exact;
        
        CGAL::Gmpq x_q(oss_x.str());
        CGAL::Gmpq y_q(oss_y.str());
        
        // Extract numerators and denominators
        CGAL::Gmpz x_num = x_q.numerator();
        CGAL::Gmpz x_den = x_q.denominator();
        CGAL::Gmpz y_num = y_q.numerator();
        CGAL::Gmpz y_den = y_q.denominator();
        
        // Homogeneous point (hx, hy, hw) represents Cartesian (hx/hw, hy/hw)
        // To represent (x_num/x_den, y_num/y_den):
        // hx = x_num * y_den
        // hy = y_num * x_den  
        // hw = x_den * y_den
        CGAL::Gmpz hx = x_num * y_den;
        CGAL::Gmpz hy = y_num * x_den;
        CGAL::Gmpz hw = x_den * y_den;
        
        return Standard_point(hx, hy, hw);
     }

     static void sanitize_polygon(Polygon_2& p) {
         if (p.size() < 2) return;

         // Remove last if equals first
         while (p.size() >= 2 && *p.vertices_begin() == *std::prev(p.vertices_end())) {
             p.erase(std::prev(p.vertices_end()));
         }

         // Remove consecutive duplicates
         auto it = p.vertices_begin();
         while (it != p.vertices_end() && p.size() > 1) {
             auto next = std::next(it);
             if (next == p.vertices_end()) next = p.vertices_begin();
             if (*it == *next && next != p.vertices_begin()) {
                 it = p.erase(next);
             }
             else {
                 ++it;
             }
         }
     }


    // THE CORE PIPELINE
     Nef_polyhedron get_nfp(ShapeID idA, Polygon_2& polyA,
        ShapeID idB, Polygon_2& polyB) {

        
        std::vector<Polygon_2> partsA = decompose(polyA);
        
        Polygon_2 polyB_reflected = reflect_polygon(polyB);

        std::vector<Polygon_2> partsB = decompose(polyB_reflected);

        //Handling for convex cases.
        if (partsA.size() == 1 && partsB.size() == 1) {
            // Simple case: both convex, just one Minkowski sum
            Polygon_2 sum = CGAL::minkowski_sum_2(partsA[0], partsB[0]).outer_boundary();

            // Return boundary directly as closed polygon (no EXCLUDED trick needed)
            std::vector<Standard_point> pts;
            for (auto v = sum.vertices_begin(); v != sum.vertices_end(); ++v)
                pts.push_back(convert_to_standard_point(*v));

            Nef_polyhedron result(pts.begin(), pts.end(), Nef_polyhedron::INCLUDED);
            return result;
        }

        // Store individual Minkowski sums
        Nef_polyhedron U_Open(Nef_polyhedron::EMPTY);

        for (const auto& subA : partsA) {
            for (const auto& subB : partsB) {

                auto sum = CGAL::minkowski_sum_2(subA, subB).outer_boundary();
                std::vector<Standard_point> pts;
                for (auto v = sum.vertices_begin(); v != sum.vertices_end(); ++v)
                    pts.push_back(convert_to_standard_point(*v));

                Nef_polyhedron nef_sub(pts.begin(), pts.end(), Nef_polyhedron::EXCLUDED);
                U_Open = U_Open.join(nef_sub);
            }
        }

        return U_Open;
     }
     
    struct ExactNefPoint {
        CGAL::Gmpz hx;
        CGAL::Gmpz hy;
        CGAL::Gmpz hw;
    };

    struct NFPPointCell {
        bool marked;
        ExactNefPoint p;
    };

    struct NFPEdgeCell {
        bool edge_marked;

        bool left_face_marked;
        bool right_face_marked;

        bool left_face_unbounded;
        bool right_face_unbounded;

        ExactNefPoint a;
        ExactNefPoint b;
    };

    struct NFPFaceCycle {
        bool face_marked;
        bool is_hole;
        bool unbounded;

        std::vector<ExactNefPoint> points;
    };

    struct NFPTopology {
        std::vector<NFPFaceCycle> face_cycles;
        std::vector<NFPFaceCycle> hole_cycles;
        std::vector<NFPEdgeCell> edges;
        std::vector<NFPPointCell> points;
    };

    struct DirectedBoundaryEdge {
        ExactNefPoint a;
        ExactNefPoint b;

        std::string a_key;
        std::string b_key;

        bool used = false;
    };

    struct NFPSlitCandidate {
        ExactNefPoint a;
        ExactNefPoint b;
    };

    struct NFPPointCandidate {
        ExactNefPoint p;
    };

    struct BuiltNFP {
        std::vector<std::vector<ExactNefPoint>> boundary_cycles;
        std::vector<std::vector<ExactNefPoint>> feasible_pocket_cycles;
        std::vector<NFPSlitCandidate> slit_candidates;
        std::vector<NFPPointCandidate> point_candidates;
    };

    static CGAL::Gmpz abs_gmpz(CGAL::Gmpz v) {
        return v < CGAL::Gmpz(0) ? -v : v;
    }

    static CGAL::Gmpz gcd_gmpz(CGAL::Gmpz a, CGAL::Gmpz b) {
        a = abs_gmpz(a);
        b = abs_gmpz(b);

        while (b != CGAL::Gmpz(0)) {
            CGAL::Gmpz r = a % b;
            a = b;
            b = r;
        }

        return a;
    }

    static CGAL::Gmpz gcd3_gmpz(
        const CGAL::Gmpz& a,
        const CGAL::Gmpz& b,
        const CGAL::Gmpz& c
    ) {
        return gcd_gmpz(gcd_gmpz(a, b), c);
    }

    static ExactNefPoint normalize_exact_point(ExactNefPoint p) {
        if (p.hw < CGAL::Gmpz(0)) {
            p.hx = -p.hx;
            p.hy = -p.hy;
            p.hw = -p.hw;
        }

        CGAL::Gmpz g = gcd3_gmpz(p.hx, p.hy, p.hw);

        if (g != CGAL::Gmpz(0) && g != CGAL::Gmpz(1)) {
            p.hx /= g;
            p.hy /= g;
            p.hw /= g;
        }

        return p;
    }

    template <typename NefPoint>
    static ExactNefPoint nef_point_to_exact(const NefPoint& pt) {
        ExactNefPoint p;
        p.hx = pt.hx();
        p.hy = pt.hy();
        p.hw = pt.hw();

        return normalize_exact_point(p);
    }

    static bool same_exact_point(
        const ExactNefPoint& a,
        const ExactNefPoint& b
    ) {
        return a.hx == b.hx &&
            a.hy == b.hy &&
            a.hw == b.hw;
    }

    static std::string exact_point_key(const ExactNefPoint& p) {
        std::ostringstream oss;
        oss << p.hx << "/" << p.hw
            << ","
            << p.hy << "/" << p.hw;

        return oss.str();
    }

    static std::string directed_edge_key_exact(
        const ExactNefPoint& a,
        const ExactNefPoint& b
    ) {
        return exact_point_key(a) + "->" + exact_point_key(b);
    }

    static std::string undirected_edge_key_exact(
        const ExactNefPoint& a,
        const ExactNefPoint& b
    ) {
        std::string ka = exact_point_key(a);
        std::string kb = exact_point_key(b);

        if (ka < kb) {
            return ka + "--" + kb;
        }

        return kb + "--" + ka;
    }

    static std::pair<double, double> exact_point_to_double(
        const ExactNefPoint& p
    ) {
        double hw = CGAL::to_double(p.hw);

        return {
            CGAL::to_double(p.hx) / hw,
            CGAL::to_double(p.hy) / hw
        };
    }

    template <typename NefPoint>
    static std::pair<double, double> nef_point_to_double_pair(
        const NefPoint& pt
    ) {
        double hw = CGAL::to_double(pt.hw());

        return {
            CGAL::to_double(pt.hx()) / hw,
            CGAL::to_double(pt.hy()) / hw
        };
    }

    static NFPFaceCycle extract_face_cycle_from_circulator(
        const Nef_polyhedron::Explorer& explorer,
        Nef_polyhedron::Explorer::Halfedge_around_face_const_circulator circ,
        bool face_marked,
        bool is_hole,
        bool unbounded
    ) {
        NFPFaceCycle cycle;
        cycle.face_marked = face_marked;
        cycle.is_hole = is_hole;
        cycle.unbounded = unbounded;

        if (unbounded) {
            return cycle;
        }

        Nef_polyhedron::Explorer::Halfedge_around_face_const_circulator start = circ;

        do {
            auto vh = explorer.source(circ);

            if (explorer.is_standard(vh)) {
                auto pt = explorer.point(vh);
                cycle.points.push_back(nef_point_to_exact(pt));
            }

            ++circ;
        } while (circ != start);

        return cycle;
    }

    static NFPTopology extract_nfp_topology(const Nef_polyhedron& nfp) {
        NFPTopology out;

        Nef_polyhedron::Explorer explorer = nfp.explorer();

        // ------------------------------------------------------------
        // 1. Extract all faces and all face-hole cycles.
        // ------------------------------------------------------------
        for (auto fit = explorer.faces_begin();
             fit != explorer.faces_end();
             ++fit) {

            bool face_marked = explorer.mark(fit);

            auto fc = explorer.face_cycle(fit);

            if (fc == Nef_polyhedron::Explorer::Halfedge_const_handle()) {
                NFPFaceCycle c;
                c.face_marked = face_marked;
                c.is_hole = false;
                c.unbounded = true;
                out.face_cycles.push_back(std::move(c));
            }
            else {
                Nef_polyhedron::Explorer::Halfedge_around_face_const_circulator circ(fc);
                NFPFaceCycle c = extract_face_cycle_from_circulator(
                    explorer,
                    circ,
                    face_marked,
                    false,
                    false
                );

                if (!c.points.empty()) {
                    out.face_cycles.push_back(std::move(c));
                }
            }

            // Extract holes of this face, regardless of whether the face itself
            // is marked or unmarked.
            for (auto hole_it = explorer.holes_begin(fit);
                 hole_it != explorer.holes_end(fit);
                 ++hole_it) {
                Nef_polyhedron::Explorer::Halfedge_around_face_const_circulator hcirc(hole_it);
                NFPFaceCycle h = extract_face_cycle_from_circulator(
                    explorer,
                    hcirc,
                    face_marked,
                    true,
                    false
                );

                if (!h.points.empty()) {
                    out.hole_cycles.push_back(std::move(h));
                }
            }
        }

        // ------------------------------------------------------------
        // 2. Extract all finite halfedges.
        //
        // For a halfedge e:
        //   explorer.face(e)       = face on the left side of e
        //   explorer.face(twin(e)) = face on the right side of e
        //
        // left_mark != right_mark means this edge separates forbidden/free
        // regions and is a strong NFP-boundary candidate.
        // ------------------------------------------------------------
        for (auto eit = explorer.halfedges_begin();
             eit != explorer.halfedges_end();
             ++eit) {

            auto s = explorer.source(eit);
            auto t = explorer.target(eit);

            // Skip rays/frame/non-finite geometry for now.
            if (!explorer.is_standard(s) || !explorer.is_standard(t)) {
                continue;
            }

            auto twin = explorer.twin(eit);

            auto left_face = explorer.face(eit);
            auto right_face = explorer.face(twin);

            auto left_fc = explorer.face_cycle(left_face);
            auto right_fc = explorer.face_cycle(right_face);

            NFPEdgeCell e;

            e.edge_marked = explorer.mark(eit);

            e.left_face_marked = explorer.mark(left_face);
            e.right_face_marked = explorer.mark(right_face);

            e.left_face_unbounded =
                (left_fc == Nef_polyhedron::Explorer::Halfedge_const_handle());

            e.right_face_unbounded =
                (right_fc == Nef_polyhedron::Explorer::Halfedge_const_handle());

            e.a = nef_point_to_exact(explorer.point(s));
            e.b = nef_point_to_exact(explorer.point(t));

            out.edges.push_back(std::move(e));
        }

        // ------------------------------------------------------------
        // 3. Extract all finite vertices.
        //
        // These are where point placements / lock-and-key placements may appear.
        // ------------------------------------------------------------
        for (auto vit = explorer.vertices_begin();
             vit != explorer.vertices_end();
             ++vit) {

            if (!explorer.is_standard(vit)) {
                continue;
            }

            NFPPointCell p;
            p.marked = explorer.mark(vit);
            p.p = nef_point_to_exact(explorer.point(vit));

            out.points.push_back(std::move(p));
        }

        return out;
    }

    static NFPPoint to_public_point(const ExactNefPoint& p) {
        auto d = exact_point_to_double(p);
        return { d.first, d.second };
    }

    static NFPSegment to_public_segment(const NFPSlitCandidate& s) {
        return {
            to_public_point(s.a),
            to_public_point(s.b)
        };
    }

    static double signed_area_of_exact_cycle(
        const std::vector<ExactNefPoint>& cycle
    ) {
        double area = 0.0;

        for (size_t i = 0; i < cycle.size(); ++i) {
            auto p = exact_point_to_double(cycle[i]);
            auto q = exact_point_to_double(cycle[(i + 1) % cycle.size()]);

            area += p.first * q.second - q.first * p.second;
        }

        return 0.5 * area;
    }

    static NFPCycleRole classify_boundary_cycle_role(
        const std::vector<ExactNefPoint>& cycle
    ) {
        return signed_area_of_exact_cycle(cycle) >= 0.0
            ? NFPCycleRole::OuterBoundary
            : NFPCycleRole::HoleBoundary;
    }

    static NFPResult make_public_result(const BuiltNFP& built) {
        NFPResult result;

        for (const auto& cycle : built.boundary_cycles) {
            NFPCycle out_cycle;
            out_cycle.role = classify_boundary_cycle_role(cycle);
            if (out_cycle.role == NFPCycleRole::OuterBoundary) {
                for (const auto& p : cycle) {
                    out_cycle.points.push_back(to_public_point(p));
                }
                result.cycles.push_back(std::move(out_cycle));
            }
        }

        for (const auto& pocket : built.feasible_pocket_cycles) {
            NFPCycle out_pocket;
            out_pocket.role = NFPCycleRole::FeasiblePocket;

            for (const auto& p : pocket) {
                out_pocket.points.push_back(to_public_point(p));
            }

            result.cycles.push_back(std::move(out_pocket));
        }

        for (const auto& slit : built.slit_candidates) {
            result.slits.push_back(to_public_segment(slit));
        }

        for (const auto& p : built.point_candidates) {
            result.isolated_points.push_back(to_public_point(p.p));
        }

        return result;
    }

    static void print_nfp_topology_summary(const NFPTopology& topo) {
        std::cout << "\nNFP topology summary\n";
        std::cout << "============================================================\n";

        std::cout << "Face cycles: " << topo.face_cycles.size() << "\n";
        for (size_t i = 0; i < topo.face_cycles.size(); ++i) {
            const auto& c = topo.face_cycles[i];

            std::cout << "  FaceCycle " << i
                      << " marked=" << c.face_marked
                      << " hole=" << c.is_hole
                      << " unbounded=" << c.unbounded
                      << " points=" << c.points.size()
                      << "\n";
        }

        std::cout << "Hole cycles: " << topo.hole_cycles.size() << "\n";
        for (size_t i = 0; i < topo.hole_cycles.size(); ++i) {
            const auto& h = topo.hole_cycles[i];

            std::cout << "  HoleCycle " << i
                      << " owner_face_marked=" << h.face_marked
                      << " points=" << h.points.size()
                      << "\n";
        }

        std::cout << "Halfedges: " << topo.edges.size() << "\n";
        for (size_t i = 0; i < topo.edges.size(); ++i) {
            const auto& e = topo.edges[i];

            bool mark_transition = (e.left_face_marked != e.right_face_marked);
            auto a = exact_point_to_double(e.a);
            auto b = exact_point_to_double(e.b);
            std::cout << "  Edge " << i
                      << " edge_marked=" << e.edge_marked
                      << " left_face_marked=" << e.left_face_marked
                      << " right_face_marked=" << e.right_face_marked
                      << " transition=" << mark_transition
                      << " a=(" << a.first << ", " << a.second << ")"
                      << " b=(" << b.first << ", " << b.second << ")"
                      << "\n";
        }

        std::cout << "Vertices: " << topo.points.size() << "\n";
        for (size_t i = 0; i < topo.points.size(); ++i) {
            const auto& p = topo.points[i];
            const auto& p1 = exact_point_to_double(p.p);
            std::cout << "  Vertex " << i
                      << " marked=" << p.marked
                      << " p=(" << p1.first << ", " << p1.second << ")"
                      << "\n";
        }

        std::cout << "============================================================\n\n";
    }

 

    static std::string point_key(const std::pair<double, double>& p) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(12)
            << p.first << "," << p.second;
        return oss.str();
    }

    static double dist2(
        const std::pair<double, double>& a,
        const std::pair<double, double>& b
    ) {
        double dx = a.first - b.first;
        double dy = a.second - b.second;
        return dx * dx + dy * dy;
    }

    static std::string point_key_double(const std::pair<double, double>& p) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(12)
            << p.first << "," << p.second;
        return oss.str();
    }

    static std::string directed_edge_key_double(
        const std::pair<double, double>& a,
        const std::pair<double, double>& b
    ) {
        return point_key_double(a) + "->" + point_key_double(b);
    }

    static bool same_point_double(
        const std::pair<double, double>& a,
        const std::pair<double, double>& b
    ) {
        return point_key_double(a) == point_key_double(b);
    }

    static std::string undirected_edge_key_double(
        const std::pair<double, double>& a,
        const std::pair<double, double>& b
    ) {
        std::string ka = point_key_double(a);
        std::string kb = point_key_double(b);

        if (ka < kb) {
            return ka + "--" + kb;
        }

        return kb + "--" + ka;
    }

    static std::vector<std::vector<ExactNefPoint>>
        stitch_boundary_cycles_from_edges(const NFPTopology& topo) {
        std::vector<DirectedBoundaryEdge> edges;
        std::set<std::string> seen_edges;

        edges.reserve(topo.edges.size());

        for (const auto& e : topo.edges) {
            if (e.left_face_unbounded || e.right_face_unbounded) {
                continue;
            }

            if (e.left_face_marked == e.right_face_marked) {
                continue;
            }

            DirectedBoundaryEdge de;

            if (e.left_face_marked && !e.right_face_marked) {
                de.a = e.a;
                de.b = e.b;
            }
            else {
                de.a = e.b;
                de.b = e.a;
            }

            if (same_exact_point(de.a, de.b)) {
                continue;
            }

            de.a_key = exact_point_key(de.a);
            de.b_key = exact_point_key(de.b);

            std::string ekey = directed_edge_key_exact(de.a, de.b);

            if (seen_edges.count(ekey)) {
                continue;
            }

            seen_edges.insert(ekey);

            edges.push_back(std::move(de));
        }

        std::map<std::string, std::vector<size_t>> outgoing;

        for (size_t i = 0; i < edges.size(); ++i) {
            outgoing[edges[i].a_key].push_back(i);
        }

        std::vector<std::vector<ExactNefPoint>> cycles;

        for (size_t i = 0; i < edges.size(); ++i) {
            if (edges[i].used) {
                continue;
            }

            std::vector<ExactNefPoint> cycle;

            size_t current = i;
            std::string start_key = edges[current].a_key;

            while (true) {
                if (edges[current].used) {
                    break;
                }

                edges[current].used = true;

                cycle.push_back(edges[current].a);

                std::string next_key = edges[current].b_key;

                if (next_key == start_key) {
                    break;
                }

                auto found = outgoing.find(next_key);

                if (found == outgoing.end()) {
                    cycle.push_back(edges[current].b);
                    break;
                }

                bool advanced = false;

                for (size_t candidate_idx : found->second) {
                    if (!edges[candidate_idx].used) {
                        current = candidate_idx;
                        advanced = true;
                        break;
                    }
                }

                if (!advanced) {
                    cycle.push_back(edges[current].b);
                    break;
                }
            }

            if (cycle.size() >= 2) {
                cycles.push_back(std::move(cycle));
            }
        }

        return cycles;
    }


    static std::vector<std::vector<ExactNefPoint>>
    extract_feasible_pocket_cycles(const NFPTopology& topo) {
        std::vector<std::vector<ExactNefPoint>> pockets;

        for (const auto& c : topo.face_cycles) {
            if (c.unbounded) {
                continue;
            }

            if (c.face_marked) {
                continue;
            }

            if (c.points.size() < 3) {
                continue;
            }

            pockets.push_back(c.points);
        }

        return pockets;
    }

    static std::vector<NFPSlitCandidate>
        extract_slit_candidates(const NFPTopology& topo) {
        std::vector<NFPSlitCandidate> slits;
        std::set<std::string> seen_slits;

        for (const auto& e : topo.edges) {
            if (e.left_face_unbounded || e.right_face_unbounded) {
                continue;
            }

            if (same_exact_point(e.a, e.b)) {
                continue;
            }

            bool embedded_in_forbidden_area =
                e.left_face_marked && e.right_face_marked;

            bool edge_itself_available =
                !e.edge_marked;

            if (embedded_in_forbidden_area && edge_itself_available) {
                std::string key = undirected_edge_key_exact(e.a, e.b);

                if (seen_slits.count(key)) {
                    continue;
                }

                seen_slits.insert(key);

                NFPSlitCandidate s;
                s.a = e.a;
                s.b = e.b;

                slits.push_back(std::move(s));
            }
        }

        return slits;
    }

    static std::vector<NFPPointCandidate>
        extract_point_candidates(
            const NFPTopology& topo,
            const std::vector<std::vector<ExactNefPoint>>& boundary_cycles,
            const std::vector<NFPSlitCandidate>& slit_candidates
        ) {
        std::set<std::string> boundary_vertices;

        for (const auto& cycle : boundary_cycles) {
            for (const auto& p : cycle) {
                boundary_vertices.insert(exact_point_key(p));
            }
        }

        std::set<std::string> slit_vertices;

        for (const auto& s : slit_candidates) {
            slit_vertices.insert(exact_point_key(s.a));
            slit_vertices.insert(exact_point_key(s.b));
        }

        std::vector<NFPPointCandidate> points;

        for (const auto& p : topo.points) {
            if (p.marked) {
                continue;
            }

            std::string key = exact_point_key(p.p);

            if (boundary_vertices.count(key)) {
                continue;
            }

            if (slit_vertices.count(key)) {
                continue;
            }

            NFPPointCandidate q;
            q.p = p.p;

            points.push_back(std::move(q));
        }

        return points;
    }

    static BuiltNFP build_nfp_from_topology(const NFPTopology& topo) {
        BuiltNFP result;

        result.boundary_cycles = stitch_boundary_cycles_from_edges(topo);

        result.feasible_pocket_cycles = extract_feasible_pocket_cycles(topo);

        result.slit_candidates = extract_slit_candidates(topo);

        result.point_candidates = extract_point_candidates(
            topo,
            result.boundary_cycles,
            result.slit_candidates
        );

        return result;
    }

    static void print_built_nfp_summary(const BuiltNFP& built) {
        std::cout << "\nBuilt NFP summary\n";
        std::cout << "============================================================\n";

        std::cout << "Boundary cycles: " << built.boundary_cycles.size() << "\n";
        for (size_t i = 0; i < built.boundary_cycles.size(); ++i) {
            std::cout << "  Boundary cycle " << i
                      << " points=" << built.boundary_cycles[i].size()
                      << "\n";
        }

        std::cout << "Feasible pocket cycles: "
                  << built.feasible_pocket_cycles.size()
                  << "\n";
        for (size_t i = 0; i < built.feasible_pocket_cycles.size(); ++i) {
            std::cout << "  Pocket cycle " << i
                      << " points=" << built.feasible_pocket_cycles[i].size()
                      << "\n";
        }

        std::cout << "Slit candidates: " << built.slit_candidates.size() << "\n";
        for (size_t i = 0; i < built.slit_candidates.size(); ++i) {
            const auto& s = built.slit_candidates[i];
            const auto sa = exact_point_to_double(s.a);
            const auto sb = exact_point_to_double(s.b);
            std::cout << "  Slit " << i
                      << " a=(" << sa.first << ", " << sa.second << ")"
                      << " b=(" << sb.first << ", " << sb.second << ")"
                      << "\n";
        }

        std::cout << "Point candidates: " << built.point_candidates.size() << "\n";
        for (size_t i = 0; i < built.point_candidates.size(); ++i) {
            const auto& p = built.point_candidates[i];
            const auto pp = exact_point_to_double(p.p);
            std::cout << "  Point " << i
                      << " p=(" << pp.first << ", " << pp.second << ")"
                      << "\n";
        }

        std::cout << "============================================================\n\n";
    }

    NFPResult processNFP(
            const std::vector<std::pair<double,double>>& polyA, 
            const std::vector<std::pair<double,double>>& polyB) {
        
        // Convert input to CGAL Polygon_2
        Polygon_2 cgalPolyA, cgalPolyB;
        for (const auto& p : polyA) {
            cgalPolyA.push_back(Point_2(p.first, p.second));
        }
        for (const auto& p : polyB) {
            cgalPolyB.push_back(Point_2(p.first, p.second));
        }

        //trying to sanitize polygon

        sanitize_polygon(cgalPolyA);
        sanitize_polygon(cgalPolyB);


        if (!cgalPolyA.is_simple())
        {
            std::cerr << "ERROR: Polygon A (fixed Polygon) is not simple!" << std::endl;
            for (auto &p:cgalPolyA){
                std::cout << CGAL::to_double(p.x()) << " " << CGAL::to_double(p.y()) << "\n";
            }
            exit(-1);
        }


        if (!cgalPolyB.is_simple())
        {
            std::cerr << "ERROR: Polygon B (rotating Polygon) is not simple!" << std::endl;
            for (auto &p:cgalPolyB){
                std::cout << CGAL::to_double(p.x()) << " " << CGAL::to_double(p.y()) << "\n";
            }
            exit(-1);
        }
        

        //ensure orientation is counter clockwise
        if (cgalPolyA.is_clockwise_oriented()) {
            cgalPolyA.reverse_orientation();
        }

        if (cgalPolyB.is_clockwise_oriented()) {
            cgalPolyB.reverse_orientation();
        }

        // Get NFP as Nef Polyhedron
        Nef_polyhedron nfp = get_nfp(0, cgalPolyA, 1, cgalPolyB);

        NFPTopology topo = extract_nfp_topology(nfp);

        //print_nfp_topology_summary(topo);

        BuiltNFP built = build_nfp_from_topology(topo);

        //print_built_nfp_summary(built);

        // Current compatibility return:
        // Return boundary cycles first, then positive-area feasible pockets.
        // Slits and points are printed but not returned because your current return
        // type cannot represent them.
        /*
        std::vector<std::vector<std::pair<double,double>>> nfp_polygons;

        for (const auto& c : built.boundary_cycles) {
            if (!c.empty()) {
                nfp_polygons.push_back(c);
            }
        }
        */

        /*
        for (const auto& p : built.feasible_pocket_cycles) {
            if (!p.empty()) {
                nfp_polygons.push_back(p);
            }
        }
        */

        return make_public_result(built);


        /*
        // Extract polygons by converting Nef back to polygon set
        std::vector<std::vector<std::pair<double,double>>> nfp_polygons;
        
        Nef_polyhedron::Explorer explorer = nfp.explorer();
        
        return extract_polygon(explorer);
        */
        
    }

    std::vector<std::vector<std::pair<double,double>>> extract_polygon( Nef_polyhedron::Explorer &explorer){

        std::vector<std::vector<std::pair<double,double>>> nfp_polygons;
        for (auto fit = explorer.faces_begin(); fit != explorer.faces_end(); ++fit) {
            // We want faces that are MARKED (part of the NFP interior)
            if (!explorer.mark(fit)) {
                continue;
            }
            
            // Get the face cycle (outer boundary of this face)
            auto fc = explorer.face_cycle(fit);
            
            // Check if this face has a valid cycle (not the unbounded face)
            if (fc == Nef_polyhedron::Explorer::Halfedge_const_handle()) {
                continue;
            }
            
            // Create circulator from the face cycle halfedge
            Nef_polyhedron::Explorer::Halfedge_around_face_const_circulator 
                circ(fc), start(fc);
            
            std::vector<std::pair<double,double>> polygon;
            
            // Walk the cycle
            do {
                auto vh = explorer.source(circ);
                
                // Only process standard (finite) vertices
                if (explorer.is_standard(vh)) {
                    auto pt = explorer.point(vh);
                    double hw = CGAL::to_double(pt.hw());
                    double x = CGAL::to_double(pt.hx()) / hw;
                    double y = CGAL::to_double(pt.hy()) / hw;
                    polygon.emplace_back(x, y);
                }
                
                ++circ;
            } while (circ != start);
            
            if (!polygon.empty()) {
                nfp_polygons.push_back(polygon);
            }
            
            // Also process holes of this face just in case if exists (inner boundaries)
            for (auto hole_it = explorer.holes_begin(fit);
                 hole_it != explorer.holes_end(fit); ++hole_it) {
                
                Nef_polyhedron::Explorer::Halfedge_around_face_const_circulator
                    hcirc(hole_it), hstart(hole_it);
                
                std::vector<std::pair<double,double>> hole_polygon;
                
                do {
                    auto vh = explorer.source(hcirc);
                    if (explorer.is_standard(vh)) {
                        auto pt = explorer.point(vh);
                        double hw = CGAL::to_double(pt.hw());
                        double x = CGAL::to_double(pt.hx()) / hw;
                        double y = CGAL::to_double(pt.hy()) / hw;
                        hole_polygon.emplace_back(x, y);
                    }
                    ++hcirc;
                } while (hcirc != hstart);
                
                if (!hole_polygon.empty()) {
                    nfp_polygons.push_back(hole_polygon);
                }
            }
        }
        return nfp_polygons;
    }
}
#endif // NFP_GENERATOR_H

