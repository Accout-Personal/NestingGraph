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

#include <vector>
#include <map>
#include <iostream>
#include <sstream>

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


class NFPGeneratorCGAL {
private:
    // Map key: Pair of (Stationary_ID, Orbiting_ID)
    // Value: The robust Nef Polyhedron representing the NFP
    std::map<std::pair<ShapeID, ShapeID>, Nef_polyhedron> cache;

    // Helper: Decompose a concave polygon into convex pieces
    std::vector<Polygon_2> decompose(const Polygon_2& p) {
        // Use partition traits polygons for output (they use std::list internally)
        std::vector<Partition_polygon> partition_parts;
        
        CGAL::optimal_convex_partition_2(
            p.vertices_begin(), 
            p.vertices_end(), 
            std::back_inserter(partition_parts),
            Partition_traits());  // Pass the traits!
        
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

public:
    // THE CORE PIPELINE
    Nef_polyhedron get_nfp(ShapeID idA, Polygon_2& polyA,
        ShapeID idB, Polygon_2& polyB) {

        std::pair<ShapeID, ShapeID> key = { idA, idB };

        std::vector<Polygon_2> partsA = decompose(polyA);
        Polygon_2 polyB_reflected = reflect_polygon(polyB);

        if (polyB_reflected.is_clockwise_oriented()) {
            polyB_reflected.reverse_orientation();
        }

        std::vector<Polygon_2> partsB = decompose(polyB_reflected);

        // Store individual Minkowski sums
        Nef_polyhedron U_Open(Nef_polyhedron::EMPTY);

        for (const auto& subA : partsA) {
            for (const auto& subB : partsB) {
                Polygon_2 sum = CGAL::minkowski_sum_2(subA, subB).outer_boundary();

                std::vector<Standard_point> pts;
                for (auto v = sum.vertices_begin(); v != sum.vertices_end(); ++v)
                    pts.push_back(convert_to_standard_point(*v));

                Nef_polyhedron nef_sub(pts.begin(), pts.end(), Nef_polyhedron::EXCLUDED);
                U_Open = U_Open.join(nef_sub);
            }
        }

        cache[key] = U_Open;
        return U_Open;
    }

    // Helper to clear cache if needed (e.g. low memory)
    void clear() {
        cache.clear();
    }


    std::vector<std::vector<std::pair<double,double>>> processNFP(
            std::vector<std::pair<double,double>>& polyA, 
            std::vector<std::pair<double,double>>& polyB) {
        
        // Convert input to CGAL Polygon_2
        Polygon_2 cgalPolyA, cgalPolyB;
        for (const auto& p : polyA) {
            cgalPolyA.push_back(Point_2(p.first, p.second));
        }
        for (const auto& p : polyB) {
            cgalPolyB.push_back(Point_2(p.first, p.second));
        }

        // Get NFP as Nef Polyhedron
        Nef_polyhedron nfp = get_nfp(0, cgalPolyA, 1, cgalPolyB);

        
        // Extract polygons by converting Nef back to polygon set
        std::vector<std::vector<std::pair<double,double>>> nfp_polygons;
        
        Nef_polyhedron::Explorer explorer = nfp.explorer();
        
        
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
        cache.clear();
        return nfp_polygons;
    }
};

#endif // NFP_GENERATOR_H
