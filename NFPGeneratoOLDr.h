#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Nef_polyhedron_2.h>
#include <CGAL/minkowski_sum_2.h>
#include <CGAL/partition_2.h>
#include <CGAL/Aff_transformation_2.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Extended_homogeneous.h>
#include <CGAL/Homogeneous.h>  // ADD THIS

#include <vector>
#include <map>
#include <iostream>

// Kernel typedefs
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2                                   Point_2;
typedef CGAL::Polygon_2<Kernel>                           Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel>                Polygon_with_holes_2;

// Partition traits - this ensures output polygons match our Polygon_2 type
typedef CGAL::Partition_traits_2<Kernel>                  Partition_traits;
typedef Partition_traits::Polygon_2                       Partition_polygon;  // Uses std::list internally

// Nef kernel typedefs
typedef CGAL::Extended_homogeneous<CGAL::Exact_rational>  Nef_kernel;
typedef CGAL::Nef_polyhedron_2<Nef_kernel>                Nef_polyhedron;
typedef CGAL::Homogeneous<CGAL::Exact_rational>           Standard_kernel;
typedef Standard_kernel::Point_2                          Standard_point;

typedef int ShapeID; // Simple alias for shape identifiers
class NFPCache {
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
        CGAL::Aff_transformation_2<Kernel> t(CGAL::SCALING, -1.0);
        Polygon_2 p_reflected;
        for(auto v : p) p_reflected.push_back(t.transform(v));
        
        // Reflection reverses orientation, fix it to CCW
        if (p_reflected.is_clockwise_oriented()) p_reflected.reverse_orientation();
        return p_reflected;
    }

public:
    // THE CORE PIPELINE
    Nef_polyhedron get_nfp(ShapeID idA, Polygon_2& polyA, 
                           ShapeID idB, Polygon_2& polyB) {
        
        // 1. Check Cache
        std::pair<ShapeID, ShapeID> key = {idA, idB};
        if (cache.find(key) != cache.end()) {
            // Cache Hit! Return instantly.
            return cache[key];
        }

        std::cout << "[NFP] Computing new NFP for " << idA << " vs " << idB << "..." << std::endl;
        if (polyA.is_clockwise_oriented()) {
            polyA.reverse_orientation();
        }

    
        // 2. Prepare Geometry
        // Decompose A (Stationary)
        std::vector<Polygon_2> partsA = decompose(polyA);
        
        // Reflect and Decompose B (Orbiting) -> (-B)
        Polygon_2 polyB_reflected = reflect_polygon(polyB);

        if (polyB_reflected.is_clockwise_oriented()) {
            polyB_reflected.reverse_orientation();
        }
        
        std::vector<Polygon_2> partsB = decompose(polyB_reflected);

        // 3. Initialize Accumulator (Empty Nef Polyhedron)
        Nef_polyhedron nfp_accumulator; // Starts empty (false everywhere)
        
        // 4. Pairwise Sum & Union Loop
        // Complexity: N_parts_A * N_parts_B
        for (const auto& subA : partsA) {
            for (const auto& subB : partsB) {
                
                // A. Compute Convex Minkowski Sum
                Polygon_with_holes_2 convex_sum = CGAL::minkowski_sum_2(subA, subB);
                            
                // B. Convert to Nef Polyhedron (with Kernel Conversion)
                const auto& outer = convex_sum.outer_boundary();

                // Use Standard_kernel points, NOT Nef_kernel points
                std::vector<Standard_point> nef_points;
                nef_points.reserve(outer.size());

                for (auto vit = outer.vertices_begin(); vit != outer.vertices_end(); ++vit) {
                    // Convert EPECK coordinates to Exact_rational for Homogeneous kernel
                    double x_d = CGAL::to_double(vit->x());
                    double y_d = CGAL::to_double(vit->y());
                    nef_points.emplace_back(Standard_kernel::RT(x_d), Standard_kernel::RT(y_d));
                }

                // Now construction works
                Nef_polyhedron nef_sub(nef_points.begin(), nef_points.end());
                
                // C. Union with total result
                if (nfp_accumulator.is_empty()) {
                    nfp_accumulator = nef_sub;
                } else {
                    nfp_accumulator += nef_sub; // += is usually more efficient than +
                }
            }
        }

        // 5. Store in Cache
        cache[key] = nfp_accumulator;
        return nfp_accumulator;
    }

    // Helper to clear cache if needed (e.g. low memory)
    void clear() {
        cache.clear();
    }

    std::vector<std::vector<std::pair<double,double>>> processNFP(std::vector<std::pair<double,double>> &polyA, std::vector<std::pair<double,double>> &polyB) {
        // Convert input to CGAL Polygon_2
        Polygon_2 cgalPolyA, cgalPolyB;
        for (const auto& p : polyA) {
            cgalPolyA.push_back(Point_2(p.first, p.second));
        }
        for (const auto& p : polyB) {
            cgalPolyB.push_back(Point_2(p.first, p.second));
        }

        // Get NFP as Nef Polyhedron
        Nef_polyhedron::Explorer explorer = get_nfp(0, cgalPolyA, 1, cgalPolyB).explorer();

        // Extract outer boundary of NFP
        std::vector<std::vector<std::pair<double,double>>> nfp_polygons;

        // 2. Iterate Faces using the Explorer
        for (auto fit = explorer.faces_begin(); fit != explorer.faces_end(); ++fit) {
        
            // 3. Check if the face is marked (Part of the NFP)
            if (explorer.mark(fit)) { 
                
                // 4. Iterate over the "Holes" (Boundary Cycles) of this face
                auto hole_it = explorer.holes_begin(fit);
                auto hole_end = explorer.holes_end(fit);
            
                for (; hole_it != hole_end; ++hole_it) {
                    
                    // The iterator points to a Halfedge that starts the cycle
                    // We create a circulator to walk around the polygon loop
                    Nef_polyhedron::Explorer::Halfedge_around_face_const_circulator circ(hole_it);
                    Nef_polyhedron::Explorer::Halfedge_around_face_const_circulator start = circ;
                
                    std::vector<std::pair<double,double>> polygon;
                
                    // 5. Walk the cycle
                    do {
                        // ACCESSOR SYNTAX: explorer.point(explorer.source(edge))
                        auto p = explorer.point(explorer.source(circ));
                        
                        polygon.emplace_back(CGAL::to_double(p.x()), CGAL::to_double(p.y()));
                        
                        ++circ; // Move to next edge in the loop
                    } while (circ != start);
                
                    nfp_polygons.push_back(polygon);
                }
            }
        }

        // For simplicity, return only the first polygon found
        if (!nfp_polygons.empty()) {
            return nfp_polygons;
        } else {
            return {};
        }
    }
};