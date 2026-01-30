#include <iostream>
#include <iomanip>
#include "common/geometry_util.hpp"
#include "common/NFPGenerator.h"



using Point = std::pair<double, double>;
using Polygon = std::vector<Point>;

//using namespace GeometryUtil;
//
//
//void printPolygon(const std::string& name, const Polygon& poly) {
//    std::cout << name << " (" << poly.vertices.size() << " vertices):\n";
//    for (size_t i = 0; i < poly.vertices.size(); i++) {
//        std::cout << "  v" << i << ": (" << std::fixed << std::setprecision(2)
//                  << poly.vertices[i].x << ", " << poly.vertices[i].y << ")\n";
//    }
//    std::cout << std::endl;
//}
//
//void printNFP(const std::vector<std::vector<Point>>& nfp) {
//    std::cout << "NFP: " << nfp.size() << " cycle(s)\n";
//    std::cout << std::string(60, '=') << "\n\n";
//
//    for (size_t i = 0; i < nfp.size(); i++) {
//        std::cout << "Cycle " << (i + 1) << " (" << nfp[i].size() << " points):\n";
//        for (size_t j = 0; j < nfp[i].size(); j++) {
//            std::cout << "  p" << j << ": (" << std::fixed << std::setprecision(4)
//                      << nfp[i][j].x << ", " << nfp[i][j].y << ")\n";
//        }
//        std::cout << std::endl;
//    }
//}
//
//void test1_convex_convex() {
//    std::cout << "\n" << std::string(70, '=') << "\n";
//    std::cout << "TEST 1: Convex-Convex (Rectangle + Triangle)\n";
//    std::cout << std::string(70, '=') << "\n\n";
//
//    std::vector<Point> verts_a = {
//        Point(0, 0), Point(4, 0), Point(4, 3), Point(0, 3)
//    };
//
//    std::vector<Point> verts_b = {
//        Point(0, 0), Point(2, 0), Point(1, 2)
//    };
//
//    Polygon poly_a(verts_a);
//    Polygon poly_b(verts_b);
//
//    printPolygon("Polygon A (rectangle)", poly_a);
//    printPolygon("Polygon B (triangle)", poly_b);
//
//    NFPGenerator gen;
//    auto nfp = gen.generate(poly_a, poly_b);
//
//    printNFP(nfp);
//}
//

//

void printPolygon(const std::string& name, const Polygon& poly) {
    std::cout << name << " (" << poly.size() << " vertices):\n";
    for (size_t i = 0; i < poly.size(); i++) {
        std::cout << "  v" << i << ": (" << std::fixed << std::setprecision(2)
            << poly[i].first << ", " << poly[i].second << ")\n";
    }
    std::cout << std::endl;
}

void printNFP(const std::vector<std::vector<Point>>& nfp) {
    std::cout << "NFP: " << nfp.size() << " cycle(s)\n";
    std::cout << std::string(60, '=') << "\n\n";

    for (size_t i = 0; i < nfp.size(); i++) {
        std::cout << "Cycle " << (i + 1) << " (" << nfp[i].size() << " points):\n";
        for (size_t j = 0; j < nfp[i].size(); j++) {
            std::cout << "  p" << j << ": (" << std::fixed << std::setprecision(4)
                << nfp[i][j].first << ", " << nfp[i][j].second << ")\n";
        }
        std::cout << std::endl;
    }
}


void test2_L_shape_square() {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "TEST 2: L-Shape + Square (CRITICAL TEST)\n";
    std::cout << std::string(70, '=') << "\n\n";

    // L-shape
    std::vector<Point> verts_a = {
        Point(0, 0),
        Point(3, 0),
        Point(3, 2),
        Point(2, 2),
        Point(2, 3),
        Point(0, 3)
    };

    // Square
    std::vector<Point> verts_b = {
        Point(0, 0),
        Point(1, 0),
        Point(1, 1),
        Point(0, 1)
    };

    Polygon poly_a(verts_a);
    Polygon poly_b(verts_b);

    printPolygon("Polygon A (L-shape)", poly_a);
    printPolygon("Polygon B (square)", poly_b);

    NFPGeneratorCGAL gen;
    auto nfp = gen.processNFP(poly_a, poly_b);

    printNFP(nfp);

    // Verify no self-intersection
    
    if (!nfp.empty()) {
        std::cout << "NFP generated successfully with " << nfp[0].size() << " vertices\n";
        std::cout << "Visual inspection needed to verify correctness.\n";
    }
}


void test3_L_shape_L() {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "TEST 3: L-Shape + L-Shape (CRITICAL TEST)\n";
    std::cout << std::string(70, '=') << "\n\n";

    // L-1
    std::vector<Point> verts_a = {
        Point(0, 0),
        Point(3, 0),
        Point(3, 2),
        Point(2, 2),
        Point(2, 3),
        Point(0, 3)
    };

    // L-2
    std::vector<Point> verts_b = {
        Point(0, 0),
        Point(3, 0),
        Point(3, 2),
        Point(2, 2),
        Point(2, 3),
        Point(0, 3)
    };

    Polygon poly_a(verts_a);
    Polygon poly_b(verts_b);

    printPolygon("Polygon A (L-shape)", poly_a);
    printPolygon("Polygon B (L-shape)", poly_b);

    NFPGeneratorCGAL gen;
    auto nfp = gen.processNFP(poly_a, poly_b);

    printNFP(nfp);

    // Verify no self-intersection
    
    if (!nfp.empty()) {
        std::cout << "NFP generated successfully with " << nfp[0].size() << " vertices\n";
        std::cout << "Visual inspection needed to verify correctness.\n";
    }
}
//
//
//
//void test4() {
//    std::cout << "\n" << std::string(70, '=') << "\n";
//    std::cout << "TEST 4: Complex Non-Convex polygon\n";
//    std::cout << std::string(70, '=') << "\n\n";
//
//    //PolyA
//    //std::vector<Point> verts_a = {
//    //    Point(0, 0),
//    //    Point(3, 0),
//    //    Point(2, -2),
//    //    Point(3, -4),
//    //    Point(3, -5),
//    //    Point(1, -5),
//    //    Point(-1, -3),
//    //    Point(-1, -1)
//    //};
//
//    //PolyA
//    int square_size = 12;
//    int concavLength = 4;
//    int concavheight = 7;
//    int concavStart = 3;
//    std::vector<Point> verts_a = {
//        Point(0, 0),
//        Point(square_size, 0),
//        Point(square_size, concavStart),
//        Point(square_size-concavLength, concavStart),
//        Point(square_size-concavLength, concavStart+concavheight),
//        Point(square_size, concavStart+concavheight),
//        Point(square_size, square_size),
//        Point(0, square_size),
//    };
//
//    //polyB
//    std::vector<Point> verts_b = {
//        Point(0, 0),
//        Point(2, 1),
//        Point(4, 0),
//        Point(3, 2),
//        Point(4, 5),
//        Point(2, 4),
//        Point(0, 5),
//        Point(1, 3)
//    };
//
//    Polygon poly_a(verts_a);
//    Polygon poly_b(verts_b);
//
//    printPolygon("Polygon A ", poly_a);
//    printPolygon("Polygon B ", poly_b);
//
//    NFPGenerator gen;
//    auto nfp = gen.generate(poly_a, poly_b);
//
//    printNFP(nfp);
//
//    // Verify no self-intersection
//    
//    if (!nfp.empty()) {
//        std::cout << "NFP generated successfully with " << nfp[0].size() << " vertices\n";
//        std::cout << "Visual inspection needed to verify correctness.\n";
//    }
//}
//
//int main() {
//    std::cout << "\n";
//    std::cout << "╔══════════════════════════════════════════════════════════════════╗\n";
//    std::cout << "║        NFP Algorithm - Faithful Implementation Test             ║\n";
//    std::cout << "║           Bennell & Song (2008) - Complete Version              ║\n";
//    std::cout << "╚══════════════════════════════════════════════════════════════════╝\n";
//
//    try {
//        //test1_convex_convex();
//        test2_L_shape_square();
//        test3_L_shape_L();
//        //test4();
//
//        std::cout << "\n" << std::string(70, '=') << "\n";
//        std::cout << "Tests completed.\n";
//        std::cout << std::string(70, '=') << "\n\n";
//
//    } catch (const std::exception& e) {
//        std::cerr << "\nError: " << e.what() << std::endl;
//        return 1;
//    }
//
//    return 0;
//}
//

void test5() {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "TEST 4: Complex Non-Convex polygon with degeneracy\n";
    std::cout << std::string(70, '=') << "\n\n";
//

    //PolyA
    //std::vector<Point> verts_a = {
    //    Point(0, 0),
    //    Point(3, 0),
    //    Point(2, -2),
    //    Point(3, -4),
    //    Point(3, -5),
    //    Point(1, -5),
    //    Point(-1, -3),
    //    Point(-1, -1)
    //};
//

    //PolyA
    int square_size = 12;
    int concavLength = 4;
    int concavheight = 5;
    int concavStart = 3;
    std::vector<Point> verts_a = {
        Point(0, 0),
        Point(square_size, 0),
        Point(square_size, concavStart),
        Point(square_size-concavLength, concavStart),
        Point(square_size-concavLength, concavStart+concavheight),
        Point(square_size, concavStart+concavheight),
        Point(square_size, square_size),
        Point(0, square_size),
    };
//

    //polyB
    std::vector<Point> verts_b = {
        Point(0, 0),
        Point(2, 1),
        Point(4, 0),
        Point(3, 2),
        Point(4, 5),
        Point(2, 4),
        Point(0, 5),
        Point(1, 3)
    };
//

    Polygon poly_a(verts_a);
    Polygon poly_b(verts_b);
//

    printPolygon("Polygon A ", poly_a);
    printPolygon("Polygon B ", poly_b);
//

    NFPGeneratorCGAL gen;
    auto nfp = gen.processNFP(poly_a, poly_b);
//

    printNFP(nfp);
//

    // Verify no self-intersection
    
    if (!nfp.empty()) {
        std::cout << "NFP generated successfully with " << nfp[0].size() << " vertices\n";
        std::cout << "Visual inspection needed to verify correctness.\n";
    }
}


void test4() {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "TEST 4: Complex Non-Convex polygon\n";
    std::cout << std::string(70, '=') << "\n\n";
    

    //PolyA
    std::vector<Point> verts_a = {
        Point(0, 0),
        Point(3, 0),
        Point(2, -2),
        Point(3, -4),
        Point(3, -5),
        Point(1, -5),
        Point(-1, -3),
        Point(-1, -1)
    };
    

    //polyB
    std::vector<Point> verts_b = {
        Point(0, 0),
        Point(2, 1),
        Point(4, 0),
        Point(3, 2),
        Point(4, 5),
        Point(2, 4),
        Point(0, 5),
        Point(1, 3)
    };
    //

    Polygon poly_a(verts_a);
    Polygon poly_b(verts_b);
    //

    printPolygon("Polygon A ", poly_a);
    printPolygon("Polygon B ", poly_b);
    //

    NFPGeneratorCGAL gen;
    auto nfp = gen.processNFP(poly_a, poly_b);
    //

    printNFP(nfp);
    //

    if (!nfp.empty()) {
        std::cout << "NFP generated successfully with " << nfp[0].size() << " vertices\n";
        std::cout << "Visual inspection needed to verify correctness.\n";
    }
}


void test6() {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "TEST 6: Complex Non-Convex polygon with degeneracy 2\n";
    std::cout << std::string(70, '=') << "\n\n";

    //PolyA
    int square_size = 15;
    int concavLength = 8;
    int concavheight = 5;
    int concavStart = 3;
    std::vector<Point> verts_a = {
        Point(0, 0),
        Point(square_size, 0),
        Point(square_size, concavStart),
        Point(square_size-concavLength, concavStart),
        Point(square_size-concavLength, concavStart+2*concavheight),
        Point(square_size-concavLength/2, concavStart+2*concavheight),
        Point(square_size-concavLength/2, concavStart+concavheight),
        Point(square_size, concavStart+concavheight),
        Point(square_size, square_size),
        Point(0, square_size),
    };
    

    //polyB
    std::vector<Point> verts_b = {
        Point(0, 0),
        Point(2, 1),
        Point(4, 0),
        Point(3, 2),
        Point(4, 5),
        Point(2, 4),
        Point(0, 5),
        Point(1, 3)
    };
//

    Polygon poly_a(verts_a);
    Polygon poly_b(verts_b);
//

    printPolygon("Polygon A ", poly_a);
    printPolygon("Polygon B ", poly_b);
//

    NFPGeneratorCGAL gen;
    auto nfp = gen.processNFP(poly_a, poly_b);
//

    printNFP(nfp);
//

    // Verify no self-intersection
    
    if (!nfp.empty()) {
        std::cout << "NFP generated successfully with " << nfp[0].size() << " vertices\n";
        std::cout << "Visual inspection needed to verify correctness.\n";
    }
}


int main() {
    std::cout << "NestingGraph NFP Test Program\n";
    test2_L_shape_square();
    test3_L_shape_L();
    test4();
    test5();
    test6();
    return 0;
}