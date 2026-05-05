#include <iostream>
#include <iomanip>
#include "common/geometry_util.hpp"
#include "nfp/NFPGenerator.hpp"


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
//    auto nfp = nfp::generate(poly_a, poly_b);
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

static const char* cycleRoleToString(nfp::NFPCycleRole role) {
    switch (role) {
    case nfp::NFPCycleRole::OuterBoundary:
        return "OuterBoundary";
    case nfp::NFPCycleRole::HoleBoundary:
        return "HoleBoundary";
    case nfp::NFPCycleRole::FeasiblePocket:
        return "FeasiblePocket";
    default:
        return "Unknown";
    }
}

void printNFP(const nfp::NFPResult& nfpResult) {
    std::cout << "NFP Result\n";
    std::cout << std::string(60, '=') << "\n\n";

    std::cout << "Cycles: " << nfpResult.cycles.size() << "\n";
    std::cout << std::string(60, '-') << "\n";

    for (size_t i = 0; i < nfpResult.cycles.size(); i++) {
        const auto& cycle = nfpResult.cycles[i];

        std::cout << "Cycle " << (i + 1)
            << " [" << cycleRoleToString(cycle.role) << "]"
            << " (" << cycle.points.size() << " points):\n";

        for (size_t j = 0; j < cycle.points.size(); j++) {
            const auto& p = cycle.points[j];

            std::cout << "  p" << j << ": ("
                << std::fixed << std::setprecision(4)
                << p.x << ", " << p.y << ")\n";
        }

        std::cout << "\n";
    }

    std::cout << "Slits: " << nfpResult.slits.size() << "\n";
    std::cout << std::string(60, '-') << "\n";

    for (size_t i = 0; i < nfpResult.slits.size(); i++) {
        const auto& s = nfpResult.slits[i];

        std::cout << "Slit " << (i + 1) << ":\n";
        std::cout << "  a: ("
            << std::fixed << std::setprecision(4)
            << s.a.x << ", " << s.a.y << ")\n";
        std::cout << "  b: ("
            << std::fixed << std::setprecision(4)
            << s.b.x << ", " << s.b.y << ")\n\n";
    }

    std::cout << "Isolated points: " << nfpResult.isolated_points.size() << "\n";
    std::cout << std::string(60, '-') << "\n";

    for (size_t i = 0; i < nfpResult.isolated_points.size(); i++) {
        const auto& p = nfpResult.isolated_points[i];

        std::cout << "Point " << (i + 1) << ": ("
            << std::fixed << std::setprecision(4)
            << p.x << ", " << p.y << ")\n";
    }

    std::cout << "\n";
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


    auto nfp = nfp::processNFP(poly_a, poly_b);

    printNFP(nfp);

    // Verify no self-intersection
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

    
    auto nfp = nfp::processNFP(poly_a, poly_b);

    printNFP(nfp);


}

void test4_anomalyFalseSlit() {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "TEST 4:square and simple\n";
    std::cout << std::string(70, '=') << "\n\n";

    // Polygon 1
    std::vector<Point> verts_a = {
    Point(0, 0),
    Point(3, 0),
    Point(3, 1),
    Point(5, 2),
    Point(5, 3),
    Point(6, 4),
    Point(4, 4),
    Point(4, 3),
    Point(2, 2),
    Point(2, 4),
    Point(0, 4),
    Point(1, 1)
    };

    // Polygon 2
    std::vector<Point> verts_b = {
        Point(0, 0),
        Point(2, 0),
        Point(2, 2),
        Point(0, 2)
    };

    Polygon poly_a(verts_a);
    Polygon poly_b(verts_b);

    printPolygon("Polygon A (square )", poly_a);
    printPolygon("Polygon B (simple)", poly_b);


    auto nfp = nfp::processNFP(poly_a, poly_b);

    printNFP(nfp);


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
//    auto nfp = nfp::generate(poly_a, poly_b);
//
//    printNFP(nfp);
//
//    // Verify no self-intersection
//    
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

    
    auto nfp = nfp::processNFP(poly_a, poly_b);
//

    printNFP(nfp);
//

    // Verify no self-intersection
    

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

    
    auto nfp = nfp::processNFP(poly_a, poly_b);
    //

    printNFP(nfp);
    //


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

    
    auto nfp = nfp::processNFP(poly_a, poly_b);
//

    printNFP(nfp);
//

    // Verify no self-intersection

}


void testSimple() {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "TEST Simple polygon with itself.\n";
    std::cout << std::string(70, '=') << "\n\n";

    // L-1
    std::vector<Point> verts_a = {
        Point(0, 0),
        Point(2, 2),
        Point(4, 0),
        Point(2, -2)
    };

    // L-2
    std::vector<Point> verts_b = {
        Point(0, 0),
        Point(2, 2),
        Point(4, 0),
        Point(2, -2)
    };

    Polygon poly_a(verts_a);
    Polygon poly_b(verts_b);

    printPolygon("Polygon A (L-shape)", poly_a);
    printPolygon("Polygon B (L-shape)", poly_b);

    
    auto nfp = nfp::processNFP(poly_a, poly_b);

    printNFP(nfp);

    // Verify no self-intersection

}

void Test_interlocking_point(){
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "TEST Simple polygon with interlocking point hole in NFP.\n";
    std::cout << std::string(70, '=') << "\n\n";
    //
    


    // L-1
    std::vector<Point> verts_a = {
        Point(0,0),
        Point(5,0),
        Point(5,4),
        Point(4,3),
        Point(4,8),
        Point(5,7),
        Point(5,10),
        Point(0,10)
    };

    // L-2
    std::vector<Point> verts_b = {
        Point(5,0),
        Point(10,0),
        Point(10,10),
        Point(5,10),
        Point(5,7),
        Point(4,8),
        Point(4,3),
        Point(5,4)
    };

    Polygon poly_a(verts_a);
    Polygon poly_b(verts_b);

    printPolygon("Polygon A (The die)", poly_a);
    printPolygon("Polygon B (the tap)", poly_b);

    
    auto nfp = nfp::processNFP(poly_a, poly_b);

    printNFP(nfp);

    // Verify no self-intersection


}

void Test_interlocking_line() {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "TEST Simple polygon with interlocking line hole in NFP.\n";
    std::cout << std::string(70, '=') << "\n\n";

    // L-1
    std::vector<Point> verts_a = {
        Point(0,0),
        Point(5,0),
        Point(5,4),
        Point(4,3),
        Point(2,3),
        Point(2,7),
        Point(4,7),
        Point(5,6),
        Point(5,10),
        Point(0,10)
    };

    // L-2
    std::vector<Point> verts_b = {
        Point(5,0),
        Point(10,0),
        Point(10,10),
        Point(5,10),
        Point(5,6),
        Point(3,6),
        Point(3,7),
        Point(2,7),
        Point(2,3),
        Point(3,3),
        Point(3,4),
        Point(5,4)
    };

    Polygon poly_a(verts_a);
    Polygon poly_b(verts_b);

    printPolygon("Polygon A (The die)", poly_a);
    printPolygon("Polygon B (The tap)", poly_b);

    
    auto nfp = nfp::processNFP(poly_a, poly_b);

    printNFP(nfp);

    // Verify no self-intersection

}

void Test_interlocking_hole() {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "TEST Simple polygon with interlocking area hole in NFP.\n";
    std::cout << std::string(70, '=') << "\n\n";

    // L-1
    std::vector<Point> verts_a = {
        Point(0,0),
        Point(5,0),
        Point(5,4),
        Point(4,3),
        Point(2,3),
        Point(2,7),
        Point(4,7),
        Point(5,6),
        Point(5,10),
        Point(0,10)
    };

    // L-2
    std::vector<Point> verts_b = {
        Point(5,0),
        Point(10,0),
        Point(10,10),
        Point(5,10),
        Point(5,5),
        Point(3,5),
        Point(3,6),
        Point(2,6),
        Point(2,3),
        Point(3,3),
        Point(3,4),
        Point(5,4)
    };

    Polygon poly_a(verts_a);
    Polygon poly_b(verts_b);

    printPolygon("Polygon A (the die)", poly_a);
    printPolygon("Polygon B (the tap)", poly_b);

    
    auto nfp = nfp::processNFP(poly_a, poly_b);

    printNFP(nfp);

    // Verify no self-intersection

}


int main() {
    std::cout << "NestingGraph NFP Test Program\n";
    Test_interlocking_point();
    Test_interlocking_line();
    Test_interlocking_hole();
    test4_anomalyFalseSlit();
    //test3_L_shape_L();
    //test4();
    //test5();
    //test6();
    //testSimple();
    return 0;
}