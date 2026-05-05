// nfp_bridge.h
#pragma once
#include <utility>
#include <vector>

#if defined(_WIN32)
  #if defined(NFP_BRIDGE_BUILD)
    #define NFP_API __declspec(dllexport)
  #else
    #define NFP_API __declspec(dllimport)
  #endif
#else
  #define NFP_API
#endif

namespace nfp {
using Point   = std::pair<double,double>;
using Polygon = std::vector<Point>;
using Polys   = std::vector<Polygon>;

enum class NFPCycleRole {
    OuterBoundary,
    HoleBoundary,
    FeasiblePocket
};

struct NFPPoint {
    double x;
    double y;
};

struct NFPSegment {
    NFPPoint a;
    NFPPoint b;
};

struct NFPCycle {
    NFPCycleRole role;
    std::vector<NFPPoint> points;
};

struct NFPResult {
    std::vector<NFPCycle> cycles;
    std::vector<NFPSegment> slits;
    std::vector<NFPPoint> isolated_points;
};

// Keep this frozen forever if the client must never touch NFP again:

NFP_API NFPResult processNFP(const Polygon& polyA, const Polygon& polyB);

}