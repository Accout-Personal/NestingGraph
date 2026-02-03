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

// Keep this frozen forever if the client must never touch NFP again:

NFP_API Polys processNFP(const Polygon& polyA, const Polygon& polyB);

}