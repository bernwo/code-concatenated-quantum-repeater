#ifndef _TREEHPP_
#define _TREEHPP_

#include <algorithm>
#include <cassert>
#include <cmath> // Provides std::floor
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

class Tree {
  private:
    bool verbose;
    uint16_t aMin;
    uint16_t bMin;
    uint16_t cMin;
    uint16_t NMax;
    uint16_t aMax;
    uint16_t cMax;
    uint32_t vecSize;
    std::vector<std::vector<uint16_t>> treeVec;
    void calcVecSize_();
    void calcVec_();

  public:
    Tree() = default;
    Tree(const uint16_t NMax, const uint16_t aMin, const uint16_t bMin, const uint16_t cMin, const int verbosity);
    uint32_t getSize();
    std::vector<std::vector<uint16_t>> getVec();
    void printVec();
    static uint16_t getPhotonNumber(std::vector<uint16_t> const &treeVec);
    static std::vector<std::vector<uint16_t>> getSortedTreeVecs(const double L0, const double Latt, const double eta_d,
                                                                std::vector<std::vector<uint16_t>> const &treeVecs,
                                                                int length);
};

#endif