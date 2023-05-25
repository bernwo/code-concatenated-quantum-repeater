#ifndef _AUXMATHHPP_
#define _AUXMATHHPP_
#include <cstdint>
#include <vector>

namespace AuxMath{
    double nchoosek(const int n, const int k);
    std::vector<int> arange(const int start, const int end, const int d);
    std::vector<double> linspace(const double start, const double end, const int n);
    std::vector<double> logspace(const double start, const double end, const int n);
    double f_(const double i, const double start, const double end, const int n);
    template <typename T>
    void print1DVec(std::vector<T> const &v);
}

#endif