#include "AuxMath.hpp"
#include <cstdint>
#include <vector>
#include <cmath>
#include <iostream>
#include <cassert>

double AuxMath::nchoosek(const int n, const int k){ // Copied from https://stackoverflow.com/questions/55421835/c-binomial-coefficient-is-too-slow
    assert(n >= k);
    double tmp = n - k + 1;
    for (int i = 1; i < k; ++i) {
        tmp = tmp * (n - k + 1 + i) / (i + 1);
    }
    return tmp;
}

double AuxMath::f_(const double i, const double start, const double end, const int n){
    return start + i*(end - start)/(static_cast<double>(n) - 1);
}

std::vector<int> AuxMath::arange(const int start, const int end, const int d){ // Inclusive of 'end'.
    assert(d > 0);
    const int n = (end-start)/d + 1; // equivalent to std::ceil((end-start)/d).
    std::vector<int> out(n);
    for (int i = 0; i < n; i++){
        out[i] = start+i*d;
    }
    return out;
}

std::vector<double> AuxMath::linspace(const double start, const double end, const int n){
    std::vector<double> out(n);
    for (int i = 0; i < n; i++){
        out[i] = AuxMath::f_(i,start,end,n);
    }
    return out;
}

std::vector<double> AuxMath::logspace(const double start, const double end, const int n){
    std::vector<double> out(n);
    for (int i = 0; i < n; i++){
        out[i] = std::pow(10.0,AuxMath::f_(i,start,end,n));
    }
    return out;
}

template <typename T>
void AuxMath::print1DVec(std::vector<T> const &v){
    std::cout << "\033[0;34m";
    for (auto const &el: v){
        std::cout << el << " ";
    }
    std::cout << "\033[0m\n";
}
template void AuxMath::print1DVec(std::vector<uint16_t> const &v);
template void AuxMath::print1DVec(std::vector<int> const &v);
template void AuxMath::print1DVec(std::vector<double> const &v);