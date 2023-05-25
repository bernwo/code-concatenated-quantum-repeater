#include "FiveQubitCode.hpp"

std::array<double, 32> FiveQubitCode::getLogicalKet(const double theta) {
    const double s = std::sin(0.5 * theta);
    const double c = std::cos(0.5 * theta);
    std::array<double, 32> output = {0.25 * c, -0.25 * s, -0.25 * s, -0.25 * c, -0.25 * s, 0.25 * c, -0.25 * c, -0.25 * s, -0.25 * s, 0.25 * c, 0.25 * c, 0.25 * s, -0.25 * c, 0.25 * s, -0.25 * s, -0.25 * c, -0.25 * s, -0.25 * c, 0.25 * c, -0.25 * s, 0.25 * c, 0.25 * s, 0.25 * s, -0.25 * c, -0.25 * c, -0.25 * s, 0.25 * s, -0.25 * c, -0.25 * s, -0.25 * c, -0.25 * c, 0.25 * s};
    return output;
}

std::array<double, FastPauli::CArrLen<5>> FiveQubitCode::getLogicalDensityMatrix(const double theta) {
    std::array<double, FastPauli::CArrLen<5>> cArrCodeword;
    const std::array<double, 32> ketL = FiveQubitCode::getLogicalKet(theta);
    FastPauli::ketToCArr(cArrCodeword.data(), ketL.data(), 5);
    return cArrCodeword;
}
