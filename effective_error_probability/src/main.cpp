#include <array>
#include <chrono>  // https://stackoverflow.com/questions/22387586/measuring-execution-time-of-a-function-in-c
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>  // Provides std::ofstream
#include <iomanip>
#include <iostream>
#include <numbers>
#include <sstream>
#include <vector>

#include "AuxMath.hpp"
#include "FiveQubitCode.hpp"
#include "fpauli.hpp"
#include "npauli.hpp"

#ifndef NUMPOINTS
    #define NUMPOINTS 10000
#endif
#define NUMTHREADS 16  // Number of CPU threads

#define STRCOMPARE(cstring1, cstring2) static_cast<std::string>(cstring1) == cstring2

size_t getRounds(int x, const char* mode);
std::string datetime();
std::string to_string_with_precision(const double a_value, const int n);

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Too few arguments! Exiting!\n";
        std::cerr << "Usage: main.exe [mode (string)] [reencoding error (float)]\n";
        std::cerr << "\twhere mode can be '1erasure' or 'flag'\n";
        return EXIT_FAILURE;
    }
    const char* mode = argv[1];  // mode ∈ {"1erasure", "flag"
    std::istringstream os(argv[2]);
    double err;
    os >> err;

    if (!(STRCOMPARE(mode, "1erasure") || STRCOMPARE(mode, "flag"))) {
        std::cerr << "Usage: main.exe [mode (string)] [reencoding error (float)]\n";
        std::cerr << "\twhere mode can be '1erasure' or 'flag'\n";
        return EXIT_FAILURE;
    }
    if (argc > 3) {
        std::cerr << "Warning: Extra arguments ignored\n";
    }

    const auto t1 = std::chrono::high_resolution_clock::now();
    constexpr double errTeleFactor = 3;        // errTeleport = errTeleFactor * eps0.
    constexpr double errReencodingFactor = 3;  // eps0 = errReencoding / errReencodingFactor with eps0 being the error of non-teleported 2-qubit gates.
    constexpr double theta = std::numbers::pi / 2;
    auto nLinksArray = AuxMath::linspace(1, 10000, NUMPOINTS);
    const size_t rounds = getRounds(1, mode);

    std::cout << "  \x1b[0;94mmode: " << mode << "\x1b[0m" << std::endl;
    std::cout << "  \x1b[0;94mNUMPOINTS: " << NUMPOINTS << "\x1b[0m" << std::endl;
    std::cout << "  \x1b[0;94merr: " << to_string_with_precision(err, 0) << "\x1b[0m" << std::endl;
    std::cout << "  \x1b[0;94merrTeleFactor: " << errTeleFactor << "\x1b[0m" << std::endl;
    std::cout << "  \x1b[0;94merrReencodingFactor: " << errReencodingFactor << "\x1b[0m" << std::endl;
    std::cout << "  \x1b[0;94mtheta: " << theta << "\x1b[0m" << std::endl;
    std::cout << "  \x1b[0;94mrounds: " << rounds << "\x1b[0m" << std::endl;

    const std::string dataPath = "./binary_data/";
    std::filesystem::create_directory(dataPath);
    const std::string dataSpecificPath = dataPath + mode + (std::string) "_err" + to_string_with_precision(err, 0) + (std::string) "_reencodingFactor" + std::to_string((int)errReencodingFactor) + (std::string) "_teleFactor" + std::to_string((int)errTeleFactor) + (std::string) "/";
    std::filesystem::create_directory(dataSpecificPath);

    FiveQubitCode::ECInfo* ecInfo = new FiveQubitCode::ECInfo[rounds];
    for (size_t i = 0; i < rounds; i++) {
        ecInfo[i].resize(NUMPOINTS);
    }

    constexpr uint64_t cArrLen5 = FastPauli::CArrLen<5>;
    std::array<double, 528> cArrCodeword = FiveQubitCode::getLogicalDensityMatrix(theta);
    std::array<double, 32> logicalKet = FiveQubitCode::getLogicalKet(theta);
    std::array<double, 528> tempCArrCodeword;

    double etrans;

    std::array<double, 528> correctedDM;

    int i, j;
#pragma omp parallel shared(cArrCodeword, err, nLinksArray, ecInfo, rounds, logicalKet, errTeleFactor) private(etrans, j, correctedDM, i, tempCArrCodeword)  // Reference: https://stackoverflow.com/questions/2352895/how-to-ensure-a-dynamically-allocated-array-is-private-in-openmp
    {
#pragma omp for schedule(static, NUMPOINTS / NUMTHREADS)
        for (i = 0; i < nLinksArray.size(); i++) {
            tempCArrCodeword = cArrCodeword;
            for (j = 0; j < rounds; j++) {
                etrans = 1 - std::pow(1 - err, nLinksArray[i]) * (1 - err / errReencodingFactor);
                NoisyPauli::SingleQubitDepolarise_r(tempCArrCodeword.data(), 5, 0, etrans);
                NoisyPauli::SingleQubitDepolarise_r(tempCArrCodeword.data(), 5, 1, etrans);
                NoisyPauli::SingleQubitDepolarise_r(tempCArrCodeword.data(), 5, 2, etrans);
                NoisyPauli::SingleQubitDepolarise_r(tempCArrCodeword.data(), 5, 3, etrans);
                NoisyPauli::SingleQubitDepolarise_r(tempCArrCodeword.data(), 5, 4, etrans);

                if (STRCOMPARE(mode, "flag")) {
                    correctedDM = FiveQubitCode::OneRoundSequential_flag(ecInfo[j],
                                                                         logicalKet.data(),
                                                                         tempCArrCodeword.data(),
                                                                         err / errReencodingFactor,
                                                                         errTeleFactor * err / errReencodingFactor,
                                                                         errTeleFactor * err / errReencodingFactor,
                                                                         i);
                } else {
                    correctedDM = FiveQubitCode::OneRoundSequential(tempCArrCodeword.data(),
                                                                    err / errReencodingFactor,
                                                                    errTeleFactor * err / errReencodingFactor,
                                                                    FiveQubitCode::correctionMode::oneLoss,
                                                                    0);
                    ecInfo[j].grandErrors[i] = (1 - FastPauli::fidel_r(logicalKet.data(), correctedDM.data(), 5));
                    ecInfo[j].grandTraces[i] = (FastPauli::trace_r(correctedDM.data(), 5));
                }
                tempCArrCodeword = correctedDM;
            }
        }
    }

    const auto instantaneous_datetime = datetime();
    for (size_t i = 0; i < rounds; i++) {
        ecInfo[i].save(dataSpecificPath + (std::string) "ecInfo" + std::to_string(i) + (std::string) "_" + instantaneous_datetime + (std::string) ".dat");
    }

    std::ofstream outfile_nLinksArray((dataSpecificPath + (std::string) "nLinksArray_" + instantaneous_datetime + (std::string) ".dat"), std::ios::binary);  // This binary file can be read in Wolfram Mathematica by writing BinaryReadList["err.dat", "Real64"] and then you can proceed to plot it.
    outfile_nLinksArray.write((char*)nLinksArray.data(), NUMPOINTS * 8);
    outfile_nLinksArray.close();

    delete[] ecInfo;

    const auto t2 = std::chrono::high_resolution_clock::now();
    const auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "\n\033[1;34mElapsed time: \033[0m\033[1;4;33m" << ms_int.count() << "ms\033[0m\n";
    std::cout << "instantaneous_datetime: " << instantaneous_datetime << std::endl;
    std::cout << '\a' << std::endl;
    return EXIT_SUCCESS;
}

std::string datetime() {
    time_t rawtime;
    struct tm* timeinfo;
    char buffer[80];

    time(&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(buffer, 80, "%Y-%m-%d %H-%M-%S", timeinfo);
    return std::string(buffer);
}

std::string to_string_with_precision(const double a_value, const int n) {
    // Reference: https://stackoverflow.com/questions/16605967/set-precision-of-stdto-string-when-converting-floating-point-values
    std::ostringstream out;
    out.precision(n);
    out << std::scientific << a_value;
    return out.str();
}

size_t getRounds(int x, const char* mode) {
    return STRCOMPARE(mode, "flag") ? std::max(static_cast<size_t>(x), static_cast<size_t>(2)) : static_cast<size_t>(x);
}