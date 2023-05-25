#ifndef _BINARYDATAHPP_
#define _BINARYDATAHPP_

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

namespace BinaryData {
struct ECFlagInfo {
    std::array<std::vector<double>, 5> traces;  // These act like probabilities.
    std::array<double, 10000> errs0;
    std::array<double, 10000> errs1;
    void read(std::string const &pathFile0, std::string const &pathFile1) {
        std::ifstream file0(pathFile0, std::ios::binary);
        // file0.seekg(0, std::ios::end);
        // auto size0 = file0.tellg();
        file0.seekg(0, std::ios::beg);
        file0.ignore(800000);  // Skip 800000 bytes. Reference: https://stackoverflow.com/questions/18640001/read-several-bytes-jump-over-n-bytes-and-then-read-several-bytes-again-how
        // std::cout << "Size of file: " << size0 << " bytes." << std::endl;
        file0.read(reinterpret_cast<char *>(&errs0[0]), errs0.size() * sizeof(double));  // Reference: https://stackoverflow.com/questions/15143670/how-can-i-use-fread-on-a-binary-file-to-read-the-data-into-a-stdvector
        file0.close();

        std::ifstream file1(pathFile1, std::ios::binary);
        // file0.seekg(0, std::ios::end);
        // auto size0 = file0.tellg();
        file1.seekg(0, std::ios::beg);
        file1.ignore(400000);  // Skip 400000 bytes. Reference: https://stackoverflow.com/questions/18640001/
        traces[0].resize(10000);
        file1.read(reinterpret_cast<char *>(&traces[0][0]), traces[0].size() * sizeof(double));
        traces[1].resize(10000);
        file1.read(reinterpret_cast<char *>(&traces[1][0]), traces[1].size() * sizeof(double));
        traces[2].resize(10000);
        file1.read(reinterpret_cast<char *>(&traces[2][0]), traces[2].size() * sizeof(double));
        traces[3].resize(10000);
        file1.read(reinterpret_cast<char *>(&traces[3][0]), traces[3].size() * sizeof(double));
        traces[4].resize(10000);
        file1.read(reinterpret_cast<char *>(&traces[4][0]), traces[4].size() * sizeof(double));
        // file1.ignore(800000);  // Skip 800000 bytes. Reference: https://stackoverflow.com/questions/18640001/read-several-bytes-jump-over-n-bytes-and-then-read-several-bytes-again-how
        // std::cout << "Size of file: " << size0 << " bytes." << std::endl;
        file1.read(reinterpret_cast<char *>(&errs1[0]), errs1.size() * sizeof(double));  // Reference: https://stackoverflow.com/questions/15143670/how-can-i-use-fread-on-a-binary-file-to-read-the-data-into-a-stdvector
        file1.close();
    }
    void println(std::array<double, 10000> const &array) const {
        for (size_t i = 0; i < 10000; i++) {
            std::cout << std::setprecision(std::numeric_limits<double>::digits10) << array[i] << '\n';
        }
    }
    double getRecursiveFidelity(const size_t n, const uint64_t c) const {
        double denom;
        const double fidelity0 = 1 - errs0[n];
        const double fidelity1 = 1 - errs1[n];
        const double zeta = std::sqrt(4 * fidelity1 - 3 * (fidelity0 * fidelity0));
        switch (c) {
            case UINT64_C(1):
                return fidelity0;
                break;
            case UINT64_C(2):
                return fidelity1;
                break;
            default:
                denom = std::pow(2.0, c + 1);
                if (!std::isnan(denom)) {
                    return (std::pow(fidelity0 + zeta, c + 1) - std::pow(fidelity0 - zeta, c + 1)) / zeta / denom;
                } else {
                    return 0.0;
                }
                break;
        }
    }
    void clear() {
        for (size_t i = 0; i < 5; i++) {
            if (!traces[i].empty()) {
                traces[i].clear();
            }
        }
        errs0.fill(0);
        errs1.fill(0);
    }
};
struct ECErasureInfo {
    std::array<double, 10000> errs;
    void read(std::string const &pathFile) {
        std::ifstream input_file;
        if (std::filesystem::exists(pathFile)) {
            input_file.open(pathFile, std::ios::binary);
            // file0.seekg(0, std::ios::end);
            // auto size0 = file0.tellg();
            input_file.seekg(0, std::ios::beg);
            input_file.ignore(800000);  // Skip 800000 bytes. Reference: https://stackoverflow.com/questions/18640001/read-several-bytes-jump-over-n-bytes-and-then-read-several-bytes-again-how
            // std::cout << "Size of file: " << size0 << " bytes." << std::endl;
            input_file.read(reinterpret_cast<char *>(&this->errs[0]), this->errs.size() * sizeof(double));  // Reference: https://stackoverflow.com/questions/15143670/how-can-i-use-fread-on-a-binary-file-to-read-the-data-into-a-stdvector
            input_file.close();
        } else {
            std::cerr << "\x1b[0;91mFile not found!\n";
        }
    }
};

}  // namespace BinaryData

#endif