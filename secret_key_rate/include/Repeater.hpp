#ifndef _REPEATERHPP_
#define _REPEATERHPP_
#include <boost/math/special_functions/binomial.hpp>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <filesystem>
#ifdef __APPLE__
    #ifndef FMT_HEADER_ONLY
        #define FMT_HEADER_ONLY
    #endif
    #include <fmt/core.h>  // Use `homebrew` and run `brew install fmt` on your terminal to install fmt. We have to use `fmt` instead because `std::format` is not supported by Apple's Clang.
    #include <fmt/os.h>
using fmt::format;
#else
    #include <format>
using std::format;
#endif
#include <fstream>
#include <tuple>
#include <vector>

#include "AuxMath.hpp"
#include "BinaryData.hpp"

namespace Repeater {
struct RepeaterInfo {
    double reencodingError;     // Denoted as ϵ_r. It is the reencoding error rate.
    int reencodingErrorFactor;  // reencodingErrorFactor is such that ϵ_r = reencodingErrorFactor × ϵ_0, where ϵ_0 is the
                                // depolarising error of each single qubit in the generated tree. ϵ_0 is also the error
                                // rate of (non-teleported) spin-spin gates.
    int teleportedErrorFactor;  // teleportedErrorFactor is such that ϵ_tele = teleportedErrorFactor × ϵ_0, where ϵ_tele
                                // is the error rate of the teleported spin-spin gates.
    double Latt_val;            // Attenuation length.
    double eta_d_val;           // Overall efficiency.
    double tau_ph_val;          // Single photon emission time.
    double tau_ss_val;          // Spin-spin gate time.
    std::vector<double> SKR;
    std::vector<double> cost;
    std::vector<double> L0;
    std::vector<int> ntot;
    std::vector<double> Ltot;
    std::vector<std::vector<uint16_t>> treeVecs;
    std::vector<int> c;
    std::vector<double> pTransSum;
    std::vector<double> logicalError;
    std::vector<double> eta_e;
    std::vector<double> f_val;
    size_t size;
    size_t treeVecSize;
    void resize(const size_t input_size, const size_t eachTreeVecSize) {
        this->SKR.resize(input_size);
        this->cost.resize(input_size);
        this->L0.resize(input_size);
        this->ntot.resize(input_size);
        this->Ltot.resize(input_size);
        this->treeVecs.resize(input_size);
        for (size_t i = 0; i < input_size; i++) {
            this->treeVecs[i].resize(eachTreeVecSize);
        }
        this->c.resize(input_size);
        this->pTransSum.resize(input_size);
        this->logicalError.resize(input_size);
        this->eta_e.resize(input_size);
        this->f_val.resize(input_size);
        this->size = input_size;
        this->treeVecSize = eachTreeVecSize;
    }
    void save(std::string const &filename) const {
        std::ofstream outfile_noappend;
        std::ofstream outfile_append;
        if (!this->SKR.empty() && !this->cost.empty() && !this->L0.empty() && !this->Ltot.empty() &&
            !this->treeVecs.empty() && !this->c.empty()) {
            if (this->doAppend(filename)) {
                outfile_append.open(filename, std::ios::binary | std::ios::app);
                outfile_append.write((char *)&this->reencodingError, sizeof(double));
                outfile_append.write((char *)&this->reencodingErrorFactor, sizeof(int));
                outfile_append.write((char *)&this->teleportedErrorFactor, sizeof(int));
                outfile_append.write((char *)&this->Latt_val, sizeof(double));
                outfile_append.write((char *)&this->eta_d_val, sizeof(double));
                outfile_append.write((char *)&this->tau_ph_val, sizeof(double));
                outfile_append.write((char *)&this->tau_ss_val, sizeof(double));
                outfile_append.write((char *)&this->size, sizeof(size_t));
                outfile_append.write((char *)&this->treeVecSize, sizeof(size_t));
                outfile_append.write((char *)this->SKR.data(), this->SKR.size() * sizeof(double));
                outfile_append.write((char *)this->cost.data(), this->cost.size() * sizeof(double));
                outfile_append.write((char *)this->L0.data(), this->L0.size() * sizeof(double));
                outfile_append.write((char *)this->ntot.data(), this->ntot.size() * sizeof(int));
                outfile_append.write((char *)this->Ltot.data(), this->Ltot.size() * sizeof(double));
                for (size_t i = 0; i < this->treeVecs.size(); i++) {
                    outfile_append.write((char *)this->treeVecs[i].data(), this->treeVecs[i].size() * sizeof(uint16_t));
                }
                outfile_append.write((char *)this->c.data(), this->c.size() * sizeof(int));
                outfile_append.write((char *)this->pTransSum.data(), this->pTransSum.size() * sizeof(double));
                outfile_append.write((char *)this->logicalError.data(), this->logicalError.size() * sizeof(double));
                outfile_append.write((char *)this->eta_e.data(), this->eta_e.size() * sizeof(double));
                outfile_append.write((char *)this->f_val.data(), this->f_val.size() * sizeof(double));
                outfile_append.close();
                std::cout << format(
                                 "\x1b[0;92mRepeater::RepeaterInfo - Data appended at \x1b[0m\x1b[1;4;92m{}!\x1b[0m",
                                 filename)
                          << std::endl;
            } else if (!std::filesystem::exists(filename)) {
                outfile_noappend.open(filename, std::ios::binary);
                outfile_noappend.write((char *)&this->reencodingError, sizeof(double));
                outfile_noappend.close();

                outfile_append.open(filename, std::ios::binary | std::ios::app);
                outfile_append.write((char *)&this->reencodingErrorFactor, sizeof(int));
                outfile_append.write((char *)&this->teleportedErrorFactor, sizeof(int));
                outfile_append.write((char *)&this->Latt_val, sizeof(double));
                outfile_append.write((char *)&this->eta_d_val, sizeof(double));
                outfile_append.write((char *)&this->tau_ph_val, sizeof(double));
                outfile_append.write((char *)&this->tau_ss_val, sizeof(double));
                outfile_append.write((char *)&this->size, sizeof(size_t));
                outfile_append.write((char *)&this->treeVecSize, sizeof(size_t));
                outfile_append.write((char *)this->SKR.data(), this->SKR.size() * sizeof(double));
                outfile_append.write((char *)this->cost.data(), this->cost.size() * sizeof(double));
                outfile_append.write((char *)this->L0.data(), this->L0.size() * sizeof(double));
                outfile_append.write((char *)this->ntot.data(), this->ntot.size() * sizeof(int));
                outfile_append.write((char *)this->Ltot.data(), this->Ltot.size() * sizeof(double));
                for (size_t i = 0; i < this->treeVecs.size(); i++) {
                    outfile_append.write((char *)this->treeVecs[i].data(), this->treeVecs[i].size() * sizeof(uint16_t));
                }
                outfile_append.write((char *)this->c.data(), this->c.size() * sizeof(int));
                outfile_append.write((char *)this->pTransSum.data(), this->pTransSum.size() * sizeof(double));
                outfile_append.write((char *)this->logicalError.data(), this->logicalError.size() * sizeof(double));
                outfile_append.write((char *)this->eta_e.data(), this->eta_e.size() * sizeof(double));
                outfile_append.write((char *)this->f_val.data(), this->f_val.size() * sizeof(double));
                outfile_append.close();
                std::cout << format(
                                 "\x1b[0;92mRepeater::RepeaterInfo - File written at \x1b[0m\x1b[1;4;92m{}!\x1b[0m",
                                 filename)
                          << std::endl;
            } else {
                std::cout << format(
                                 "\x1b[0;93mRepeater::RepeaterInfo - File not written because the data for "
                                 "this particular error rate exist already in \x1b[0m\x1b[1;4;93m{}!\x1b[0m",
                                 filename)
                          << std::endl;
            }
        } else {
            std::cerr << "\x1b[0;91mRepeater::RepeaterInfo - File not written because at least one of the "
                         "std::vector's are empty!\x1b[0m\n";
        }
    }
    void load(std::string const &filename) {
        std::ifstream inputfile(filename, std::ios::binary);
        inputfile.seekg(0, std::ios::beg);
        inputfile.read(reinterpret_cast<char *>(&this->reencodingError), sizeof(double));
        inputfile.read(reinterpret_cast<char *>(&this->reencodingErrorFactor), sizeof(int));
        inputfile.read(reinterpret_cast<char *>(&this->teleportedErrorFactor), sizeof(int));
        inputfile.read(reinterpret_cast<char *>(&this->Latt_val), sizeof(double));
        inputfile.read(reinterpret_cast<char *>(&this->eta_d_val), sizeof(double));
        inputfile.read(reinterpret_cast<char *>(&this->tau_ph_val), sizeof(double));
        inputfile.read(reinterpret_cast<char *>(&this->tau_ss_val), sizeof(double));
        inputfile.read(reinterpret_cast<char *>(&this->size), sizeof(size_t));
        inputfile.read(reinterpret_cast<char *>(&this->treeVecSize), sizeof(size_t));
        resize(this->size, this->treeVecSize);
        inputfile.read(reinterpret_cast<char *>(&this->SKR[0]), this->SKR.size() * sizeof(double));
        inputfile.read(reinterpret_cast<char *>(&this->cost[0]), this->cost.size() * sizeof(double));
        inputfile.read(reinterpret_cast<char *>(&this->L0[0]), this->L0.size() * sizeof(double));
        inputfile.read(reinterpret_cast<char *>(&this->ntot[0]), this->ntot.size() * sizeof(int));
        inputfile.read(reinterpret_cast<char *>(&this->Ltot[0]), this->Ltot.size() * sizeof(double));
        for (size_t i = 0; i < this->treeVecs.size(); i++) {
            inputfile.read(reinterpret_cast<char *>(&this->treeVecs[i][0]),
                           this->treeVecs[i].size() * sizeof(uint16_t));
        }
        inputfile.read(reinterpret_cast<char *>(&this->c[0]), this->c.size() * sizeof(int));
        inputfile.read(reinterpret_cast<char *>(&this->pTransSum[0]), this->pTransSum.size() * sizeof(double));
        inputfile.read(reinterpret_cast<char *>(&this->logicalError[0]), this->logicalError.size() * sizeof(double));
        inputfile.read(reinterpret_cast<char *>(&this->eta_e[0]), this->eta_e.size() * sizeof(double));
        inputfile.read(reinterpret_cast<char *>(&this->f_val[0]), this->f_val.size() * sizeof(double));
        inputfile.close();
        std::cout << format("\x1b[0;92mRepeater::RepeaterInfo - File loaded from \x1b[0m\x1b[1;4;92m{}!\x1b[0m",
                            filename)
                  << std::endl;
    }
    bool doSKRCalculation(std::string const &filename) const {
        if (!std::filesystem::exists(filename) || this->doAppend(filename)) {
            return true;
        } else {
            return false;
        }
    }

   private:
    bool doAppend(std::string const &filename) const {
        double reencodingError_tmp;
        size_t L_points_tmp;
        size_t treeVecSize_tmp;
        std::ifstream file_input;
        bool notfoundMatch = true;
        if (std::filesystem::exists(filename)) {
            file_input.open(filename, std::ios::binary);

            while (file_input.peek() != EOF) {
                file_input.read(reinterpret_cast<char *>(&reencodingError_tmp), sizeof(double));
                // std::cout << "reencodingError: " << reencodingError_tmp << std::endl;
                if (reencodingError_tmp == this->reencodingError) {
                    notfoundMatch = false;
                    break;
                }
                file_input.ignore(4 * sizeof(double) + 2 * sizeof(int));
                file_input.read(reinterpret_cast<char *>(&L_points_tmp), sizeof(size_t));
                file_input.read(reinterpret_cast<char *>(&treeVecSize_tmp), sizeof(size_t));
                file_input.ignore(8 * L_points_tmp * sizeof(double) + 2 * L_points_tmp * sizeof(int) +
                                  L_points_tmp * treeVecSize_tmp * sizeof(uint16_t));
            }
            // std::cout << "while loop ended." << std::endl;
        } else {
            std::cerr << "\x1b[0;91mError: File not found!\x1b[0m\n";
        }
        return notfoundMatch;
    }
};

std::vector<std::pair<double, int>> getL0AndntotArray(const double Ltot, const double L0_min, const double L0_max,
                                                      const int numPerArray, const bool verbose);

/**
 * @brief Computes the generation time τ_tree of each given tree vector and returns it in units of seconds.
 *
 * @param tau_ph The emission time for a single photon in units of seconds.
 * @param tau_ss The spin-spin (non-teleported) gate time in units of seconds.
 * @param treeVecs The vector of tree vectors.
 * @return std::vector<double>
 */
std::vector<double> getTauTreeArray(const double tau_ph, const double tau_ss,
                                    std::vector<std::vector<uint16_t>> const &treeVecs);

/**
 * @brief Compute the binary entropy function h(x) = -x⋅log₂(x)-(1-x)⋅log₂(1-x).
 *
 * @param x A real number between 0 and 1, inclusive.
 * @return double
 */
double h(const double x);

/**
 * @brief Computes the secret key fraction for the asymptotic six-state variant of the BB84 protocol.
 *
 * @param error The effective logical error rate at the end of the repeater network chain.
 * @return double
 */
double f(const double error);

double eta(const double L0, const double Latt);

double mu(const double L0, const double Latt, const double eta_d);

double eta_e(const double L0, const double Latt, const double eta_d, std::vector<uint16_t> const &treeVec);

/**
 * @brief Computes the tree generation time τ_tree of a given tree vector and returns it in units of seconds.
 *
 * @param tau_ph The emission time for a single photon in units of seconds.
 * @param tau_ss The spin-spin (non-teleported) gate time in units of seconds.
 * @param treeVec The vector of tree vectors.
 * @return double
 */
double tauTree(const double tau_ph, const double tau_ss, std::vector<uint16_t> const &treeVec);

double pTrans(const double eta_e_val, const double n, const int c, const int i);

std::tuple<double, double, double, double, double, double>
SKR(BinaryData::ECFlagInfo const &ecFlagInfo, BinaryData::ECErasureInfo const &ecErasureInfo, const double Ltot,
    const int c, const int cMax, const double L0, const double Latt, const int ntot, const double eta_d,
    const double tau_ph, const double tau_ss, const double tau_tele, const double tau_meas, const double tauTree,
    std::vector<uint16_t> const &treeVec, const double kappa, const bool verbose);

std::tuple<double, double, double, double> SKRHomogeneous(const double Ltot, const double L0, const double Latt,
                                                          const int ntot, const double eta_d, const double tau_ph,
                                                          const double tau_ss, const double tauTree,
                                                          std::vector<uint16_t> const &treeVec,
                                                          const double errReencoding);

}  // namespace Repeater

#endif