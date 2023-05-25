#ifndef _FIVEQUBITCODEHPP_
#define _FIVEQUBITCODEHPP_

#include <algorithm>  // Provide std::transform. Reference: https://stackoverflow.com/questions/3376124/how-to-add-element-by-element-of-two-stl-vectors
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>  // Provides definition for uint64_t.
#include <fstream>
#include <functional>
#include <iostream>
#include <numbers>
#include <string>
#include <tuple>
#include <vector>

#include "fpauli.hpp"
#include "npauli.hpp"

class FiveQubitCode {
   private:
    /**
     * @brief Applies correction to the input density matrix according to the syndrome `i`. Not fault-tolerant.
     *
     * @param i i ∈ {0,1,2,...,15} is all the possibilities of the ancilla qubits on/off combinations, i.e., the syndromes.
     * @param cArr5 The real density matrix in the form of a C-style array of length `pow(2,n-1)+pow(2,2n-1)` where `n=5` is the number of qubits.
     */
    static void applyCorrection(const uint16_t i, double *cArr5);  // For non-fault tolerant implementation of the [[5,1,3]] code, i.e., no flag qubits.

    /**
     * @brief Applies correction to the input density matrix according to the syndrome `i` and which `flag ∈ {0,1,2,3}` was triggered. Fault-tolerant. See Bernard's MSc thesis.
     *
     * @param i i ∈ {0,1,2,...,15} is all the possibilities of the ancilla qubits on/off combinations, i.e., the syndromes.
     * @param flag flag ∈ {0,1,2,3} specifies which flag was triggered.
     * @param cArr5 The real density matrix in the form of a C-style array of length `pow(2,n-1)+pow(2,2n-1)` where `n=5` is the number of qubits.
     */
    static void applyCorrection_flag(const uint16_t i, const uint16_t flag, double *cArr5);

    /**
     * @brief Applies 1-erasure error correction to the input density matrix according to the syndrome `i`. Not fault-tolerant by nature.
     *
     * @param i i ∈ {0,1,2,...,15} is all the possibilities of the ancilla qubits on/off combinations, i.e., the syndromes.
     * @param cArr5 The real density matrix in the form of a C-style array of length `pow(2,n-1)+pow(2,2n-1)` where `n=5` is the number of qubits.
     * @param k i ∈ {0,1,2,3,4} is qubit lost index.
     */
    static void applyCorrection_1erasure(const uint16_t i, double *cArr5, const int k);

   public:
    enum correctionMode { noLossW1,
                          noLossW2,
                          twoLoss,
                          oneLoss };
    struct ECInfo {
        std::array<std::vector<double>, 5> errors;  // 1 - fidelity
        std::array<std::vector<double>, 5> traces;  // These act like probabilities.
        std::vector<double> grandErrors;
        std::vector<double> grandTraces;
        void resize(size_t numElems) {
            for (int i = 0; i < 5; i++) {
                errors[i].resize(numElems);
                traces[i].resize(numElems);
            }
            grandErrors.resize(numElems);
            grandTraces.resize(numElems);
        }
        void save(std::string const &filename) {
            std::ofstream outfile_noappend;
            std::ofstream outfile_append;
            outfile_noappend.open(filename, std::ios::binary);
            outfile_noappend.write((char *)errors[0].data(), errors[0].size() * 8);
            outfile_noappend.close();

            outfile_append.open(filename, std::ios::binary | std::ios::app);
            for (int i = 1; i < 5; i++) {
                outfile_append.write((char *)errors[i].data(), errors[i].size() * 8);
            }
            for (int i = 0; i < 5; i++) {
                outfile_append.write((char *)traces[i].data(), traces[i].size() * 8);
            }
            outfile_append.write((char *)grandErrors.data(), grandErrors.size() * 8);
            outfile_append.write((char *)grandTraces.data(), grandTraces.size() * 8);
            outfile_append.close();
        }
        void print() {
            std::cout << "\x1b[0;95mecInfo.grandTraces = {";
            for (size_t i = 0; i < grandTraces.size(); i++) {
                if (i == grandTraces.size() - 1) {
                    std::cout << grandTraces[i];
                } else {
                    std::cout << grandTraces[i] << ", ";
                }
            }
            std::cout << "};\x1b[0m" << std::endl;
            std::cout << "\x1b[0;92mecInfo.traces[0] = {";
            for (size_t i = 0; i < traces[0].size(); i++) {
                if (i == traces[0].size() - 1) {
                    std::cout << traces[0][i];
                } else {
                    std::cout << traces[0][i] << ", ";
                }
            }
            std::cout << "};\x1b[0m" << std::endl;
            std::cout << "\x1b[0;93mecInfo.traces[1] = {";
            for (size_t i = 0; i < traces[1].size(); i++) {
                if (i == traces[1].size() - 1) {
                    std::cout << traces[1][i];
                } else {
                    std::cout << traces[1][i] << ", ";
                }
            }
            std::cout << "};\x1b[0m" << std::endl;
            std::cout << "\x1b[0;95mecInfo.traces[2] = {";
            for (size_t i = 0; i < traces[2].size(); i++) {
                if (i == traces[2].size() - 1) {
                    std::cout << traces[2][i];
                } else {
                    std::cout << traces[2][i] << ", ";
                }
            }
            std::cout << "};\x1b[0m" << std::endl;
            std::cout << "\x1b[0;96mecInfo.traces[3] = {";
            for (size_t i = 0; i < traces[3].size(); i++) {
                if (i == traces[3].size() - 1) {
                    std::cout << traces[3][i];
                } else {
                    std::cout << traces[3][i] << ", ";
                }
            }
            std::cout << "};\x1b[0m" << std::endl;
            std::cout << "\x1b[0;91mecInfo.traces[4] = {";
            for (size_t i = 0; i < traces[4].size(); i++) {
                if (i == traces[4].size() - 1) {
                    std::cout << traces[4][i];
                } else {
                    std::cout << traces[4][i] << ", ";
                }
            }
            std::cout << "};\x1b[0m" << std::endl;
            std::cout << "====================================" << std::endl;
            std::cout << "\x1b[0;94mecInfo.grandErrors = {";
            for (size_t i = 0; i < grandErrors.size(); i++) {
                if (i == grandErrors.size() - 1) {
                    std::cout << grandErrors[i];
                } else {
                    std::cout << grandErrors[i] << ", ";
                }
            }
            std::cout << "};\x1b[0m" << std::endl;

            std::cout << "\x1b[0;92mecInfo.errors[0] = {";
            for (size_t i = 0; i < errors[0].size(); i++) {
                if (i == errors[0].size() - 1) {
                    std::cout << errors[0][i];
                } else {
                    std::cout << errors[0][i] << ", ";
                }
            }
            std::cout << "};\x1b[0m" << std::endl;
            std::cout << "\x1b[0;93mecInfo.errors[1] = {";
            for (size_t i = 0; i < errors[1].size(); i++) {
                if (i == errors[1].size() - 1) {
                    std::cout << errors[1][i];
                } else {
                    std::cout << errors[1][i] << ", ";
                }
            }
            std::cout << "};\x1b[0m" << std::endl;
            std::cout << "\x1b[0;95mecInfo.errors[2] = {";
            for (size_t i = 0; i < errors[2].size(); i++) {
                if (i == errors[2].size() - 1) {
                    std::cout << errors[2][i];
                } else {
                    std::cout << errors[2][i] << ", ";
                }
            }
            std::cout << "};\x1b[0m" << std::endl;
            std::cout << "\x1b[0;96mecInfo.errors[3] = {";
            for (size_t i = 0; i < errors[3].size(); i++) {
                if (i == errors[3].size() - 1) {
                    std::cout << errors[3][i];
                } else {
                    std::cout << errors[3][i] << ", ";
                }
            }
            std::cout << "};\x1b[0m" << std::endl;
            std::cout << "\x1b[0;91mecInfo.errors[4] = {";
            for (size_t i = 0; i < errors[4].size(); i++) {
                if (i == errors[4].size() - 1) {
                    std::cout << errors[4][i];
                } else {
                    std::cout << errors[4][i] << ", ";
                }
            }
            std::cout << "};\x1b[0m" << std::endl;
        }
    };

    /**
     * @brief The first set of stabiliser operations corresponding to `XZZXI`. Fault-tolerant.
     *
     * @param cArray The real density matrix in the form of a C-style array of length `pow(2,n-1)+pow(2,2n-1)` where `n=7` is the number of qubits (including the ancilla and flag qubit).
     * @param err The error rate of the two-qubit gates in the stabiliser operation between the ancilla and the data qubits, i.e., the CNOT and CPHASE gates.
     * @param errTeleport The error rate of the teleported two-qubit gates in 3 of the 4 two-qubit gates. This is typically higher than `err`.
     * @param ancillaFlagErr The error rate of the two-qubit gates between the ancilla and the flag qubit. Usually same as `err`.
     */
    static void S1_flag(double *cArray, const double err, const double errTeleport, const double ancillaFlagErr);

    /**
     * @brief The first set of stabiliser operations corresponding to `IXZZX`. Fault-tolerant.
     *
     * @param cArray The real density matrix in the form of a C-style array of length `pow(2,n-1)+pow(2,2n-1)` where `n=7` is the number of qubits (including the ancilla and flag qubit).
     * @param err The error rate of the two-qubit gates in the stabiliser operation between the ancilla and the data qubits, i.e., the CNOT and CPHASE gates.
     * @param errTeleport The error rate of the teleported two-qubit gates in 3 of the 4 two-qubit gates. This is typically higher than `err`.
     * @param ancillaFlagErr The error rate of the two-qubit gates between the ancilla and the flag qubit. Usually same as `err`.
     */
    static void S2_flag(double *cArray, const double err, const double errTeleport, const double ancillaFlagErr);

    /**
     * @brief The first set of stabiliser operations corresponding to `XIXZZ`. Fault-tolerant.
     *
     * @param cArray The real density matrix in the form of a C-style array of length `pow(2,n-1)+pow(2,2n-1)` where `n=7` is the number of qubits (including the ancilla and flag qubit).
     * @param err The error rate of the two-qubit gates in the stabiliser operation between the ancilla and the data qubits, i.e., the CNOT and CPHASE gates.
     * @param errTeleport The error rate of the teleported two-qubit gates in 3 of the 4 two-qubit gates. This is typically higher than `err`.
     * @param ancillaFlagErr The error rate of the two-qubit gates between the ancilla and the flag qubit. Usually same as `err`.
     */
    static void S3_flag(double *cArray, const double err, const double errTeleport, const double ancillaFlagErr);

    /**
     * @brief The first set of stabiliser operations corresponding to `ZXIXZ`. Fault-tolerant.
     *
     * @param cArray The real density matrix in the form of a C-style array of length `pow(2,n-1)+pow(2,2n-1)` where `n=7` is the number of qubits (including the ancilla and flag qubit).
     * @param err The error rate of the two-qubit gates in the stabiliser operation between the ancilla and the data qubits, i.e., the CNOT and CPHASE gates.
     * @param errTeleport The error rate of the teleported two-qubit gates in 3 of the 4 two-qubit gates. This is typically higher than `err`.
     * @param ancillaFlagErr The error rate of the two-qubit gates between the ancilla and the flag qubit. Usually same as `err`.
     */
    static void S4_flag(double *cArray, const double err, const double errTeleport, const double ancillaFlagErr);

    /**
     * @brief The first set of stabiliser operations corresponding to `XZZXI`. Not fault-tolerant.
     *
     * @param cArray The real density matrix in the form of a C-style array of length `pow(2,n-1)+pow(2,2n-1)` where `n=6` is the number of qubits.
     * @param err The error rate of the two-qubit gates in the stabiliser operation, i.e., the CNOT and CPHASE gates.
     * @param errTeleport The error rate of the teleported two-qubit gates in 3 of the 4 two-qubit gates. This is typically higher than `err`.
     */
    static void S1(double *cArray, const double err, const double errTeleport);

    /**
     * @brief The second set of stabiliser operations corresponding to `IXZZX`. Not fault-tolerant.
     *
     * @param cArray The real and symmetric density matrix in the form of a C-style array of length `pow(2,n-1)+pow(2,2n-1)` where `n=7` is the number of qubits.
     * @param err The error rate of the two-qubit gates in the stabiliser operation, i.e., the CNOT and CPHASE gates.
     * @param errTeleport The error rate of the teleported two-qubit gates in 3 of the 4 two-qubit gates. This is typically higher than `err`.
     */
    static void S2(double *cArray, const double err, const double errTeleport);

    /**
     * @brief The third set of stabiliser operations corresponding to `XIXZZ`. Not fault-tolerant.
     *
     * @param cArray The real and symmetric density matrix in the form of a C-style array of length `pow(2,n-1)+pow(2,2n-1)` where `n=8` is the number of qubits.
     * @param err The error rate of the two-qubit gates in the stabiliser operation, i.e., the CNOT and CPHASE gates.
     * @param errTeleport The error rate of the teleported two-qubit gates in 3 of the 4 two-qubit gates. This is typically higher than `err`.
     */
    static void S3(double *cArray, const double err, const double errTeleport);

    /**
     * @brief The fourth set of stabiliser operations corresponding to `ZXIXZ`. Not fault-tolerant.
     *
     * @param cArray The real and symmetric density matrix in the form of a C-style array of length `pow(2,n-1)+pow(2,2n-1)` where `n=9` is the number of qubits.
     * @param err The error rate of the two-qubit gates in the stabiliser operation, i.e., the CNOT and CPHASE gates.
     * @param errTeleport The error rate of the teleported two-qubit gates in 3 of the 4 two-qubit gates. This is typically higher than `err`.
     */
    static void S4(double *cArray, const double err, const double errTeleport);

    /**
     * @brief The first set of stabiliser operations corresponding to `XZZXI`. Not fault-tolerant.
     *
     * @param cArray The real density matrix in the form of a C-style array of length `pow(2,n-1)+pow(2,2n-1)` where `n=6` is the number of qubits.
     * @param err The error rate of the two-qubit gates in the stabiliser operation, i.e., the CNOT and CPHASE gates.
     * @param errTeleport The error rate of the teleported two-qubit gates in 3 of the 4 two-qubit gates. This is typically higher than `err`.
     */
    static void S1Sequential(double *cArray, const double err, const double errTeleport);

    static std::pair<FastPauli::CArr<5>, FastPauli::CArr<5>> S1Sequential(FastPauli::CArr<5> const &cArr5, const double err, const double errTeleport);

    /**
     * @brief The second set of stabiliser operations corresponding to `IXZZX`. Not fault-tolerant.
     *
     * @param cArray The real density matrix in the form of a C-style array of length `pow(2,n-1)+pow(2,2n-1)` where `n=6` is the number of qubits.
     * @param err The error rate of the two-qubit gates in the stabiliser operation, i.e., the CNOT and CPHASE gates.
     * @param errTeleport The error rate of the teleported two-qubit gates in 3 of the 4 two-qubit gates. This is typically higher than `err`.
     */
    static void S2Sequential(double *cArray, const double errTeleport);

    static std::pair<FastPauli::CArr<5>, FastPauli::CArr<5>> S2Sequential(FastPauli::CArr<5> const &cArr5, const double errTeleport);

    /**
     * @brief The second set of stabiliser operations corresponding to `XIXZZ`. Not fault-tolerant.
     *
     * @param cArray The real density matrix in the form of a C-style array of length `pow(2,n-1)+pow(2,2n-1)` where `n=6` is the number of qubits.
     * @param err The error rate of the two-qubit gates in the stabiliser operation, i.e., the CNOT and CPHASE gates.
     * @param errTeleport The error rate of the teleported two-qubit gates in 3 of the 4 two-qubit gates. This is typically higher than `err`.
     */
    static void S3Sequential(double *cArray, const double errTeleport);

    static std::pair<FastPauli::CArr<5>, FastPauli::CArr<5>> S3Sequential(FastPauli::CArr<5> const &cArr5, const double errTeleport);

    /**
     * @brief The second set of stabiliser operations corresponding to `ZXIXZ`. Not fault-tolerant.
     *
     * @param cArray The real density matrix in the form of a C-style array of length `pow(2,n-1)+pow(2,2n-1)` where `n=6` is the number of qubits.
     * @param err The error rate of the two-qubit gates in the stabiliser operation, i.e., the CNOT and CPHASE gates.
     * @param errTeleport The error rate of the teleported two-qubit gates in 3 of the 4 two-qubit gates. This is typically higher than `err`.
     */
    static void S4Sequential(double *cArray, const double errTeleport);

    static std::pair<FastPauli::CArr<5>, FastPauli::CArr<5>> S4Sequential(FastPauli::CArr<5> const &cArr5, const double errTeleport);

    /**
     * @brief Apply one round of non fault tolerant [[5,1,3]] code correction to the input 5-qubit density matrix `cArrCodeword` with error rate of `err`.
     *
     * @param cArr6 The 6-qubit real density matrix in the form of a C-style array input to the stabiliser operation `XZZXI`.
     * @param cArr7 The 7-qubit real density matrix in the form of a C-style array input to the stabiliser operation `IXZZX`.
     * @param cArr8 The 8-qubit real density matrix in the form of a C-style array input to the stabiliser operation `XIXZZ`.
     * @param cArr9 The 9-qubit real density matrix in the form of a C-style array input to the stabiliser operation `ZXIXZ`.
     * @param cArrCodeword The 5-qubit real density matrix in the form of a C-style array which acts as the starting point of the error correction.
     * @param cArrH1 The 1-qubit density matrix in the form of a C-style array which acts as the state of each ancilla qubit, i.e., |+⟩=(|0⟩+|1⟩)/√2.
     * @param err The error rate of the two-qubit gates in the stabiliser operation, i.e., the CNOT and CPHASE gates.
     * @param errTeleport The error rate of the teleported two-qubit gates in 3 of the 4 two-qubit gates. This is typically higher than `err`.
     */
    static void OneRound(double *cArr6, double *cArr7, double *cArr8, double *cArr9, const double *cArrCodeword, const double *cArrH1, const double err, const double errTeleport);

    static FastPauli::CArr<5> OneRoundSequential(const double *cArrCodeword, const double err, const double errTeleport, const FiveQubitCode::correctionMode correctionMode, const int qubitLostIndex);

    static FastPauli::CArr<5> OneRoundSequential_flag(FiveQubitCode::ECInfo &ecInfo, double *logicalKet, double *cArrCodeword, const double err, const double errTeleport, const double ancillaFlagErr, const int ecInfoIndex);

    static FastPauli::CArr<5> OneRoundSequential_flag(double *cArrCodeword, const double err, const double errTeleport, const double ancillaFlagErr);

    /**
     * @brief Returns an array that corresponds to the ket vector |+⟩^(⊗numQubits).
     *
     * @tparam numQubits
     * @tparam Size
     * @return std::array<double, Size>
     */
    template <int numQubits, size_t Size = UINT64_C(1) << numQubits>
    static std::array<double, Size> getKetH() {  // Gets a ket vector.
        const double elemVal = 1.0 / std::pow(std::sqrt(2), numQubits);
        std::array<double, Size> array1D;
        for (auto &elem : array1D) {
            elem = elemVal;
        }
        return array1D;
    }

    /**
     * @brief Gets a C-style array which represents a real density matrix (|+⟩⟨+|)^(⊗numQubits).
     *
     * @tparam numQubits
     * @tparam Size
     * @return std::array<double, Size>
     */
    template <int numQubits, size_t Size = FastPauli::CArrLen<numQubits>>
    static std::array<double, Size> getH() {  // Gets a symmetric density matrix C-style.
        std::array<double, Size> cArray;
        const double elemVal = 1.0 / static_cast<double>(UINT64_C(1) << numQubits);
        for (auto &elem : cArray) {
            elem = elemVal;
        }
        return cArray;
    }

    static std::array<std::array<double, 16>, 16> getCorrectionArrays();

    static FastPauli::CArr<5> applyCorrectionArrays(const double *cArr9, std::array<std::array<double, 16>, 16> const &ancillaKets);

    /**
     * @brief Get the std::array<double, 10> result which represents DensityMatrix(ancilla,flag). For example, if `isAncillaTriggered = true` and `isFlagTriggered = true`, then the result is DensityMatrix(+,0).
     *
     * @tparam isAncillaTriggered true if ancilla triggered, else false.
     * @tparam isFlagTriggered true if flag triggered, else false.
     * @return constexpr std::array<double, FastPauli::CArrLen<2>>
     */
    template <bool isAncillaTriggered, bool isFlagTriggered>
    static constexpr std::array<double, 10> getAncillaFlagDensityMatrix() {
        if (!isAncillaTriggered && !isFlagTriggered) {  // DensityMatrix(+,0)
            return {0.5, 0, 0.5, 0, 0, 0, 0, 0.5, 0, 0};
        } else if (isAncillaTriggered && !isFlagTriggered) {  // DensityMatrix(-,0)
            return {0.5, 0, -0.5, 0, 0, 0, 0, 0.5, 0, 0};
        } else if (!isAncillaTriggered && isFlagTriggered) {  // DensityMatrix(+,1)
            return {0, 0, 0, 0, 0.5, 0, 0.5, 0, 0, 0.5};
        } else {  // DensityMatrix(-,1)
            return {0, 0, 0, 0, 0.5, 0, -0.5, 0, 0, 0.5};
        }
    }

    template <bool isAncillaTriggered, bool isFlagTriggered>
    static constexpr std::array<double, 4> getAncillaFlagKet() {
        if (!isAncillaTriggered && !isFlagTriggered) {  // Ket(+,0)
            return {std::numbers::sqrt2 / 2, 0, std::numbers::sqrt2 / 2, 0};
        } else if (isAncillaTriggered && !isFlagTriggered) {  // Ket(-,0)
            return {std::numbers::sqrt2 / 2, 0, -std::numbers::sqrt2 / 2, 0};
        } else if (!isAncillaTriggered && isFlagTriggered) {  // Ket(+,1)
            return {0, std::numbers::sqrt2 / 2, 0, std::numbers::sqrt2 / 2};
        } else {  // Ket(-,1)
            return {0, std::numbers::sqrt2 / 2, 0, -std::numbers::sqrt2 / 2};
        }
    }

    /**
     * @brief Get the Logical Ket vector for the [[5,1,3]] code given by |ψ_L⟩ = cos(θ/2)|0_L⟩ + sin(θ/2)|1_L⟩ which has length pow(2,5) = 32.
     *
     * @param theta A real number θ such that |ψ_L⟩ = cos(θ/2)|0_L⟩ + sin(θ/2)|1_L⟩.
     * @return std::array<double, 32>
     */
    static std::array<double, 32> getLogicalKet(const double theta);

    static FastPauli::CArr<5> getLogicalDensityMatrix(const double theta);

    static FastPauli::CArr<5> OneRound_flag(ECInfo &ecInfo, double *logicalKet, double *cArrCodeword, const double err, const double errTeleport, const double ancillaFlagErr, std::array<std::array<double, 16>, 16> const &ancillaKets, const int ecInfoIndex);
    // static void OneRound_flag(double *logicalKet, double *cArrCodeword, const double err, const double errTeleport, const double ancillaFlagErr, std::array<std::array<double, 16>, 16> const &ancillaKets);

    static FastPauli::CArr<5> OneRound_1erasure(double *logicalKet, const double *cArrCodeword, const double err, const double errTeleport, const int k);

    static FastPauli::CArr<5> applyCorrectionArrays_flag(const double *cArr9, std::array<std::array<double, 16>, 16> const &ancillaKets, const uint16_t flag);

    static FastPauli::CArr<5> applyCorrectionArrays_1erasure(const double *cArr9, const int k);

    /**
     * @brief This function basically does `(⟨ancilla|⊗𝕀...⊗𝕀).ρ.(|ancilla⟩⊗𝕀...⊗𝕀)` where ρ is the input real density matrix and 𝕀 is a 2×2 identity matrix.
     *
     * @param cArr5 The output density matrix with reduced dimension of length `pow(2,numQubits-1)+pow(2,2numQubits-1)`.
     * @param numAncillaQubits The number of ancilla qubits.
     * @param numQubits The number of qubits in the resulting reduced density matrix.
     * @param ancillaKet The ket vector of the ancilla qubits, i.e., |ancilla⟩.
     * @param cArr9 The input density matrix with the full dimension of length `pow(2,(numQubits+numAncillaQubits)-1)+pow(2,2(numQubits+numAncillaQubits)-1)`.
     */
    static void applyAncillaKet(double *cArr5, const int numAncillaQubits, const int numQubits, const double *ancillaKet, const double *cArr9);

    template <int numAncillaQubits, int numQubits, size_t Size = FastPauli::CArrLen<numQubits>>
    static std::array<double, Size> applyAncillaKet(const double *ancillaKet, const double *cArr9) {
        std::array<double, Size> cArr5;
        FiveQubitCode::applyAncillaKet(cArr5.data(), numAncillaQubits, numQubits, ancillaKet, cArr9);
        return cArr5;
    }

    static FastPauli::CArr<5> getDensityMatrixW0(const double *cArray);
    static FastPauli::CArr<5> getDensityMatrixW1(const double *cArray);
    static FastPauli::CArr<5> getDensityMatrixW2(const double *cArray);
};

#endif