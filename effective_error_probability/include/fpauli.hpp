#ifndef _FPAULIHPP_
#define _FPAULIHPP_

#include <algorithm>  // Provides std::copy
#include <array>
#include <cassert>
#include <cmath>  // Provides std::max
#include <cstddef>
#include <cstdint>  // Provides definition for uint64_t.
#include <iomanip>
#include <iostream>
#include <vector>
// #include <iterator> // std::begin() and std::end()

class FastPauli {
   private:
    static void apply_permutation(double *v, int *indices, const int vSize);
    // static int getCol(int i, int j);
    static int fastMod(const int i, const int t);

   public:
    template <int numQubits>
    static constexpr uint64_t CArrLen = (UINT64_C(1) << (numQubits - 1)) + (UINT64_C(1) << (2 * numQubits - 1));

    template <int numQubits>
    using CArr = std::array<double, CArrLen<numQubits>>;

    static size_t getCycleLeftInteger(const uint64_t i, const int rotateAmount, const int numQubits);
    static size_t getCycleRightInteger(const uint64_t i, const int rotateAmount, const int numQubits);
    static std::vector<double> getCycledDM(const double *cArray, const int rotateAmount, const int numQubits);

    template <int numQubits, size_t Size = CArrLen<numQubits>, uint64_t n = UINT64_C(1) << numQubits>
    static std::array<double, Size> getCycledDM(const double *cArray, const int rotateAmount) {
        // Initialise variables.
        int cIndex, cIndex2, tmpI, tmpJ;

        std::array<double, Size> result;
        std::copy(cArray, cArray + Size, result.data());

        // Determine whether to rotate left (negative) or right (positive).
        if (rotateAmount < 0) {
            for (size_t i = 0; i < n; i++) {
                for (size_t j = i; j < n; j++) {
                    cIndex = FastPauli::getRow(i, j, n);
                    tmpI = FastPauli::getCycleLeftInteger(i, -rotateAmount, numQubits);
                    tmpJ = FastPauli::getCycleLeftInteger(j, -rotateAmount, numQubits);
                    cIndex2 = tmpI < tmpJ ? FastPauli::getRow(tmpI, tmpJ, n) : FastPauli::getRow(tmpJ, tmpI, n);
                    result[cIndex] = cArray[cIndex2];
                }
            }
        } else if (rotateAmount > 0) {
            for (size_t i = 0; i < n; i++) {
                for (size_t j = i; j < n; j++) {
                    cIndex = FastPauli::getRow(i, j, n);
                    tmpI = FastPauli::getCycleRightInteger(i, rotateAmount, numQubits);
                    tmpJ = FastPauli::getCycleRightInteger(j, rotateAmount, numQubits);
                    cIndex2 = tmpI < tmpJ ? FastPauli::getRow(tmpI, tmpJ, n) : FastPauli::getRow(tmpJ, tmpI, n);
                    result[cIndex] = cArray[cIndex2];
                }
            }
        }

        return result;
    }

    /**
     * @brief Numbers congruent to mod 2^t. In Matlab, it would be '2*i - 2^t + mod(-i,2^t)'.
     * Note that this uses index starting from 1 as opposed to the C/C++ function here which uses index starting from 0.
     * congMod(i,t) gives the numbers that are congruent to {0,1,...,2^t-1} mod 2^(t+1).
     *
     * In Mathematica, you can do:
     * CongMod[i_,exponent_]:=2*i - 2^exponent + Mod[-i,2^exponent];
     * nq=4;
     * Do[Print[Table[CongMod[i,nq-nt],{i,1,2^(nq-1)}]];,{nt,1,nq}]
     *
     * @param i Integer.
     * @param t Exponent.
     * @return int
     */
    static int congMod(const int i, const int t);
    static int getRow(int i, int j, const int n);
    // static int log2_64(uint64_t value);
    // static constexpr uint64_t getCArrLen();

    static const uint64_t getCArrLen(const int numQubits);

    /**
     * @brief This function basically does `DensityMatrixNew = DensityMatrixOld1⊗DensityMatrixOld2` where `DensityMatrixXXX` is a density matrix corresponding to `cArrayXXX`.
     *
     * @param cArrayNew The results will be stored in this. (Modified in-place)
     * @param cArrayOld1 Density matrix (in cArray form) to be processed.
     * @param cArrayOld2 Density matrix (in cArray form) to be processed.
     * @param numQubits1 Number of qubits in cArrayOld1.
     * @param numQubits2 Number of qubits in cArrayOld2.
     */
    static void kron(double *cArrayNew, const double *cArrayOld1, const double *cArrayOld2, const int numQubits1, const int numQubits2);

    template <int numQubits1, int numQubits2, size_t Size = CArrLen<numQubits1 + numQubits2>>
    static std::array<double, Size> kron(const double *cArrayOld1, const double *cArrayOld2) {
        std::array<double, Size> result;
        kron(result.data(), cArrayOld1, cArrayOld2, numQubits1, numQubits2);
        return result;
    }

    static void ketToCArr(double *cArray, const double *ket, const int numQubits);
    static void ketKron(double *array1DNew, const double *array1DOld1, const double *array1DOld2, const int numQubits1, const int numQubits2);
    static void ketX_r(double *array1D, const int numQubits, const int targetQubit);
    static void ketZ_r(double *array1D, const int numQubits, const int targetQubit);
    static void ketY_r(double *array1D, const int numQubits, const int targetQubit);  // Up to some global phase!
    static void X_r(double *cArray, const int numQubits, const int targetQubit);
    static void Z_r(double *cArray, const int numQubits, const int targetQubit);
    static void Y_r(double *cArray, const int numQubits, const int targetQubit);
    static void H_r(double *cArray, const int numQubits, const int targetQubit);
    static std::vector<double> Reset_r(const double *cArray, const int numQubits, const int targetQubit);
    template <int numQubits, size_t Size = CArrLen<numQubits>>
    static std::array<double, Size> Reset_r(const double *cArray, const int targetQubit) {
        std::array<double, Size> outCArray = {};  // Initialise to zero.
        const int n = 1 << (numQubits - 1);
        const int n2 = 1 << (numQubits - targetQubit - 1);
        int tmpI, tmpJ, tmpIndex, tmpIndex2;
        for (int i = 0; i < n; i++) {
            for (int j = i; j < n; j++) {
                tmpI = congMod(i, numQubits - targetQubit - 1);
                tmpJ = congMod(j, numQubits - targetQubit - 1);
                tmpIndex = getRow(tmpI, tmpJ, n << 1);
                tmpIndex2 = getRow(tmpI + n2, tmpJ + n2, n << 1);
                outCArray[tmpIndex] = cArray[tmpIndex] + cArray[tmpIndex2];
            }
        }
        return outCArray;
    }

    static void CZ_r(double *cArray, const int numQubits, const int controlQubit, const int targetQubit);
    static void CX_r(double *cArray, const int numQubits, const int controlQubit, const int targetQubit);
    static void toDM(const double *Array1D, double *cArray, const int numQubits);
    static double fidel_r(const double *Array1D, const double *cArray, const int numQubits);
    static double trace_r(const double *cArray, const int numQubits);
    static double trace2_r(const double *cArray, const int numQubits);
    static void printUTri(const double *cArray, const int numQubits);
};

#endif