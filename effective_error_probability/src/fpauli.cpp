#include "fpauli.hpp"

// #define NUMTHREADS 8  // Number of CPU threads

// int FastPauli::log2_64(uint64_t value){ // See: https://stackoverflow.com/questions/11376288/fast-computing-of-log2-for-64-bit-integers
//     const int tab64[64] = {
//         63,  0, 58,  1, 59, 47, 53,  2,
//         60, 39, 48, 27, 54, 33, 42,  3,
//         61, 51, 37, 40, 49, 18, 28, 20,
//         55, 30, 34, 11, 43, 14, 22,  4,
//         62, 57, 46, 52, 38, 26, 32, 41,
//         50, 36, 17, 19, 29, 10, 13, 21,
//         56, 45, 25, 31, 35, 16,  9, 12,
//         44, 24, 15,  8, 23,  7,  6,  5};
//     value |= value >> 1;
//     value |= value >> 2;
//     value |= value >> 4;
//     value |= value >> 8;
//     value |= value >> 16;
//     value |= value >> 32;
//     return tab64[((uint64_t)((value - (value >> 1))*0x07EDD5E59A4E28C2)) >> 58];
// }

// Should probably put in OEIS. Reference: https://oeis.org/A080414
size_t FastPauli::getCycleLeftInteger(const uint64_t i, const int rotateAmount, const int numQubits) {
    const uint64_t selectBits = (UINT64_C(1) << numQubits) - 1;
    const int modulo = (rotateAmount % numQubits);
    return ((i >> modulo) & selectBits) | ((i << (numQubits - modulo)) & selectBits);
}

// Should probably put in OEIS. Reference: https://oeis.org/A080414
size_t FastPauli::getCycleRightInteger(const uint64_t i, const int rotateAmount, const int numQubits) {
    const uint64_t selectBits = (UINT64_C(1) << numQubits) - 1;
    const int modulo = (rotateAmount % numQubits);
    return ((i << modulo) & selectBits) | ((i >> (numQubits - modulo)) & selectBits);
}

std::vector<double> FastPauli::getCycledDM(const double *cArray, const int rotateAmount, const int numQubits) {
    // Initialise variables.
    int cIndex, cIndex2, tmpI, tmpJ;
    const int n = 1 << numQubits;
    const uint64_t cArrLen = FastPauli::getCArrLen(numQubits);

    // Copy dynamic C-style array to std::array. Reference: https://stackoverflow.com/questions/259297/how-do-you-copy-the-contents-of-an-array-to-a-stdvector-in-c-without-looping
    std::vector<double> result;
    result.insert(result.end(), &cArray[0], &cArray[cArrLen]);

    // Determine whether to rotate left (negative) or right (positive).
    if (rotateAmount < 0) {
        for (int i = 0; i < n; i++) {
            for (int j = i; j < n; j++) {
                cIndex = FastPauli::getRow(i, j, n);
                tmpI = FastPauli::getCycleLeftInteger(i, -rotateAmount, numQubits);
                tmpJ = FastPauli::getCycleLeftInteger(j, -rotateAmount, numQubits);
                cIndex2 = tmpI < tmpJ ? FastPauli::getRow(tmpI, tmpJ, n) : FastPauli::getRow(tmpJ, tmpI, n);
                result[cIndex] = cArray[cIndex2];
            }
        }
    } else if (rotateAmount > 0) {
        for (int i = 0; i < n; i++) {
            for (int j = i; j < n; j++) {
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

const uint64_t FastPauli::getCArrLen(const int numQubits) {
    return (UINT64_C(1) << (numQubits - 1)) + (UINT64_C(1) << (2 * numQubits - 1));
}

void FastPauli::kron(double *cArrayNew, const double *cArrayOld1, const double *cArrayOld2, const int numQubits1, const int numQubits2) {  // Original. Derived myself.
    const int numQubitsTot = numQubits1 + numQubits2;
    const int ntot = 1 << numQubitsTot;
    const int n1 = 1 << numQubits1;
    const int n2 = 1 << numQubits2;
    int i, j;
    int tmp, tmp0, tmp1;
    int mi, mj;
    int newInd, oldInd1, oldInd2;
    // #if defined(_OPENMP)
    //     #pragma omp parallel for shared(cArrayNew, cArrayOld1, cArrayOld2) private(i, j, mi, mj, tmp, tmp0, tmp1, newInd, oldInd1, oldInd2) schedule(static, std::max(1, ntot / NUMTHREADS))
    // #endif
    for (i = 0; i < ntot; i++) {
        for (j = i; j < ntot; j++) {
            mi = FastPauli::fastMod(i, numQubits2);
            mj = FastPauli::fastMod(j, numQubits2);
            tmp = FastPauli::fastMod(i, numQubits2) - FastPauli::fastMod(j, numQubits2);
            tmp0 = tmp < 0 ? 0 : tmp;
            tmp = FastPauli::fastMod(i - mj, numQubits2);
            tmp1 = tmp > (n2 - mj - 1) ? 0 : tmp;
            newInd = FastPauli::getRow(i, j, ntot);
            oldInd1 = FastPauli::getRow(i / n2, j / n2, n1);
            oldInd2 = FastPauli::getRow(mi - tmp1, mj + tmp0, n2);
            cArrayNew[newInd] = cArrayOld1[oldInd1] * cArrayOld2[oldInd2];
        }
    }
}

void FastPauli::ketKron(double *array1DNew, const double *array1DOld1, const double *array1DOld2, const int numQubits1, const int numQubits2) {
    const int numQubitsTot = numQubits1 + numQubits2;
    const int ntot = 1 << numQubitsTot;
    // const int n1 = 1 << numQubits1; // Not used.
    const int n2 = 1 << numQubits2;
    for (int i = 0; i < ntot; i++) {
        array1DNew[i] = array1DOld1[i / n2] * array1DOld2[FastPauli::fastMod(i, numQubits2)];
    }
}

void FastPauli::ketToCArr(double *cArray, const double *ket, const int numQubits) {
    const int n = 1 << numQubits;
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            cArray[FastPauli::getRow(i, j, n)] = ket[i] * ket[j];
        }
    }
}

void FastPauli::apply_permutation(double *v, int *indices, const int vSize) {  // Reference: https://stackoverflow.com/questions/67751784/how-to-do-in-place-sorting-a-list-according-to-a-given-index-in-c.
    for (int i = 0; i < vSize; i++) {
        auto current = i;
        while (i != indices[current]) {
            auto next = indices[current];
            std::swap(v[current], v[next]);
            indices[current] = current;
            current = next;
        }
        indices[current] = current;
    }
}

int FastPauli::getRow(int i, int j, const int n) {
    // Note: 1 ≤ i, j ≤ n and i ≤ j.
    // Note: n here refers to the dimension of the square matrix, not the number of qubits.
    i = i + 1;
    j = j + 1;
    assert(1 <= i && j <= n && i <= j);
    return j + n * (i - 1) - ((i * (i - 1)) >> 1) - 1;  // Reference: https://www.geeksforgeeks.org/convert-given-upper-triangular-matrix-to-1d-array/
}

// int getCol(int i, int j){
//     i = i + 1;
//     j = j + 1;
//     // std::cout << "i: " << i << ", j: " << j << "\n";
//     assert(1 <= i && i <= j);
//     // Note: 1 ≤ i, j ≤ n and i ≤ j.
//     // Note: n here refers to the dimension of the square matrix, not the number of qubits.
//     return i + ((j*(j-1)) >> 1) - 1; // Reference: https://www.geeksforgeeks.org/convert-given-upper-triangular-matrix-to-1d-array/
// }

int FastPauli::fastMod(const int i, const int t) {
    return (i) & ((1 << t) - 1);
}

int FastPauli::congMod(const int i, const int t) {
    // Numbers congruent to mod 2^t;
    // In Matlab, it would be '2*i - 2^t + mod(-i,2^t)'. Note that this uses index starting from 1 as opposed to the C/C++ function here which uses index starting from 0.
    // congMod(i,t) gives the numbers that are congruent to {0,1,...,2^t-1} mod 2^(t+1).
    const int pow2t = 1 << t;
    return ((i + 1) << 1) - pow2t + ((-i - 1) & (pow2t - 1)) - 1;  // If b is a power of 2, it is possible to compute n % b as n & (b-1). Reference: https://stackoverflow.com/questions/7594508/modulo-operator-with-negative-values.
}

void FastPauli::ketX_r(double *array1D, const int numQubits, const int targetQubit) {
    int i;
    const int n = 1 << numQubits;  // Matrix dimension.
    const int tmp0 = (n >> (targetQubit + 1));
    int *ketIndex = new int[n];
    for (i = 0; i < n; i++) {
        ketIndex[i] = i ^ tmp0;
    }
    FastPauli::apply_permutation(array1D, ketIndex, n);
    delete[] ketIndex;
}

void FastPauli::ketZ_r(double *array1D, const int numQubits, const int targetQubit) {
    int i, ketIndex;
    const int n = 1 << numQubits;  // Matrix dimension.
    const int n2 = n >> 1;
    const int tmp0 = (n >> (targetQubit + 1));
    for (i = 0; i < n2; i++) {
        ketIndex = n + 1 - (2 * (i + 1) - tmp0 + FastPauli::fastMod(-i - 1, numQubits - targetQubit - 1)) - 1;
        array1D[ketIndex] = -array1D[ketIndex];
    }
}

void FastPauli::ketY_r(double *array1D, const int numQubits, const int targetQubit) {  // Up to some global phase!
    FastPauli::ketZ_r(array1D, numQubits, targetQubit);
    FastPauli::ketX_r(array1D, numQubits, targetQubit);
}

void FastPauli::X_r(double *cArray, const int numQubits, const int targetQubit) {
    const int cArrLen = (UINT64_C(1) << (numQubits - 1)) + (UINT64_C(1) << (2 * numQubits - 1));  // numQubits cannot be too high, otherwise overflow and cArrLen becomes negative!
    int i;
    int j;
    // int counter = 0;
    int col;
    int row;
    int tmp;
    int *cIndex = new int[cArrLen];
    const int n = 1 << numQubits;  // Matrix dimension.
    const int tmp0 = (n >> (targetQubit + 1));
    // #if defined(_OPENMP)
    //     #pragma omp parallel for shared(cIndex) private(i, j, col, row, tmp) schedule(static, std::max(1, n / NUMTHREADS))
    // #endif
    for (i = 0; i < n; i++) {
        col = i ^ tmp0;
        for (j = i + 1; j <= n; j++) {
            row = (j - 1) ^ tmp0;
            tmp = FastPauli::getRow(i, j - 1, n);
            cIndex[tmp] = (row > col) ? FastPauli::getRow(col, row, n) : FastPauli::getRow(row, col, n);
            // counter++;
            // cArray[cIndex] = -cArray[cIndex];
            // std::cout << "\033[1;33m(col,row): (" << col << "," << row << ")\033[0m\n"; // Debug statement.
        }
    }
    FastPauli::apply_permutation(cArray, cIndex, cArrLen);
    delete[] cIndex;
}

void FastPauli::Z_r(double *cArray, const int numQubits, const int targetQubit) {
    int i;
    int j;
    int col;
    int row;
    int cIndex;
    const int n = 1 << numQubits;  // Matrix dimension.
    const int n2 = n >> 1;
    const int tmp0 = (n >> (targetQubit + 1));
    // #if defined(_OPENMP)
    //     #pragma omp parallel for shared(cArray) private(i, j, col, row, cIndex) schedule(static, std::max(1, n2 / NUMTHREADS))
    // #endif
    for (i = 0; i < n2; i++) {
        col = 2 * (i + 1) - tmp0 + FastPauli::fastMod(-i - 1, numQubits - targetQubit - 1) - 1;
        // std::cout << "col: " << col << "\n"; // Debug statement.
        for (j = 0; j < n2; j++) {
            row = n + 1 - (2 * (j + 1) - tmp0 + FastPauli::fastMod(-j - 1, numQubits - targetQubit - 1)) - 1;
            cIndex = (row > col) ? FastPauli::getRow(col, row, n) : FastPauli::getRow(row, col, n);
            cArray[cIndex] = -cArray[cIndex];
            // std::cout << "\033[1;33m(col,row): (" << col << "," << row << ")\033[0m\n"; // Debug statement.
        }
    }
}

void FastPauli::Y_r(double *cArray, const int numQubits, const int targetQubit) {
    FastPauli::Z_r(cArray, numQubits, targetQubit);
    FastPauli::X_r(cArray, numQubits, targetQubit);
}

void FastPauli::H_r(double *cArray, const int numQubits, const int targetQubit) {
    const int cArrLen = (UINT64_C(1) << (numQubits - 1)) + (UINT64_C(1) << (2 * numQubits - 1));  // numQubits cannot be too high, otherwise overflow and cArrLen becomes negative!
    int ii, iii, jj, tmpInd;
    int tmpcounter;
    int col;
    int row;
    int cIndex;
    const int n = 1 << numQubits;  // Matrix dimension.
    const int n2 = n >> 1;
    const int tmp0 = (n >> (targetQubit + 1));
    double *zz = new double[cArrLen];
    double *zx = new double[cArrLen];
    double *xz = new double[cArrLen];

    // #if defined(_OPENMP)
    //     #pragma omp parallel for shared(zz, cArray) private(i) schedule(static, std::max(1, cArrLen / NUMTHREADS))
    // #endif
    for (int i = 0; i < cArrLen; i++) {
        zz[i] = cArray[i];
    }
    // #if defined(_OPENMP)
    //     #pragma omp simd parallel for shared(xz, cArray) private(col, tmpcounter, cIndex) schedule(static)
    // #endif
    for (int i = 0; i < n; i++) {
        col = i ^ tmp0;
        // std::cout << "col: " << col << "\n";
        for (int j = i; j < n; j++) {
            cIndex = (j > col) ? FastPauli::getRow(col, j, n) : FastPauli::getRow(j, col, n);
            tmpcounter = FastPauli::getRow(i, j, n);
            xz[tmpcounter] = cArray[cIndex];
            // counter++;
            // std::cout << "\033[1;33m(col,row): (" << col << "," << j << ")\033[0m\n"; // Debug statement.
        }
    }
    // #if defined(_OPENMP)
    //     #pragma omp parallel for shared(xz) private(i, ii, j, iii, jj, tmpInd) schedule(static, std::max(1, n2 / NUMTHREADS))
    // #endif
    for (int i = 1; i <= n2; i++) {
        ii = FastPauli::congMod(i - 1, numQubits - targetQubit - 1) + 1;
        for (int j = ii; j <= n; j++) {
            iii = n + 1 - ii;
            jj = n + 1 - j;
            tmpInd = FastPauli::getRow(jj - 1, iii - 1, n);
            xz[tmpInd] = -xz[tmpInd];
            // std::cout << "\033[1;33m(jj,iii): (" << jj-1 << "," << iii-1 << ")\033[0m\n"; // Debug statement.
        }
    }

    // counter = 0;
    // #if defined(_OPENMP)
    //     #pragma omp parallel for shared(zx, cArray) private(i, j, row, tmpcounter, cIndex) schedule(static, std::max(1, n / NUMTHREADS))
    // #endif
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            row = j ^ tmp0;
            cIndex = (row > i) ? FastPauli::getRow(i, row, n) : FastPauli::getRow(row, i, n);
            tmpcounter = FastPauli::getRow(i, j, n);
            zx[tmpcounter] = cArray[cIndex];
            // counter++;
            // std::cout << "\033[1;33m(col,row): (" << i << "," << row << ")\033[0m\n"; // Debug statement.
        }
    }
    // #if defined(_OPENMP)
    //     #pragma omp parallel for shared(zx) private(i, ii, j, tmpInd) schedule(static, std::max(1, n2 / NUMTHREADS))
    // #endif
    for (int i = 1; i <= n2; i++) {
        ii = n - FastPauli::congMod(i - 1, numQubits - targetQubit - 1);
        for (int j = ii; j <= n; j++) {
            tmpInd = FastPauli::getRow(ii - 1, j - 1, n);
            zx[tmpInd] = -zx[tmpInd];
            // std::cout << "\033[1;33m(ii,j): (" << ii-1 << "," << j-1 << ")\033[0m\n"; // Debug statement.
        }
    }

    FastPauli::X_r(cArray, numQubits, targetQubit);  // Use original cArray as xx to reduce std::copy usage.
    FastPauli::Z_r(zz, numQubits, targetQubit);
    // #if defined(_OPENMP)
    //     #pragma omp parallel for shared(cArray, xz, zx, zz) private(i) schedule(static, std::max(1, cArrLen / NUMTHREADS))
    // #endif
    for (int i = 0; i < cArrLen; i++) {
        cArray[i] = (cArray[i] + xz[i] + zx[i] + zz[i]) / 2.0;
    }
    delete[] zz;
    delete[] zx;
    delete[] xz;
}

std::vector<double> FastPauli::Reset_r(const double *cArray, const int numQubits, const int targetQubit) {
    std::vector<double> outCArray(FastPauli::getCArrLen(numQubits), 0);  // Initialise to zero.
    const int n = 1 << (numQubits - 1);
    const int n2 = 1 << (numQubits - targetQubit - 1);
    int tmpI, tmpJ, tmpIndex, tmpIndex2;
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            tmpI = FastPauli::congMod(i, numQubits - targetQubit - 1);
            tmpJ = FastPauli::congMod(j, numQubits - targetQubit - 1);
            tmpIndex = FastPauli::getRow(tmpI, tmpJ, n << 1);
            tmpIndex2 = FastPauli::getRow(tmpI + n2, tmpJ + n2, n << 1);
            outCArray[tmpIndex] = cArray[tmpIndex] + cArray[tmpIndex2];
        }
    }
    return outCArray;
}

void FastPauli::CZ_r(double *cArray, const int numQubits, const int controlQubit, const int targetQubit) {
    assert(controlQubit < targetQubit);
    int i;
    int j;
    int jj;
    int col;
    int row;
    int tmp3;
    int tmp4;
    int tmp5;
    int tmp6;
    int tmp7;
    int tmp8;
    int cIndex;
    const int distance = targetQubit - controlQubit;
    const int totalWidth = 1 << (numQubits - 2);
    const int totalHeight = 3 * totalWidth;
    const int n = 1 << numQubits;  // Matrix dimension.
    const int n1 = n >> 1;
    const int tmp0 = n1 >> distance;  // n1/(2^distance).
    const int tmp1 = n >> distance;   // n/(2^distance).

    const int h1 = n1 + tmp0;  // Height of each divided pillar.
    const int h1p = n1 + tmp1;
    const int h2 = totalHeight - h1;  // Height of the checker part (excluding white spaces).

    // std::cout << "h1: " << h1 << ", h1p: " << h1p << ", h2: " << h2 << "\n"; // Debug.

    tmp3 = 1 + controlQubit + targetQubit - numQubits;
    tmp8 = (tmp3 >= 0) ? tmp3 : -tmp3;
    tmp4 = 1 << tmp8;
    tmp5 = (1 << controlQubit) + (1 << targetQubit);
    tmp7 = n >> (targetQubit + 1);
    // std::cout << "tmp3: " << tmp3 << ", tmp8: " << tmp8 << ", tmp4: " << tmp4 << ", tmp5: " << tmp5 << ", tmp7: " << tmp7 << "\n";
    // OPENMP https://stackoverflow.com/questions/20413995/reducing-on-array-in-openmp
    // #if defined(_OPENMP)
    //     #pragma omp parallel for shared(cArray) private(i, j, col, tmp6, row, cIndex) schedule(static, std::max(1, totalWidth / NUMTHREADS))
    // #endif
    for (i = 0; i < totalWidth; i++) {  // Deal with pillar section of the matrix.
        col = n - FastPauli::congMod(FastPauli::congMod(i, numQubits - targetQubit - 1), numQubits - controlQubit - 1) - 1;
        // std::cout << "col: " << col << "\n"; // Debug statement.
        for (j = 0; j < h1; j++) {
            tmp6 = (tmp3 >= 0) ? (tmp4 * (j)) / tmp5 : (j) / (tmp4 * tmp5);
            row = j + tmp7 * (-1 + (1 << distance)) * tmp6;
            // if (i == totalWidth-1){
            //     std::cout << "row: " << row << "\n"; // Debug statement.
            // }
            cIndex = (row > col) ? FastPauli::getRow(col, row, n) : FastPauli::getRow(row, col, n);
            cArray[cIndex] = -cArray[cIndex];
        }
    }
    // #if defined(_OPENMP)
    //     #pragma omp parallel for shared(cArray) private(i, j, col, row, cIndex, jj) schedule(static, std::max(1, totalWidth / NUMTHREADS))
    // #endif
    for (i = 0; i < totalWidth; i++) {  // Deal with checker section of the matrix.
        col = n - FastPauli::congMod(FastPauli::congMod(i, numQubits - targetQubit - 1), numQubits - controlQubit - 1) - 1;
        for (j = h1 + 1; j <= totalHeight; j++) {
            jj = j - h1;
            row = -(tmp7) + 2 * (jj) + FastPauli::fastMod(-(jj), numQubits - targetQubit - 1) + (1 + ((1 << (controlQubit)) * (jj - 1)) / h2) * h1p / (1 << (controlQubit)) - 1;
            cIndex = (row > col) ? FastPauli::getRow(col, row, n) : FastPauli::getRow(row, col, n);
            cArray[cIndex] = -cArray[cIndex];
            // std::cout << "\033[1;33m(col,row): (" << col << "," << row << ")\033[0m\n"; // Debug statement.
        }
    }
}

void FastPauli::CX_r(double *cArray, const int numQubits, const int controlQubit, const int targetQubit) {
    FastPauli::H_r(cArray, numQubits, targetQubit);
    if (targetQubit > controlQubit) {
        FastPauli::CZ_r(cArray, numQubits, controlQubit, targetQubit);
    } else {
        FastPauli::CZ_r(cArray, numQubits, targetQubit, controlQubit);
    }
    FastPauli::H_r(cArray, numQubits, targetQubit);
}

void FastPauli::toDM(const double *Array1D, double *cArray, const int numQubits) {
    const int n = 1 << numQubits;
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            cArray[FastPauli::getRow(i, j, n)] = Array1D[i] * Array1D[j];
        }
    }
}

double FastPauli::trace_r(const double *cArray, const int numQubits) {
    double out = 0;
    const int n = 1 << numQubits;
    for (int i = 0; i < n; i++) {
        out += cArray[FastPauli::getRow(i, i, n)];
    }
    return out;
}

double FastPauli::trace2_r(const double *cArray, const int numQubits) {
    double out = 0;
    const int n = 1 << numQubits;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            out += 2 * cArray[FastPauli::getRow(i, j, n)] * cArray[FastPauli::getRow(i, j, n)];
        }
    }
    for (int i = 0; i < n; i++) {
        out += cArray[FastPauli::getRow(i, i, n)] * cArray[FastPauli::getRow(i, i, n)];
    }
    return out;
}

double FastPauli::fidel_r(const double *Array1D, const double *cArray, const int numQubits) {
    double out = 0;
    const int n = 1 << numQubits;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            out += 2 * Array1D[j] * Array1D[i] * cArray[FastPauli::getRow(i, j, n)];
        }
    }
    for (int i = 0; i < n; i++) {
        out += Array1D[i] * Array1D[i] * cArray[FastPauli::getRow(i, i, n)];
    }
    return out;
}

void FastPauli::printUTri(const double *cArray, const int numQubits) {
    const int n = 1 << numQubits;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (j >= i) {
                std::cout << "\033[1m";
                std::cout << std::setw(5) << cArray[FastPauli::getRow(i, j, n)];
                std::cout << "\033[0m";
                std::cout << " ";
            } else {
                std::cout << std::setw(5) << 0;
                std::cout << " ";
            }
        }
        std::cout << "\n";
    }
}