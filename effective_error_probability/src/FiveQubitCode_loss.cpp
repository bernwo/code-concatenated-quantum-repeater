#include "FiveQubitCode.hpp"

void FiveQubitCode::applyCorrection_1erasure(const uint16_t i, double *cArr5, const int k2) {
    constexpr int numQubits = 5;
    const int k = (k2 + 1) % numQubits;  // See Table 4.2 in Wo Kah Jen's MSc thesis.
    switch (i) {                         // a1 a2 a3 a4 : Pi where P ∈ {X,Y,Z} is a single-qubit Pauli operator and i is the index of the qubit to correct with i ∈ {0,1,2,3,4};
        case UINT16_C(1):                //  0  0  0  1 : X𝒦
            FastPauli::X_r(cArr5, numQubits, k2);
            break;
        case UINT16_C(2):  //  0  0  1  0 : XkZ𝒦
            FastPauli::X_r(cArr5, numQubits, k);
            FastPauli::Z_r(cArr5, numQubits, k2);
            break;
        case UINT16_C(3):  //  0  0  1  1 : XkY𝒦
            FastPauli::X_r(cArr5, numQubits, k);
            FastPauli::Y_r(cArr5, numQubits, k2);
            break;
        case UINT16_C(4):  //  0  1  0  0 : ZkX𝒦
            FastPauli::Z_r(cArr5, numQubits, k);
            FastPauli::X_r(cArr5, numQubits, k2);
            break;
        case UINT16_C(5):  //  0  1  0  1 : Zk
            FastPauli::Z_r(cArr5, numQubits, k);
            break;
        case UINT16_C(6):  //  0  1  1  0 : YkY𝒦
            FastPauli::Y_r(cArr5, numQubits, k);
            FastPauli::Y_r(cArr5, numQubits, k2);
            break;
        case UINT16_C(7):  //  0  1  1  1 : YkZ𝒦
            FastPauli::Y_r(cArr5, numQubits, k);
            FastPauli::Z_r(cArr5, numQubits, k2);
            break;
        case UINT16_C(8):  //  1  0  0  0 : Xk
            FastPauli::X_r(cArr5, numQubits, k);
            break;
        case UINT16_C(9):  //  1  0  0  1 : XkX𝒦
            FastPauli::X_r(cArr5, numQubits, k);
            FastPauli::X_r(cArr5, numQubits, k2);
            break;
        case UINT16_C(10):  //  1  0  1  0 : Z𝒦
            FastPauli::Z_r(cArr5, numQubits, k2);
            break;
        case UINT16_C(11):  //  1  0  1  1 : Y𝒦
            FastPauli::Y_r(cArr5, numQubits, k2);
            break;
        case UINT16_C(12):  //  1  1  0  0 : YkX𝒦
            FastPauli::Y_r(cArr5, numQubits, k);
            FastPauli::X_r(cArr5, numQubits, k2);
            break;
        case UINT16_C(13):  //  1  1  0  1 : Yk
            FastPauli::Y_r(cArr5, numQubits, k);
            break;
        case UINT16_C(14):  //  1  1  1  0 : ZkY𝒦
            FastPauli::Z_r(cArr5, numQubits, k);
            FastPauli::Y_r(cArr5, numQubits, k2);
            break;
        case UINT16_C(15):  //  1  1  1  1 : ZkZ𝒦
            FastPauli::Z_r(cArr5, numQubits, k);
            FastPauli::Z_r(cArr5, numQubits, k2);
            break;
        default:
            break;
    }
}

std::array<double, FastPauli::CArrLen<5>> FiveQubitCode::applyCorrectionArrays_1erasure(const double *cArr9, const int k) {
    constexpr int numAncillaQubits = 4;
    constexpr int numQubits = 5;
    constexpr uint64_t cArrLen = FastPauli::CArrLen<numQubits>;
    std::array<double, cArrLen> dm;
    dm.fill(0.0);
    std::array<double, cArrLen> tempDm;
    std::vector<double> tempDm_vec;
    double trace;
    std::array<double, 2> ketH1 = FiveQubitCode::getKetH<1>();
    std::array<double, 2> tempKetH1;
    std::array<double, FastPauli::CArrLen<8>> cArr8;
    std::array<double, FastPauli::CArrLen<7>> cArr7;
    std::array<double, FastPauli::CArrLen<6>> cArr6;
    for (uint16_t a = 0; a < 16; a++) {
        // Loop over each ancilla qubit.
        for (uint16_t j = 0; j < UINT16_C(4); j++) {
            tempKetH1 = ketH1;
            if (((a >> (3 - j)) & 0x01)) {
                FastPauli::ketZ_r(tempKetH1.data(), 1, 0);
            }
            switch (j) {
                case UINT16_C(0):
                    FiveQubitCode::applyAncillaKet(cArr8.data(), 1, 8, tempKetH1.data(), cArr9);
                    break;
                case UINT16_C(1):
                    FiveQubitCode::applyAncillaKet(cArr7.data(), 1, 7, tempKetH1.data(), cArr8.data());
                    break;
                case UINT16_C(2):
                    FiveQubitCode::applyAncillaKet(cArr6.data(), 1, 6, tempKetH1.data(), cArr7.data());
                    break;
                case UINT16_C(3):
                    FiveQubitCode::applyAncillaKet(tempDm.data(), 1, numQubits, tempKetH1.data(), cArr6.data());
                default:
                    break;
            }
        }
        // Return to original permutation.
        tempDm_vec = FastPauli::getCycledDM(tempDm.data(), k, numQubits);

        // Apply Pauli correction to tempDm.
        FiveQubitCode::applyCorrection_1erasure(a, tempDm_vec.data(), k);

        // Store corrected tempDm in dm.
        for (uint64_t kk = 0; kk < cArrLen; kk++) {
            dm[kk] += tempDm_vec[kk];
        }
    }

    trace = FastPauli::trace_r(dm.data(), 5);
    for (size_t i = 0; i < FastPauli::CArrLen<5>; i++) {
        dm[i] = dm[i] / trace;
    }

    // Return the corrected 5-qubit density matrix (C-array form).
    return dm;
}

std::array<double, FastPauli::CArrLen<5>> FiveQubitCode::OneRound_1erasure(double *logicalKet, const double *cArrCodeword, const double err, const double errTeleport, const int k) {
    // Initialise the density matrix for a single ancilla qubit.
    std::array<double, FastPauli::CArrLen<1>> cArrH1 = FiveQubitCode::getH<1>();

    // Initialise 5-qubit density matrix std::array's in which we store the results after the reset gate..
    std::vector<double> erasure_k = FastPauli::Reset_r(cArrCodeword, 5, k);

    // Initialise temporary arrays.
    std::array<double, FastPauli::CArrLen<5>> cArr5;  // Stores result.
    std::array<double, FastPauli::CArrLen<6>> cArr6;
    std::array<double, FastPauli::CArrLen<7>> cArr7;
    std::array<double, FastPauli::CArrLen<8>> cArr8;
    std::array<double, FastPauli::CArrLen<9>> cArr9;  // Need to increase stack size. Use /F option in MSVC compiler to specify large stack size.

    // Perform cyclic permutation of the data qubits depending on which qubit was lost.
    erasure_k = FastPauli::getCycledDM(erasure_k.data(), -k, 5);

    // Perform OneRound on each of the erasure_k std::array's.
    FiveQubitCode::OneRound(cArr6.data(), cArr7.data(), cArr8.data(), cArr9.data(), erasure_k.data(), cArrH1.data(), err, errTeleport);

    cArr5 = FiveQubitCode::applyCorrectionArrays_1erasure(cArr9.data(), k);

    return cArr5;
}
