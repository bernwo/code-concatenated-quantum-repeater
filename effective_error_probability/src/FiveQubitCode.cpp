#include "FiveQubitCode.hpp"

void FiveQubitCode::S1(double *cArray, const double err, const double errTeleport) {
    constexpr int numQubits = 6;
    NoisyPauli::NCX_r(cArray, numQubits, 0, 1, err);
    NoisyPauli::NCZ_r(cArray, numQubits, 0, 2, errTeleport);
    NoisyPauli::NCZ_r(cArray, numQubits, 0, 3, errTeleport);
    NoisyPauli::NCX_r(cArray, numQubits, 0, 4, errTeleport);
}

void FiveQubitCode::S2(double *cArray, const double err, const double errTeleport) {
    constexpr int numQubits = 7;
    NoisyPauli::NCX_r(cArray, numQubits, 0, 3, err);
    NoisyPauli::NCZ_r(cArray, numQubits, 0, 4, errTeleport);
    NoisyPauli::NCZ_r(cArray, numQubits, 0, 5, errTeleport);
    NoisyPauli::NCX_r(cArray, numQubits, 0, 6, errTeleport);
}

void FiveQubitCode::S3(double *cArray, const double err, const double errTeleport) {
    constexpr int numQubits = 8;
    NoisyPauli::NCX_r(cArray, numQubits, 0, 3, err);
    NoisyPauli::NCX_r(cArray, numQubits, 0, 5, errTeleport);
    NoisyPauli::NCZ_r(cArray, numQubits, 0, 6, errTeleport);
    NoisyPauli::NCZ_r(cArray, numQubits, 0, 7, errTeleport);
}

void FiveQubitCode::S4(double *cArray, const double err, const double errTeleport) {
    constexpr int numQubits = 9;
    NoisyPauli::NCZ_r(cArray, numQubits, 0, 4, err);
    NoisyPauli::NCX_r(cArray, numQubits, 0, 5, errTeleport);
    NoisyPauli::NCX_r(cArray, numQubits, 0, 7, errTeleport);
    NoisyPauli::NCZ_r(cArray, numQubits, 0, 8, errTeleport);
}

void FiveQubitCode::S1Sequential(double *cArray, const double err, const double errTeleport) {
    constexpr int numQubits = 6;
    NoisyPauli::NCX_r(cArray, numQubits, 0, 1, err);
    NoisyPauli::NCZ_r(cArray, numQubits, 0, 2, errTeleport);
    NoisyPauli::NCZ_r(cArray, numQubits, 0, 3, errTeleport);
    NoisyPauli::NCX_r(cArray, numQubits, 0, 4, errTeleport);
}

std::pair<FastPauli::CArr<5>, FastPauli::CArr<5>> FiveQubitCode::S1Sequential(FastPauli::CArr<5> const &cArr5, const double err, const double errTeleport) {
    // Initialise variables.
    std::pair<FastPauli::CArr<5>, FastPauli::CArr<5>> result;
    std::array<double, FastPauli::CArrLen<6>> tmp;

    // Get ready |+⟩⟨+|.
    std::array<double, FastPauli::CArrLen<1>> cArrH1P = FiveQubitCode::getH<1>();

    // Get ready |+⟩.
    std::array<double, 2> ketH1P = FiveQubitCode::getKetH<1>();

    // Get ready |-⟩.
    std::array<double, 2> ketH1N = ketH1P;
    FastPauli::ketZ_r(ketH1N.data(), 1, 0);

    // Perform stabiliser operations.
    tmp = FastPauli::kron<1, 5>(cArrH1P.data(), cArr5.data());
    FiveQubitCode::S1Sequential(tmp.data(), err, errTeleport);

    // Measure ancilla qubit.
    result.first = FiveQubitCode::applyAncillaKet<1, 5>(ketH1P.data(), tmp.data());
    result.second = FiveQubitCode::applyAncillaKet<1, 5>(ketH1N.data(), tmp.data());
    return result;
}

void FiveQubitCode::S2Sequential(double *cArray, const double errTeleport) {
    constexpr int numQubits = 6;
    NoisyPauli::NCX_r(cArray, numQubits, 0, 2, errTeleport);
    NoisyPauli::NCZ_r(cArray, numQubits, 0, 3, errTeleport);
    NoisyPauli::NCZ_r(cArray, numQubits, 0, 4, errTeleport);
    NoisyPauli::NCX_r(cArray, numQubits, 0, 5, errTeleport);
}

std::pair<FastPauli::CArr<5>, FastPauli::CArr<5>> FiveQubitCode::S2Sequential(FastPauli::CArr<5> const &cArr5, const double errTeleport) {
    // Initialise variables.
    std::pair<FastPauli::CArr<5>, FastPauli::CArr<5>> result;
    std::array<double, FastPauli::CArrLen<6>> tmp;

    // Get ready |+⟩⟨+|.
    std::array<double, FastPauli::CArrLen<1>> cArrH1P = FiveQubitCode::getH<1>();

    // Get ready |+⟩.
    std::array<double, 2> ketH1P = FiveQubitCode::getKetH<1>();

    // Get ready |-⟩.
    std::array<double, 2> ketH1N = ketH1P;
    FastPauli::ketZ_r(ketH1N.data(), 1, 0);

    // Perform stabiliser operations.
    tmp = FastPauli::kron<1, 5>(cArrH1P.data(), cArr5.data());
    FiveQubitCode::S2Sequential(tmp.data(), errTeleport);

    // Measure ancilla qubit.
    result.first = FiveQubitCode::applyAncillaKet<1, 5>(ketH1P.data(), tmp.data());
    result.second = FiveQubitCode::applyAncillaKet<1, 5>(ketH1N.data(), tmp.data());
    return result;
}

void FiveQubitCode::S3Sequential(double *cArray, const double errTeleport) {
    constexpr int numQubits = 6;
    NoisyPauli::NCX_r(cArray, numQubits, 0, 1, errTeleport);
    NoisyPauli::NCX_r(cArray, numQubits, 0, 3, errTeleport);
    NoisyPauli::NCZ_r(cArray, numQubits, 0, 4, errTeleport);
    NoisyPauli::NCZ_r(cArray, numQubits, 0, 5, errTeleport);
}

std::pair<FastPauli::CArr<5>, FastPauli::CArr<5>> FiveQubitCode::S3Sequential(FastPauli::CArr<5> const &cArr5, const double errTeleport) {
    // Initialise variables.
    std::pair<FastPauli::CArr<5>, FastPauli::CArr<5>> result;
    std::array<double, FastPauli::CArrLen<6>> tmp;

    // Get ready |+⟩⟨+|.
    std::array<double, FastPauli::CArrLen<1>> cArrH1P = FiveQubitCode::getH<1>();

    // Get ready |+⟩.
    std::array<double, 2> ketH1P = FiveQubitCode::getKetH<1>();

    // Get ready |-⟩.
    std::array<double, 2> ketH1N = ketH1P;
    FastPauli::ketZ_r(ketH1N.data(), 1, 0);

    // Perform stabiliser operations.
    tmp = FastPauli::kron<1, 5>(cArrH1P.data(), cArr5.data());
    FiveQubitCode::S3Sequential(tmp.data(), errTeleport);

    // Measure ancilla qubit.
    result.first = FiveQubitCode::applyAncillaKet<1, 5>(ketH1P.data(), tmp.data());
    result.second = FiveQubitCode::applyAncillaKet<1, 5>(ketH1N.data(), tmp.data());
    return result;
}

void FiveQubitCode::S4Sequential(double *cArray, const double errTeleport) {
    constexpr int numQubits = 6;
    NoisyPauli::NCZ_r(cArray, numQubits, 0, 1, errTeleport);
    NoisyPauli::NCX_r(cArray, numQubits, 0, 2, errTeleport);
    NoisyPauli::NCX_r(cArray, numQubits, 0, 4, errTeleport);
    NoisyPauli::NCZ_r(cArray, numQubits, 0, 5, errTeleport);
}

std::pair<FastPauli::CArr<5>, FastPauli::CArr<5>> FiveQubitCode::S4Sequential(FastPauli::CArr<5> const &cArr5, const double errTeleport) {
    // Initialise variables.
    std::pair<FastPauli::CArr<5>, FastPauli::CArr<5>> result;
    std::array<double, FastPauli::CArrLen<6>> tmp;

    // Get ready |+⟩⟨+|.
    std::array<double, FastPauli::CArrLen<1>> cArrH1P = FiveQubitCode::getH<1>();

    // Get ready |+⟩.
    std::array<double, 2> ketH1P = FiveQubitCode::getKetH<1>();

    // Get ready |-⟩.
    std::array<double, 2> ketH1N = ketH1P;
    FastPauli::ketZ_r(ketH1N.data(), 1, 0);

    // Perform stabiliser operations.
    tmp = FastPauli::kron<1, 5>(cArrH1P.data(), cArr5.data());
    FiveQubitCode::S4Sequential(tmp.data(), errTeleport);

    // Measure ancilla qubit.
    result.first = FiveQubitCode::applyAncillaKet<1, 5>(ketH1P.data(), tmp.data());
    result.second = FiveQubitCode::applyAncillaKet<1, 5>(ketH1N.data(), tmp.data());
    return result;
}

void FiveQubitCode::OneRound(double *cArr6, double *cArr7, double *cArr8, double *cArr9, const double *cArrCodeword, const double *cArrH1, const double err, const double errTeleport) {
    FastPauli::kron(cArr6, cArrH1, cArrCodeword, 1, 5);
    FiveQubitCode::S1(cArr6, err, errTeleport);
    FastPauli::kron(cArr7, cArrH1, cArr6, 1, 6);
    FiveQubitCode::S2(cArr7, err, errTeleport);
    FastPauli::kron(cArr8, cArrH1, cArr7, 1, 7);
    FiveQubitCode::S3(cArr8, err, errTeleport);
    FastPauli::kron(cArr9, cArrH1, cArr8, 1, 8);
    FiveQubitCode::S4(cArr9, err, errTeleport);
    // double fidel = FiveQubitCode::applyCorrectionArrays(correctionArrays, cArr9);
    // return fidel;
}

FastPauli::CArr<5> FiveQubitCode::OneRoundSequential(const double *cArrCodeword, const double err, const double errTeleport, const FiveQubitCode::correctionMode correctionMode, const int qubitLostIndex) {
    std::array<FastPauli::CArr<5>, 16> cArrAggregate;
    FastPauli::CArr<5> cArrCodewordtmp;
    std::copy(cArrCodeword, cArrCodeword + FastPauli::CArrLen<5>, cArrCodewordtmp.data());
    if (correctionMode == FiveQubitCode::correctionMode::oneLoss) {
        cArrCodewordtmp = FastPauli::getCycledDM<5>(cArrCodewordtmp.data(), -qubitLostIndex);
    }

    // Get ready |+⟩⟨+|.
    std::array<double, 3> cArrH1 = FiveQubitCode::getH<1>();

    // Get ready |+⟩.
    std::array<double, 2> ketH1P = FiveQubitCode::getKetH<1>();

    // Get ready |-⟩.
    std::array<double, 2> ketH1N = ketH1P;
    FastPauli::ketZ_r(ketH1N.data(), 1, 0);

    // Initialise arrays to store measurement results.
    FastPauli::CArr<5> cArr_0;
    FastPauli::CArr<5> cArr_1;
    FastPauli::CArr<5> cArr_00;
    FastPauli::CArr<5> cArr_10;
    FastPauli::CArr<5> cArr_01;
    FastPauli::CArr<5> cArr_11;
    FastPauli::CArr<5> cArr_000;
    FastPauli::CArr<5> cArr_100;
    FastPauli::CArr<5> cArr_010;
    FastPauli::CArr<5> cArr_110;
    FastPauli::CArr<5> cArr_001;
    FastPauli::CArr<5> cArr_101;
    FastPauli::CArr<5> cArr_011;
    FastPauli::CArr<5> cArr_111;
    FastPauli::CArr<5> cArr_0000;
    FastPauli::CArr<5> cArr_1000;
    FastPauli::CArr<5> cArr_0100;
    FastPauli::CArr<5> cArr_1100;
    FastPauli::CArr<5> cArr_0010;
    FastPauli::CArr<5> cArr_1010;
    FastPauli::CArr<5> cArr_0110;
    FastPauli::CArr<5> cArr_1110;
    FastPauli::CArr<5> cArr_0001;
    FastPauli::CArr<5> cArr_1001;
    FastPauli::CArr<5> cArr_0101;
    FastPauli::CArr<5> cArr_1101;
    FastPauli::CArr<5> cArr_0011;
    FastPauli::CArr<5> cArr_1011;
    FastPauli::CArr<5> cArr_0111;
    FastPauli::CArr<5> cArr_1111;
    FastPauli::CArr<6> cArr6;

    // ============================== Stabiliser XZZXI (1) ==============================
    FastPauli::kron(cArr6.data(), cArrH1.data(), cArrCodewordtmp.data(), 1, 5);
    FiveQubitCode::S1Sequential(cArr6.data(), err, errTeleport);
    FiveQubitCode::applyAncillaKet(cArr_0.data(), 1, 5, ketH1P.data(), cArr6.data());
    FiveQubitCode::applyAncillaKet(cArr_1.data(), 1, 5, ketH1N.data(), cArr6.data());

    // ============================== Stabiliser IXZZX (2) ==============================
    FastPauli::kron(cArr6.data(), cArrH1.data(), cArr_0.data(), 1, 5);
    FiveQubitCode::S2Sequential(cArr6.data(), errTeleport);
    FiveQubitCode::applyAncillaKet(cArr_00.data(), 1, 5, ketH1P.data(), cArr6.data());
    FiveQubitCode::applyAncillaKet(cArr_01.data(), 1, 5, ketH1N.data(), cArr6.data());

    FastPauli::kron(cArr6.data(), cArrH1.data(), cArr_1.data(), 1, 5);
    FiveQubitCode::S2Sequential(cArr6.data(), errTeleport);
    FiveQubitCode::applyAncillaKet(cArr_10.data(), 1, 5, ketH1P.data(), cArr6.data());
    FiveQubitCode::applyAncillaKet(cArr_11.data(), 1, 5, ketH1N.data(), cArr6.data());

    // ============================== Stabiliser XIXZZ (3) ==============================
    FastPauli::kron(cArr6.data(), cArrH1.data(), cArr_00.data(), 1, 5);
    FiveQubitCode::S3Sequential(cArr6.data(), errTeleport);
    FiveQubitCode::applyAncillaKet(cArr_000.data(), 1, 5, ketH1P.data(), cArr6.data());
    FiveQubitCode::applyAncillaKet(cArr_001.data(), 1, 5, ketH1N.data(), cArr6.data());

    FastPauli::kron(cArr6.data(), cArrH1.data(), cArr_01.data(), 1, 5);
    FiveQubitCode::S3Sequential(cArr6.data(), errTeleport);
    FiveQubitCode::applyAncillaKet(cArr_010.data(), 1, 5, ketH1P.data(), cArr6.data());
    FiveQubitCode::applyAncillaKet(cArr_011.data(), 1, 5, ketH1N.data(), cArr6.data());

    FastPauli::kron(cArr6.data(), cArrH1.data(), cArr_10.data(), 1, 5);
    FiveQubitCode::S3Sequential(cArr6.data(), errTeleport);
    FiveQubitCode::applyAncillaKet(cArr_100.data(), 1, 5, ketH1P.data(), cArr6.data());
    FiveQubitCode::applyAncillaKet(cArr_101.data(), 1, 5, ketH1N.data(), cArr6.data());

    FastPauli::kron(cArr6.data(), cArrH1.data(), cArr_11.data(), 1, 5);
    FiveQubitCode::S3Sequential(cArr6.data(), errTeleport);
    FiveQubitCode::applyAncillaKet(cArr_110.data(), 1, 5, ketH1P.data(), cArr6.data());
    FiveQubitCode::applyAncillaKet(cArr_111.data(), 1, 5, ketH1N.data(), cArr6.data());

    // ============================== Stabiliser ZXIXZ (4) ==============================
    FastPauli::kron(cArr6.data(), cArrH1.data(), cArr_000.data(), 1, 5);
    FiveQubitCode::S4Sequential(cArr6.data(), errTeleport);
    FiveQubitCode::applyAncillaKet(cArr_0000.data(), 1, 5, ketH1P.data(), cArr6.data());
    FiveQubitCode::applyAncillaKet(cArr_0001.data(), 1, 5, ketH1N.data(), cArr6.data());
    FastPauli::kron(cArr6.data(), cArrH1.data(), cArr_001.data(), 1, 5);
    FiveQubitCode::S4Sequential(cArr6.data(), errTeleport);
    FiveQubitCode::applyAncillaKet(cArr_0010.data(), 1, 5, ketH1P.data(), cArr6.data());
    FiveQubitCode::applyAncillaKet(cArr_0011.data(), 1, 5, ketH1N.data(), cArr6.data());
    FastPauli::kron(cArr6.data(), cArrH1.data(), cArr_010.data(), 1, 5);
    FiveQubitCode::S4Sequential(cArr6.data(), errTeleport);
    FiveQubitCode::applyAncillaKet(cArr_0100.data(), 1, 5, ketH1P.data(), cArr6.data());
    FiveQubitCode::applyAncillaKet(cArr_0101.data(), 1, 5, ketH1N.data(), cArr6.data());
    FastPauli::kron(cArr6.data(), cArrH1.data(), cArr_011.data(), 1, 5);
    FiveQubitCode::S4Sequential(cArr6.data(), errTeleport);
    FiveQubitCode::applyAncillaKet(cArr_0110.data(), 1, 5, ketH1P.data(), cArr6.data());
    FiveQubitCode::applyAncillaKet(cArr_0111.data(), 1, 5, ketH1N.data(), cArr6.data());
    FastPauli::kron(cArr6.data(), cArrH1.data(), cArr_100.data(), 1, 5);
    FiveQubitCode::S4Sequential(cArr6.data(), errTeleport);
    FiveQubitCode::applyAncillaKet(cArr_1000.data(), 1, 5, ketH1P.data(), cArr6.data());
    FiveQubitCode::applyAncillaKet(cArr_1001.data(), 1, 5, ketH1N.data(), cArr6.data());
    FastPauli::kron(cArr6.data(), cArrH1.data(), cArr_101.data(), 1, 5);
    FiveQubitCode::S4Sequential(cArr6.data(), errTeleport);
    FiveQubitCode::applyAncillaKet(cArr_1010.data(), 1, 5, ketH1P.data(), cArr6.data());
    FiveQubitCode::applyAncillaKet(cArr_1011.data(), 1, 5, ketH1N.data(), cArr6.data());
    FastPauli::kron(cArr6.data(), cArrH1.data(), cArr_110.data(), 1, 5);
    FiveQubitCode::S4Sequential(cArr6.data(), errTeleport);
    FiveQubitCode::applyAncillaKet(cArr_1100.data(), 1, 5, ketH1P.data(), cArr6.data());
    FiveQubitCode::applyAncillaKet(cArr_1101.data(), 1, 5, ketH1N.data(), cArr6.data());
    FastPauli::kron(cArr6.data(), cArrH1.data(), cArr_111.data(), 1, 5);
    FiveQubitCode::S4Sequential(cArr6.data(), errTeleport);
    FiveQubitCode::applyAncillaKet(cArr_1110.data(), 1, 5, ketH1P.data(), cArr6.data());
    FiveQubitCode::applyAncillaKet(cArr_1111.data(), 1, 5, ketH1N.data(), cArr6.data());

    cArrAggregate[0] = cArr_0000;
    cArrAggregate[1] = cArr_0001;
    cArrAggregate[2] = cArr_0010;
    cArrAggregate[3] = cArr_0011;
    cArrAggregate[4] = cArr_0100;
    cArrAggregate[5] = cArr_0101;
    cArrAggregate[6] = cArr_0110;
    cArrAggregate[7] = cArr_0111;
    cArrAggregate[8] = cArr_1000;
    cArrAggregate[9] = cArr_1001;
    cArrAggregate[10] = cArr_1010;
    cArrAggregate[11] = cArr_1011;
    cArrAggregate[12] = cArr_1100;
    cArrAggregate[13] = cArr_1101;
    cArrAggregate[14] = cArr_1110;
    cArrAggregate[15] = cArr_1111;

    switch (correctionMode) {
        case FiveQubitCode::correctionMode::noLossW1:
            for (uint16_t i = 0; i < UINT16_C(16); i++) {
                FiveQubitCode::applyCorrection(i, cArrAggregate[i].data());
            }
            break;
        case FiveQubitCode::correctionMode::oneLoss:
            for (uint16_t i = 0; i < UINT16_C(16); i++) {
                cArrAggregate[i] = FastPauli::getCycledDM<5>(cArrAggregate[i].data(), qubitLostIndex);
                FiveQubitCode::applyCorrection_1erasure(i, cArrAggregate[i].data(), qubitLostIndex);
            }
            break;
        case FiveQubitCode::correctionMode::noLossW2:
            for (uint16_t i = 0; i < UINT16_C(16); i++) {
                // TODO: Change `qubitLostIndex`'s name to something more general since it is shared by two different correction modes.
                FiveQubitCode::applyCorrection_flag(i, qubitLostIndex, cArrAggregate[i].data());
            }
            break;
        default:
            break;
    }

    // Sum all matrices and store the result in cArrAggregate[0].
    for (uint16_t i = 1; i < UINT16_C(16); i++) {
        std::transform(cArrAggregate[0].begin(), cArrAggregate[0].end(), cArrAggregate[i].begin(), cArrAggregate[0].begin(), std::plus<double>());  // Reference: https://stackoverflow.com/questions/3376124/how-to-add-element-by-element-of-two-stl-vectors
    }
    return cArrAggregate[0];
}

void FiveQubitCode::applyCorrection(const uint16_t i, double *cArr5) {
    constexpr int numQubits = 5;
    switch (i) {           // a1 a2 a3 a4 : Pi where P ∈ {X,Y,Z} is a single-qubit Pauli operator and i is the index of the qubit to correct with i ∈ {0,1,2,3,4};
        case UINT16_C(8):  //  0  0  0  1 : X1
            FastPauli::X_r(cArr5, numQubits, 1);
            break;
        case UINT16_C(4):  //  0  0  1  0 : Z4
            FastPauli::Z_r(cArr5, numQubits, 4);
            break;
        case UINT16_C(12):  //  0  0  1  1 : X2
            FastPauli::X_r(cArr5, numQubits, 2);
            break;
        case UINT16_C(2):  //  0  1  0  0 : Z2
            FastPauli::Z_r(cArr5, numQubits, 2);
            break;
        case UINT16_C(10):  //  0  1  0  1 : Z0
            FastPauli::Z_r(cArr5, numQubits, 0);
            break;
        case UINT16_C(6):  //  0  1  1  0 : X3
            FastPauli::X_r(cArr5, numQubits, 3);
            break;
        case UINT16_C(14):  //  0  1  1  1 : Y2
            FastPauli::Y_r(cArr5, numQubits, 2);
            break;
        case UINT16_C(1):  //  1  0  0  0 : X0
            FastPauli::X_r(cArr5, numQubits, 0);
            break;
        case UINT16_C(9):  //  1  0  0  1 : Z3
            FastPauli::Z_r(cArr5, numQubits, 3);
            break;
        case UINT16_C(5):  //  1  0  1  0 : Z1
            FastPauli::Z_r(cArr5, numQubits, 1);
            break;
        case UINT16_C(13):  //  1  0  1  1 : Y1
            FastPauli::Y_r(cArr5, numQubits, 1);
            break;
        case UINT16_C(3):  //  1  1  0  0 : X4
            FastPauli::X_r(cArr5, numQubits, 4);
            break;
        case UINT16_C(11):  //  1  1  0  1 : Y0
            FastPauli::Y_r(cArr5, numQubits, 0);
            break;
        case UINT16_C(7):  //  1  1  1  0 : Y4
            FastPauli::Y_r(cArr5, numQubits, 4);
            break;
        case UINT16_C(15):  //  1  1  1  1 : Y3
            FastPauli::Y_r(cArr5, numQubits, 3);
            break;
        default:
            break;
    }
}

std::array<std::array<double, 16>, 16> FiveQubitCode::getCorrectionArrays() {
    std::array<std::array<double, 16>, 16> myAncillaKets;
    std::array<double, 16> ketH4 = FiveQubitCode::getKetH<4>();  // ket(+,+,+,+)
    std::array<double, 16> tempKetH4;
    for (uint16_t i = 0; i < 16Ui16; i++) {
        tempKetH4 = ketH4;
        for (uint16_t j = 0; j < 4Ui16; j++) {
            if (((i >> (3 - j)) & 0x01)) {
                FastPauli::ketZ_r(tempKetH4.data(), 4, j);
            }
        }
        myAncillaKets[i] = tempKetH4;
    }
    return myAncillaKets;
}

// double FiveQubitCode::applyCorrectionArrays(std::array<std::array<double, 512>, 16> const &correctionArrays, const double *cArr9) {
//     // To-do: Make this function also return the resulting 5-qubit matrix so that we can use it to calculate the 2nd order fidelity using the recursion approximation.
//     double fidel = 0;
//     for (uint16_t i = 0; i < 16Ui16; i++) {
//         fidel += FastPauli::fidel_r(correctionArrays[i].data(), cArr9, 9);
//     }
//     return fidel;
// }

void FiveQubitCode::applyAncillaKet(double *cArr5, const int numAncillaQubits, const int numQubits, const double *ancillaKet, const double *cArr9) {
    const int nA = 1 << numAncillaQubits;
    const int nQ = 1 << numQubits;
    const int n = 1 << (numQubits + numAncillaQubits);
    int cIndex, cIndex2, cIndex3, cIndex4;
    double sum;
    for (int p = 0; p < nQ; p++) {
        for (int q = p; q < nQ; q++) {
            cIndex = FastPauli::getRow(p, q, nQ);
            sum = 0;
            for (int i = 0; i < nA; i++) {
                cIndex2 = FastPauli::getRow(i * nQ + p, i * nQ + q, n);
                sum += cArr9[cIndex2] * ancillaKet[i] * ancillaKet[i];
                for (int j = i + 1; j < nA; j++) {
                    cIndex3 = FastPauli::getRow(i * nQ + p, j * nQ + q, n);
                    cIndex4 = FastPauli::getRow(i * nQ + q, j * nQ + p, n);
                    sum += (cArr9[cIndex3] + cArr9[cIndex4]) * ancillaKet[i] * ancillaKet[j];
                }
            }
            cArr5[cIndex] = sum;
        }
    }
}

FastPauli::CArr<5> FiveQubitCode::applyCorrectionArrays(const double *cArr9, std::array<std::array<double, 16>, 16> const &ancillaKets) {
    constexpr int numAncillaQubits = 4;
    constexpr int numQubits = 5;
    constexpr uint64_t cArrLen = FastPauli::CArrLen<numQubits>;
    std::array<double, cArrLen> dm;
    dm.fill(0.0);
    std::array<double, cArrLen> tempDm;

    for (uint16_t a = 0; a < 16; a++) {
        FiveQubitCode::applyAncillaKet(tempDm.data(), numAncillaQubits, numQubits, ancillaKets[a].data(), cArr9);
        // Apply Pauli correction to tempDm.
        FiveQubitCode::applyCorrection(a, tempDm.data());

        // Store corrected tempDm in dm.
        for (uint64_t k = 0; k < cArrLen; k++) {
            dm[k] += tempDm[k];
        }
    }

    // Return the corrected 5-qubit density matrix (C-array form).
    return dm;
}