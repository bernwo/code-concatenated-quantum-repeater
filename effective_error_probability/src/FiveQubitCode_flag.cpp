#include "FiveQubitCode.hpp"

void FiveQubitCode::S1_flag(double *cArray, const double err, const double errTeleport, const double ancillaFlagErr) {
    constexpr int numQubits = 7;
    NoisyPauli::NCX_r(cArray, numQubits, 0, 2, err);
    NoisyPauli::NCX_r(cArray, numQubits, 0, 1, ancillaFlagErr);
    NoisyPauli::NCZ_r(cArray, numQubits, 0, 3, errTeleport);
    NoisyPauli::NCZ_r(cArray, numQubits, 0, 4, errTeleport);
    NoisyPauli::NCX_r(cArray, numQubits, 0, 1, ancillaFlagErr);
    NoisyPauli::NCX_r(cArray, numQubits, 0, 5, errTeleport);
}

void FiveQubitCode::S2_flag(double *cArray, const double err, const double errTeleport, const double ancillaFlagErr) {
    constexpr int numQubits = 7;
    NoisyPauli::NCX_r(cArray, numQubits, 0, 3, err);
    NoisyPauli::NCX_r(cArray, numQubits, 0, 1, ancillaFlagErr);
    NoisyPauli::NCZ_r(cArray, numQubits, 0, 4, errTeleport);
    NoisyPauli::NCZ_r(cArray, numQubits, 0, 5, errTeleport);
    NoisyPauli::NCX_r(cArray, numQubits, 0, 1, ancillaFlagErr);
    NoisyPauli::NCX_r(cArray, numQubits, 0, 6, errTeleport);
}

void FiveQubitCode::S3_flag(double *cArray, const double err, const double errTeleport, const double ancillaFlagErr) {
    constexpr int numQubits = 7;
    NoisyPauli::NCX_r(cArray, numQubits, 0, 2, errTeleport);
    NoisyPauli::NCX_r(cArray, numQubits, 0, 1, ancillaFlagErr);
    NoisyPauli::NCX_r(cArray, numQubits, 0, 4, err);
    NoisyPauli::NCZ_r(cArray, numQubits, 0, 5, errTeleport);
    NoisyPauli::NCX_r(cArray, numQubits, 0, 1, ancillaFlagErr);
    NoisyPauli::NCZ_r(cArray, numQubits, 0, 6, errTeleport);
}

void FiveQubitCode::S4_flag(double *cArray, const double err, const double errTeleport, const double ancillaFlagErr) {
    constexpr int numQubits = 7;
    NoisyPauli::NCZ_r(cArray, numQubits, 0, 2, errTeleport);
    NoisyPauli::NCX_r(cArray, numQubits, 0, 1, ancillaFlagErr);
    NoisyPauli::NCX_r(cArray, numQubits, 0, 3, errTeleport);
    NoisyPauli::NCX_r(cArray, numQubits, 0, 5, err);
    NoisyPauli::NCX_r(cArray, numQubits, 0, 1, ancillaFlagErr);
    NoisyPauli::NCZ_r(cArray, numQubits, 0, 6, errTeleport);
}

FastPauli::CArr<5> FiveQubitCode::getDensityMatrixW0(const double *cArray) {
    constexpr std::array<double, 4> Ket_AF_00 = FiveQubitCode::getAncillaFlagKet<false, false>();
    FastPauli::CArr<5> result = FiveQubitCode::applyAncillaKet<2, 5>(Ket_AF_00.data(), cArray);
    return result;
}

FastPauli::CArr<5> FiveQubitCode::getDensityMatrixW1(const double *cArray) {
    constexpr std::array<double, 4> Ket_AF_10 = FiveQubitCode::getAncillaFlagKet<true, false>();
    FastPauli::CArr<5> result = FiveQubitCode::applyAncillaKet<2, 5>(Ket_AF_10.data(), cArray);
    return result;
}

FastPauli::CArr<5> FiveQubitCode::getDensityMatrixW2(const double *cArray) {
    constexpr std::array<double, 4> Ket_AF_01 = FiveQubitCode::getAncillaFlagKet<false, true>();
    constexpr std::array<double, 4> Ket_AF_11 = FiveQubitCode::getAncillaFlagKet<true, true>();
    FastPauli::CArr<5> tmp = FiveQubitCode::applyAncillaKet<2, 5>(Ket_AF_11.data(), cArray);
    FastPauli::CArr<5> result = FiveQubitCode::applyAncillaKet<2, 5>(Ket_AF_01.data(), cArray);
    std::transform(result.begin(), result.end(), tmp.begin(), result.begin(), std::plus<double>());
    return result;
}

FastPauli::CArr<5> FiveQubitCode::OneRoundSequential_flag(FiveQubitCode::ECInfo &ecInfo, double *logicalKet, double *cArrCodeword, const double err, const double errTeleport, const double ancillaFlagErr, const int ecInfoIndex) {
    // Initialise the ancilla-flag pair density matrix. Only the untriggered one is needed.
    constexpr std::array<double, 10> AF_00 = FiveQubitCode::getAncillaFlagDensityMatrix<false, false>();

    // Initialise temporary cArrays for performing intermediary calculations.
    FastPauli::CArr<7> tmp;
    FastPauli::CArr<5> cArrW0, cArrW1;
    FastPauli::CArr<5> cArrW00, cArrW01;
    FastPauli::CArr<5> cArrW000, cArrW001;
    FastPauli::CArr<5> cArrW0000, cArrW0001;
    FastPauli::CArr<5> cArrW2, cArrW02, cArrW002, cArrW0002;

    // ========================== Apply XZZXI stabiliser (1) ==========================
    tmp = FastPauli::kron<2, 5>(AF_00.data(), cArrCodeword);
    FiveQubitCode::S1_flag(tmp.data(), err, errTeleport, ancillaFlagErr);
    cArrW0 = FiveQubitCode::getDensityMatrixW0(tmp.data());
    cArrW1 = FiveQubitCode::getDensityMatrixW1(tmp.data());
    cArrW2 = FiveQubitCode::getDensityMatrixW2(tmp.data());

    cArrW1 = FiveQubitCode::OneRoundSequential(cArrW1.data(), err, errTeleport, FiveQubitCode::correctionMode::noLossW1, 0);
    cArrW2 = FiveQubitCode::OneRoundSequential(cArrW2.data(), err, errTeleport, FiveQubitCode::correctionMode::noLossW2, 0);

    // ========================== Apply IXZZX stabiliser (2) ==========================
    tmp = FastPauli::kron<2, 5>(AF_00.data(), cArrW0.data());
    FiveQubitCode::S2_flag(tmp.data(), err, errTeleport, ancillaFlagErr);
    cArrW00 = FiveQubitCode::getDensityMatrixW0(tmp.data());
    cArrW01 = FiveQubitCode::getDensityMatrixW1(tmp.data());
    cArrW02 = FiveQubitCode::getDensityMatrixW2(tmp.data());

    cArrW01 = FiveQubitCode::OneRoundSequential(cArrW01.data(), err, errTeleport, FiveQubitCode::correctionMode::noLossW1, 0);
    cArrW02 = FiveQubitCode::OneRoundSequential(cArrW02.data(), err, errTeleport, FiveQubitCode::correctionMode::noLossW2, 1);

    // ========================== Apply XIXZZ stabiliser (3) ==========================
    tmp = FastPauli::kron<2, 5>(AF_00.data(), cArrW00.data());
    FiveQubitCode::S3_flag(tmp.data(), err, errTeleport, ancillaFlagErr);
    cArrW000 = FiveQubitCode::getDensityMatrixW0(tmp.data());
    cArrW001 = FiveQubitCode::getDensityMatrixW1(tmp.data());
    cArrW002 = FiveQubitCode::getDensityMatrixW2(tmp.data());

    cArrW001 = FiveQubitCode::OneRoundSequential(cArrW001.data(), err, errTeleport, FiveQubitCode::correctionMode::noLossW1, 0);
    cArrW002 = FiveQubitCode::OneRoundSequential(cArrW002.data(), err, errTeleport, FiveQubitCode::correctionMode::noLossW2, 2);

    // ========================== Apply ZXIXZ stabiliser (4) ==========================
    tmp = FastPauli::kron<2, 5>(AF_00.data(), cArrW000.data());
    FiveQubitCode::S4_flag(tmp.data(), err, errTeleport, ancillaFlagErr);
    cArrW0000 = FiveQubitCode::getDensityMatrixW0(tmp.data());
    cArrW0001 = FiveQubitCode::getDensityMatrixW1(tmp.data());
    cArrW0002 = FiveQubitCode::getDensityMatrixW2(tmp.data());
    cArrW0001 = FiveQubitCode::OneRoundSequential(cArrW0001.data(), err, errTeleport, FiveQubitCode::correctionMode::noLossW1, 0);
    cArrW0002 = FiveQubitCode::OneRoundSequential(cArrW0002.data(), err, errTeleport, FiveQubitCode::correctionMode::noLossW2, 3);

    // Aggregate all corrected matrices. Use cArrW2 to store the results to avoid declaring new std::array.
    std::transform(cArrW2.begin(), cArrW2.end(), cArrW02.begin(), cArrW2.begin(), std::plus<double>());
    std::transform(cArrW2.begin(), cArrW2.end(), cArrW002.begin(), cArrW2.begin(), std::plus<double>());
    std::transform(cArrW2.begin(), cArrW2.end(), cArrW0002.begin(), cArrW2.begin(), std::plus<double>());
    std::transform(cArrW2.begin(), cArrW2.end(), cArrW1.begin(), cArrW2.begin(), std::plus<double>());
    std::transform(cArrW2.begin(), cArrW2.end(), cArrW01.begin(), cArrW2.begin(), std::plus<double>());
    std::transform(cArrW2.begin(), cArrW2.end(), cArrW001.begin(), cArrW2.begin(), std::plus<double>());
    std::transform(cArrW2.begin(), cArrW2.end(), cArrW0001.begin(), cArrW2.begin(), std::plus<double>());
    std::transform(cArrW2.begin(), cArrW2.end(), cArrW0000.begin(), cArrW2.begin(), std::plus<double>());

    const double fidel = FastPauli::fidel_r(logicalKet, cArrW2.data(), 5);
    const double trace = FastPauli::trace_r(cArrW2.data(), 5);
    ecInfo.grandErrors[ecInfoIndex] = 1 - fidel;
    ecInfo.grandTraces[ecInfoIndex] = trace;

    // Finish and return the aggregated 5-qubit density matrix.
    return cArrW2;
}

FastPauli::CArr<5> FiveQubitCode::OneRoundSequential_flag(double *cArrCodeword, const double err, const double errTeleport, const double ancillaFlagErr) {
    // Initialise the ancilla-flag pair density matrix. Only the untriggered one is needed.
    constexpr std::array<double, 10> AF_00 = FiveQubitCode::getAncillaFlagDensityMatrix<false, false>();

    // Initialise temporary cArrays for performing intermediary calculations.
    FastPauli::CArr<7> tmp;
    FastPauli::CArr<5> cArrW0, cArrW1;
    FastPauli::CArr<5> cArrW00, cArrW01;
    FastPauli::CArr<5> cArrW000, cArrW001;
    FastPauli::CArr<5> cArrW0000, cArrW0001;
    FastPauli::CArr<5> cArrW2, cArrW02, cArrW002, cArrW0002;

    // ========================== Apply XZZXI stabiliser (1) ==========================
    tmp = FastPauli::kron<2, 5>(AF_00.data(), cArrCodeword);
    FiveQubitCode::S1_flag(tmp.data(), err, errTeleport, ancillaFlagErr);
    cArrW0 = FiveQubitCode::getDensityMatrixW0(tmp.data());
    cArrW1 = FiveQubitCode::getDensityMatrixW1(tmp.data());
    cArrW2 = FiveQubitCode::getDensityMatrixW2(tmp.data());

    cArrW1 = FiveQubitCode::OneRoundSequential(cArrW1.data(), err, errTeleport, FiveQubitCode::correctionMode::noLossW1, 0);
    cArrW2 = FiveQubitCode::OneRoundSequential(cArrW2.data(), err, errTeleport, FiveQubitCode::correctionMode::noLossW2, 0);

    // ========================== Apply IXZZX stabiliser (2) ==========================
    tmp = FastPauli::kron<2, 5>(AF_00.data(), cArrW0.data());
    FiveQubitCode::S2_flag(tmp.data(), err, errTeleport, ancillaFlagErr);
    cArrW00 = FiveQubitCode::getDensityMatrixW0(tmp.data());
    cArrW01 = FiveQubitCode::getDensityMatrixW1(tmp.data());
    cArrW02 = FiveQubitCode::getDensityMatrixW2(tmp.data());

    cArrW01 = FiveQubitCode::OneRoundSequential(cArrW01.data(), err, errTeleport, FiveQubitCode::correctionMode::noLossW1, 0);
    cArrW02 = FiveQubitCode::OneRoundSequential(cArrW02.data(), err, errTeleport, FiveQubitCode::correctionMode::noLossW2, 1);

    // ========================== Apply XIXZZ stabiliser (3) ==========================
    tmp = FastPauli::kron<2, 5>(AF_00.data(), cArrW00.data());
    FiveQubitCode::S3_flag(tmp.data(), err, errTeleport, ancillaFlagErr);
    cArrW000 = FiveQubitCode::getDensityMatrixW0(tmp.data());
    cArrW001 = FiveQubitCode::getDensityMatrixW1(tmp.data());
    cArrW002 = FiveQubitCode::getDensityMatrixW2(tmp.data());

    cArrW001 = FiveQubitCode::OneRoundSequential(cArrW001.data(), err, errTeleport, FiveQubitCode::correctionMode::noLossW1, 0);
    cArrW002 = FiveQubitCode::OneRoundSequential(cArrW002.data(), err, errTeleport, FiveQubitCode::correctionMode::noLossW2, 2);

    // ========================== Apply ZXIXZ stabiliser (4) ==========================
    tmp = FastPauli::kron<2, 5>(AF_00.data(), cArrW000.data());
    FiveQubitCode::S4_flag(tmp.data(), err, errTeleport, ancillaFlagErr);
    cArrW0000 = FiveQubitCode::getDensityMatrixW0(tmp.data());
    cArrW0001 = FiveQubitCode::getDensityMatrixW1(tmp.data());
    cArrW0002 = FiveQubitCode::getDensityMatrixW2(tmp.data());
    cArrW0001 = FiveQubitCode::OneRoundSequential(cArrW0001.data(), err, errTeleport, FiveQubitCode::correctionMode::noLossW1, 0);
    cArrW0002 = FiveQubitCode::OneRoundSequential(cArrW0002.data(), err, errTeleport, FiveQubitCode::correctionMode::noLossW2, 3);

    // Aggregate all corrected matrices. Use cArrW2 to store the results to avoid declaring new std::array.
    std::transform(cArrW2.begin(), cArrW2.end(), cArrW02.begin(), cArrW2.begin(), std::plus<double>());
    std::transform(cArrW2.begin(), cArrW2.end(), cArrW002.begin(), cArrW2.begin(), std::plus<double>());
    std::transform(cArrW2.begin(), cArrW2.end(), cArrW0002.begin(), cArrW2.begin(), std::plus<double>());
    std::transform(cArrW2.begin(), cArrW2.end(), cArrW1.begin(), cArrW2.begin(), std::plus<double>());
    std::transform(cArrW2.begin(), cArrW2.end(), cArrW01.begin(), cArrW2.begin(), std::plus<double>());
    std::transform(cArrW2.begin(), cArrW2.end(), cArrW001.begin(), cArrW2.begin(), std::plus<double>());
    std::transform(cArrW2.begin(), cArrW2.end(), cArrW0001.begin(), cArrW2.begin(), std::plus<double>());
    std::transform(cArrW2.begin(), cArrW2.end(), cArrW0000.begin(), cArrW2.begin(), std::plus<double>());

    // Finish and return the aggregated 5-qubit density matrix.
    return cArrW2;
}

FastPauli::CArr<5> FiveQubitCode::OneRound_flag(FiveQubitCode::ECInfo &ecInfo, double *logicalKet, double *cArrCodeword, const double err, const double errTeleport, const double ancillaFlagErr, std::array<std::array<double, 16>, 16> const &ancillaKets, const int ecInfoIndex) {
    // void FiveQubitCode::OneRound_flag(double *logicalKet, double *cArrCodeword, const double err, const double errTeleport, const double ancillaFlagErr, std::array<std::array<double, 16>, 16> const &ancillaKets) {
    double fidelity_cum_tmp = 0;
    double fidelity_cum = 0;
    double fidelity_tmp = 0;
    double trace_tmp = 0;

    // Initialise the density matrix for a single ancilla qubit.
    std::array<double, FastPauli::CArrLen<1>> cArrH1 = FiveQubitCode::getH<1>();

    // Initialise the ancilla-flag pair density matrix. Only the untriggered one is needed.
    constexpr std::array<double, 10> AF_00 = FiveQubitCode::getAncillaFlagDensityMatrix<false, false>();

    // Initialise the ancilla-flag pair ket vectors.
    constexpr std::array<double, 4> Ket_AF_00 = FiveQubitCode::getAncillaFlagKet<false, false>();
    constexpr std::array<double, 4> Ket_AF_01 = FiveQubitCode::getAncillaFlagKet<false, true>();
    constexpr std::array<double, 4> Ket_AF_10 = FiveQubitCode::getAncillaFlagKet<true, false>();
    constexpr std::array<double, 4> Ket_AF_11 = FiveQubitCode::getAncillaFlagKet<true, true>();

    // Initiailise 5-qubit density matrix std::array's in which we store the results after partial tracing.
    FastPauli::CArr<5> res_S1_00;
    FastPauli::CArr<5> res_S1_10;
    FastPauli::CArr<5> res_S1_01;
    FastPauli::CArr<5> res_S1_11;
    FastPauli::CArr<5> res_S2_00;
    FastPauli::CArr<5> res_S2_10;
    FastPauli::CArr<5> res_S2_01;
    FastPauli::CArr<5> res_S2_11;
    FastPauli::CArr<5> res_S3_00;
    FastPauli::CArr<5> res_S3_10;
    FastPauli::CArr<5> res_S3_01;
    FastPauli::CArr<5> res_S3_11;
    FastPauli::CArr<5> res_S4_00;
    FastPauli::CArr<5> res_S4_10;
    FastPauli::CArr<5> res_S4_01;
    FastPauli::CArr<5> res_S4_11;

    // Initialise temporary arrays.
    FastPauli::CArr<5> cArr5_cum{};  // Stores result.
    FastPauli::CArr<5> cArr5;        // Stores result.
    FastPauli::CArr<6> cArr6;
    FastPauli::CArr<7> cArr7;
    FastPauli::CArr<8> cArr8;
    FastPauli::CArr<9> cArr9{};

    // Prepare the 7-qubit density matrix (1 ancilla, 1 flag, 5 data qubits) with fresh ancilla-flag pair (i.e., not triggered).
    std::array<double, FastPauli::CArrLen<7>> S1_00;
    FastPauli::kron(S1_00.data(), AF_00.data(), cArrCodeword, 2, 5);

    // Perform the `XZZXI` fault-tolerant stabiliser operation, which is 1 out of 4 sets of stabiliser operations.
    // ===========================================================================================================
    FiveQubitCode::S1_flag(S1_00.data(), err, errTeleport, ancillaFlagErr);

    // Get the truncated density matrix (5-qubits) by tracing away the ancilla-flag pair for all possible values of the ancilla-flag pair.
    FiveQubitCode::applyAncillaKet(res_S1_00.data(), 2, 5, Ket_AF_00.data(), S1_00.data());                         // This gets to continue to the next fault tolerant stabiliser operation set.
    FiveQubitCode::applyAncillaKet(res_S1_01.data(), 2, 5, Ket_AF_01.data(), S1_00.data());                         // Stop here. Go to weight<=2 correction.
    FiveQubitCode::applyAncillaKet(res_S1_10.data(), 2, 5, Ket_AF_10.data(), S1_00.data());                         // Stop here. Go to weight<=1 correction.
    FiveQubitCode::applyAncillaKet(res_S1_11.data(), 2, 5, Ket_AF_11.data(), S1_00.data());                         // Stop here. Go to weight<=2 correction.
    std::transform(res_S1_11.begin(), res_S1_11.end(), res_S1_01.begin(), res_S1_11.begin(), std::plus<double>());  // Writes the elementwise sum back to `res_S1_11`.

    // (res_S1_10) operations - ancilla triggered.
    // -------------------------------------------
    // Apply the unflagged stabiliser operations on `res_S1_10`.
    FiveQubitCode::OneRound(cArr6.data(), cArr7.data(), cArr8.data(), cArr9.data(), res_S1_10.data(), cArrH1.data(), err, errTeleport);

    // Perform the weight<=1 correction on `res_S1_10`.
    cArr5 = FiveQubitCode::applyCorrectionArrays(cArr9.data(), ancillaKets);
    std::transform(cArr5_cum.begin(), cArr5_cum.end(), cArr5.begin(), cArr5_cum.begin(), std::plus<double>());

    // std::cout << "cArr5_cum: " << cArr5_cum[34] << std::endl;
    // Store fidelity.
    fidelity_tmp = FastPauli::fidel_r(logicalKet, cArr5.data(), 5);
    fidelity_cum_tmp += fidelity_tmp;
    fidelity_cum += fidelity_tmp;

    // Store trace.
    trace_tmp += FastPauli::trace_r(cArr5.data(), 5);

    // (res_S1_01 + res_S1_11) operations - flag triggered.
    // ----------------------------------------------------
    // Apply the unflagged stabiliser operations on `res_S1_10+res_S1_11`.
    FiveQubitCode::OneRound(cArr6.data(), cArr7.data(), cArr8.data(), cArr9.data(), res_S1_11.data(), cArrH1.data(), err, errTeleport);

    // Perform the weight<=2 correction on `res_S1_10+res_S1_11`.
    cArr5 = FiveQubitCode::applyCorrectionArrays_flag(cArr9.data(), ancillaKets, 3);  // flag = 3 means flag number 4 in Bernard's thesis.
    std::transform(cArr5_cum.begin(), cArr5_cum.end(), cArr5.begin(), cArr5_cum.begin(), std::plus<double>());

    // Store fidelity.
    fidelity_tmp = FastPauli::fidel_r(logicalKet, cArr5.data(), 5);
    fidelity_cum_tmp += fidelity_tmp;
    fidelity_cum += fidelity_tmp;
    // Store trace.
    trace_tmp += FastPauli::trace_r(cArr5.data(), 5);

    // Finish `XZZXI` fault-tolerant stabiliser operation.
    // ---------------------------------------------------
    if (ecInfoIndex >= 0) {
        ecInfo.errors[0][ecInfoIndex] = 1 - fidelity_cum_tmp;
        ecInfo.traces[0][ecInfoIndex] = trace_tmp;
    }

    trace_tmp = 0;
    fidelity_cum_tmp = 0;

    // Prepare the 7-qubit density matrix (1 ancilla, 1 flag, 5 data qubits) with fresh ancilla-flag pair (i.e., not triggered).
    std::array<double, FastPauli::CArrLen<7>> S2_00;
    FastPauli::kron(S2_00.data(), AF_00.data(), res_S1_00.data(), 2, 5);

    // Perform the `IXZZX` fault-tolerant stabiliser operation, which is 2 out of 4 sets of stabiliser operations.
    // ===========================================================================================================
    FiveQubitCode::S2_flag(S2_00.data(), err, errTeleport, ancillaFlagErr);

    // Get the truncated density matrix (5-qubits) by tracing away the ancilla-flag pair for all possible values of the ancilla-flag pair.
    FiveQubitCode::applyAncillaKet(res_S2_00.data(), 2, 5, Ket_AF_00.data(), S2_00.data());                         // This gets to continue to the next fault tolerant stabiliser operation set.
    FiveQubitCode::applyAncillaKet(res_S2_01.data(), 2, 5, Ket_AF_01.data(), S2_00.data());                         // Stop here. Go to weight<=2 correction.
    FiveQubitCode::applyAncillaKet(res_S2_10.data(), 2, 5, Ket_AF_10.data(), S2_00.data());                         // Stop here. Go to weight<=1 correction.
    FiveQubitCode::applyAncillaKet(res_S2_11.data(), 2, 5, Ket_AF_11.data(), S2_00.data());                         // Stop here. Go to weight<=2 correction.
    std::transform(res_S2_11.begin(), res_S2_11.end(), res_S2_01.begin(), res_S2_11.begin(), std::plus<double>());  // Writes the elementwise sum back to `res_S2_11`.

    // (res_S2_10) operations - ancilla triggered.
    // -------------------------------------------
    // Apply the unflagged stabiliser operations on `res_S2_10`.
    FiveQubitCode::OneRound(cArr6.data(), cArr7.data(), cArr8.data(), cArr9.data(), res_S2_10.data(), cArrH1.data(), err, errTeleport);

    // Perform the weight<=1 correction on `res_S2_10`.
    cArr5 = FiveQubitCode::applyCorrectionArrays(cArr9.data(), ancillaKets);
    std::transform(cArr5_cum.begin(), cArr5_cum.end(), cArr5.begin(), cArr5_cum.begin(), std::plus<double>());

    // Store fidelity.
    fidelity_tmp = FastPauli::fidel_r(logicalKet, cArr5.data(), 5);
    fidelity_cum_tmp += fidelity_tmp;
    fidelity_cum += fidelity_tmp;
    // Store trace.
    trace_tmp += FastPauli::trace_r(cArr5.data(), 5);

    // (res_S2_01 + res_S2_11) operations - flag triggered.
    // ----------------------------------------------------
    // Apply the unflagged stabiliser operations on `res_S2_10+res_S2_11`.
    FiveQubitCode::OneRound(cArr6.data(), cArr7.data(), cArr8.data(), cArr9.data(), res_S2_11.data(), cArrH1.data(), err, errTeleport);

    // Perform the weight<=2 correction on `res_S2_10+res_S2_11`.
    cArr5 = FiveQubitCode::applyCorrectionArrays_flag(cArr9.data(), ancillaKets, 2);  // flag = 2 means flag number 3 in Bernard's thesis.
    std::transform(cArr5_cum.begin(), cArr5_cum.end(), cArr5.begin(), cArr5_cum.begin(), std::plus<double>());

    // Store fidelity.
    fidelity_tmp = FastPauli::fidel_r(logicalKet, cArr5.data(), 5);
    fidelity_cum_tmp += fidelity_tmp;
    fidelity_cum += fidelity_tmp;
    // Store trace.
    trace_tmp += FastPauli::trace_r(cArr5.data(), 5);

    // Finish `IXZZX` fault-tolerant stabiliser operation.
    // ---------------------------------------------------
    if (ecInfoIndex >= 0) {
        ecInfo.errors[1][ecInfoIndex] = 1 - fidelity_cum_tmp;
        ecInfo.traces[1][ecInfoIndex] = trace_tmp;
    }
    trace_tmp = 0;
    fidelity_cum_tmp = 0;

    // Prepare the 7-qubit density matrix (1 ancilla, 1 flag, 5 data qubits) with fresh ancilla-flag pair (i.e., not triggered).
    std::array<double, FastPauli::CArrLen<7>> S3_00;
    FastPauli::kron(S3_00.data(), AF_00.data(), res_S2_00.data(), 2, 5);

    // Perform the `XIXZZ` fault-tolerant stabiliser operation, which is 3 out of 4 sets of stabiliser operations.
    // ===========================================================================================================
    FiveQubitCode::S3_flag(S3_00.data(), err, errTeleport, ancillaFlagErr);

    // Get the truncated density matrix (5-qubits) by tracing away the ancilla-flag pair for all possible values of the ancilla-flag pair.
    FiveQubitCode::applyAncillaKet(res_S3_00.data(), 2, 5, Ket_AF_00.data(), S3_00.data());                         // This gets to continue to the next fault tolerant stabiliser operation set.
    FiveQubitCode::applyAncillaKet(res_S3_01.data(), 2, 5, Ket_AF_01.data(), S3_00.data());                         // Stop here. Go to weight<=2 correction.
    FiveQubitCode::applyAncillaKet(res_S3_10.data(), 2, 5, Ket_AF_10.data(), S3_00.data());                         // Stop here. Go to weight<=1 correction.
    FiveQubitCode::applyAncillaKet(res_S3_11.data(), 2, 5, Ket_AF_11.data(), S3_00.data());                         // Stop here. Go to weight<=2 correction.
    std::transform(res_S3_11.begin(), res_S3_11.end(), res_S3_01.begin(), res_S3_11.begin(), std::plus<double>());  // Writes the elementwise sum back to `res_S3_11`.

    // (res_S3_10) operations - ancilla triggered.
    // -------------------------------------------
    // Apply the unflagged stabiliser operations on `res_S3_10`.
    FiveQubitCode::OneRound(cArr6.data(), cArr7.data(), cArr8.data(), cArr9.data(), res_S3_10.data(), cArrH1.data(), err, errTeleport);

    // Perform the weight<=1 correction on `res_S3_10`.
    cArr5 = FiveQubitCode::applyCorrectionArrays(cArr9.data(), ancillaKets);
    std::transform(cArr5_cum.begin(), cArr5_cum.end(), cArr5.begin(), cArr5_cum.begin(), std::plus<double>());

    // Store fidelity.
    fidelity_tmp = FastPauli::fidel_r(logicalKet, cArr5.data(), 5);
    fidelity_cum_tmp += fidelity_tmp;
    fidelity_cum += fidelity_tmp;
    // Store trace.
    trace_tmp += FastPauli::trace_r(cArr5.data(), 5);

    // (res_S3_01 + res_S3_11) operations - flag triggered.
    // ----------------------------------------------------
    // Apply the unflagged stabiliser operations on `res_S3_10+res_S3_11`.
    FiveQubitCode::OneRound(cArr6.data(), cArr7.data(), cArr8.data(), cArr9.data(), res_S3_11.data(), cArrH1.data(), err, errTeleport);

    // Perform the weight<=2 correction on `res_S3_10+res_S3_11`.
    cArr5 = FiveQubitCode::applyCorrectionArrays_flag(cArr9.data(), ancillaKets, 1);  // flag = 1 means flag number 2 in Bernard's thesis.
    std::transform(cArr5_cum.begin(), cArr5_cum.end(), cArr5.begin(), cArr5_cum.begin(), std::plus<double>());

    // Store fidelity.
    fidelity_tmp = FastPauli::fidel_r(logicalKet, cArr5.data(), 5);
    fidelity_cum_tmp += fidelity_tmp;
    fidelity_cum += fidelity_tmp;
    // Store trace.
    trace_tmp += FastPauli::trace_r(cArr5.data(), 5);

    // Finish `XIXZZ` fault-tolerant stabiliser operation.
    // ---------------------------------------------------
    if (ecInfoIndex >= 0) {
        ecInfo.errors[2][ecInfoIndex] = 1 - fidelity_cum_tmp;
        ecInfo.traces[2][ecInfoIndex] = trace_tmp;
    }
    trace_tmp = 0;
    fidelity_cum_tmp = 0;

    // Prepare the 7-qubit density matrix (1 ancilla, 1 flag, 5 data qubits) with fresh ancilla-flag pair (i.e., not triggered).
    std::array<double, FastPauli::CArrLen<7>> S4_00;
    FastPauli::kron(S4_00.data(), AF_00.data(), res_S3_00.data(), 2, 5);

    // Perform the `ZXIXZ` fault-tolerant stabiliser operation, which is 4 out of 4 sets of stabiliser operations.
    // ===========================================================================================================
    FiveQubitCode::S4_flag(S4_00.data(), err, errTeleport, ancillaFlagErr);

    // Get the truncated density matrix (5-qubits) by tracing away the ancilla-flag pair for all possible values of the ancilla-flag pair.
    FiveQubitCode::applyAncillaKet(res_S4_00.data(), 2, 5, Ket_AF_00.data(), S4_00.data());                         // This gets to continue to the next fault tolerant stabiliser operation set.
    FiveQubitCode::applyAncillaKet(res_S4_01.data(), 2, 5, Ket_AF_01.data(), S4_00.data());                         // Stop here. Go to weight<=2 correction.
    FiveQubitCode::applyAncillaKet(res_S4_10.data(), 2, 5, Ket_AF_10.data(), S4_00.data());                         // Stop here. Go to weight<=1 correction.
    FiveQubitCode::applyAncillaKet(res_S4_11.data(), 2, 5, Ket_AF_11.data(), S4_00.data());                         // Stop here. Go to weight<=2 correction.
    std::transform(res_S4_11.begin(), res_S4_11.end(), res_S4_01.begin(), res_S4_11.begin(), std::plus<double>());  // Writes the elementwise sum back to `res_S4_11`.

    // (res_S4_10) operations - ancilla triggered.
    // -------------------------------------------
    // Apply the unflagged stabiliser operations on `res_S4_10`.
    FiveQubitCode::OneRound(cArr6.data(), cArr7.data(), cArr8.data(), cArr9.data(), res_S4_10.data(), cArrH1.data(), err, errTeleport);

    // Perform the weight<=1 correction on `res_S4_10`.
    cArr5 = FiveQubitCode::applyCorrectionArrays(cArr9.data(), ancillaKets);
    std::transform(cArr5_cum.begin(), cArr5_cum.end(), cArr5.begin(), cArr5_cum.begin(), std::plus<double>());

    // Store fidelity.
    fidelity_tmp = FastPauli::fidel_r(logicalKet, cArr5.data(), 5);
    fidelity_cum_tmp += fidelity_tmp;
    fidelity_cum += fidelity_tmp;
    // Store trace.
    trace_tmp += FastPauli::trace_r(cArr5.data(), 5);

    // (res_S4_01 + res_S4_11) operations - flag triggered.
    // ----------------------------------------------------
    // Apply the unflagged stabiliser operations on `res_S4_10+res_S4_11`.
    FiveQubitCode::OneRound(cArr6.data(), cArr7.data(), cArr8.data(), cArr9.data(), res_S4_11.data(), cArrH1.data(), err, errTeleport);

    // Perform the weight<=2 correction on `res_S4_10+res_S4_11`.
    cArr5 = FiveQubitCode::applyCorrectionArrays_flag(cArr9.data(), ancillaKets, 0);  // flag = 0 means flag number 1 in Bernard's thesis.
    std::transform(cArr5_cum.begin(), cArr5_cum.end(), cArr5.begin(), cArr5_cum.begin(), std::plus<double>());

    // Store fidelity.
    fidelity_tmp = FastPauli::fidel_r(logicalKet, cArr5.data(), 5);
    fidelity_cum_tmp += fidelity_tmp;
    fidelity_cum += fidelity_tmp;
    // Store trace.
    trace_tmp += FastPauli::trace_r(cArr5.data(), 5);

    // Finish `ZXIXZ` fault-tolerant stabiliser operation.
    // ---------------------------------------------------
    if (ecInfoIndex >= 0) {
        ecInfo.errors[3][ecInfoIndex] = 1 - fidelity_cum_tmp;
        ecInfo.traces[3][ecInfoIndex] = trace_tmp;
    }
    trace_tmp = 0;
    fidelity_cum_tmp = 0;

    // If no ancillas or flags were triggered at all.
    // ==============================================
    std::transform(cArr5_cum.begin(), cArr5_cum.end(), res_S4_00.begin(), cArr5_cum.begin(), std::plus<double>());
    fidelity_cum_tmp = FastPauli::fidel_r(logicalKet, res_S4_00.data(), 5);
    fidelity_cum += fidelity_cum_tmp;
    trace_tmp = FastPauli::trace_r(res_S4_00.data(), 5);
    if (ecInfoIndex >= 0) {
        ecInfo.errors[4][ecInfoIndex] = 1 - fidelity_cum_tmp;
        ecInfo.traces[4][ecInfoIndex] = trace_tmp;
        ecInfo.grandErrors[ecInfoIndex] = 1 - fidelity_cum;
        // ecInfo.grandErrors.push_back(1 - FastPauli::fidel_r(logicalKet, cArr5_cum.data(), 5)); // Good test of cArr5_cum.
        ecInfo.grandTraces[ecInfoIndex] = ecInfo.traces[0][ecInfoIndex] + ecInfo.traces[1][ecInfoIndex] + ecInfo.traces[2][ecInfoIndex] + ecInfo.traces[3][ecInfoIndex] + ecInfo.traces[4][ecInfoIndex];
    }
    return cArr5_cum;
}

void FiveQubitCode::applyCorrection_flag(const uint16_t i, const uint16_t flag, double *cArr5) {
    constexpr int numQubits = 5;
    switch (flag) {
        case UINT16_C(3):
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
                case UINT16_C(2):  //  0  1  0  0 : X3Z4
                    FastPauli::X_r(cArr5, numQubits, 3);
                    FastPauli::Z_r(cArr5, numQubits, 4);
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
                case UINT16_C(13):  //  1  0  1  1 : Z3Z4
                    FastPauli::Z_r(cArr5, numQubits, 3);
                    FastPauli::Z_r(cArr5, numQubits, 4);
                    break;
                case UINT16_C(3):  //  1  1  0  0 : X4
                    FastPauli::X_r(cArr5, numQubits, 4);
                    break;
                case UINT16_C(11):  //  1  1  0  1 : Y3Z4
                    FastPauli::Y_r(cArr5, numQubits, 3);
                    FastPauli::Z_r(cArr5, numQubits, 4);
                    break;
                case UINT16_C(7):  //  1  1  1  0 : Z0Y1
                    FastPauli::Z_r(cArr5, numQubits, 0);
                    FastPauli::Y_r(cArr5, numQubits, 1);
                    break;
                case UINT16_C(15):  //  1  1  1  1 : Z0Z1
                    FastPauli::Z_r(cArr5, numQubits, 0);
                    FastPauli::Z_r(cArr5, numQubits, 1);
                    break;
                default:
                    break;
            }
            break;
        case UINT16_C(2):
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
                case UINT16_C(2):  //  0  1  0  0 : X3Z4
                    FastPauli::X_r(cArr5, numQubits, 3);
                    FastPauli::Z_r(cArr5, numQubits, 4);
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
                case UINT16_C(13):  //  1  0  1  1 : Z3Z4
                    FastPauli::Z_r(cArr5, numQubits, 3);
                    FastPauli::Z_r(cArr5, numQubits, 4);
                    break;
                case UINT16_C(3):  //  1  1  0  0 : X0Z2
                    FastPauli::X_r(cArr5, numQubits, 0);
                    FastPauli::Z_r(cArr5, numQubits, 2);
                    break;
                case UINT16_C(11):  //  1  1  0  1 : Y3Z4
                    FastPauli::Y_r(cArr5, numQubits, 3);
                    FastPauli::Z_r(cArr5, numQubits, 4);
                    break;
                case UINT16_C(7):  //  1  1  1  0 : Y4
                    FastPauli::Y_r(cArr5, numQubits, 4);
                    break;
                case UINT16_C(15):  //  1  1  1  1 : Z0Z1
                    FastPauli::Z_r(cArr5, numQubits, 0);
                    FastPauli::Z_r(cArr5, numQubits, 1);
                    break;
                default:
                    break;
            }
            break;
        case UINT16_C(1):
            switch (i) {           // a1 a2 a3 a4 : Pi where P ∈ {X,Y,Z} is a single-qubit Pauli operator and i is the index of the qubit to correct with i ∈ {0,1,2,3,4};
                case UINT16_C(8):  //  0  0  0  1 : X1
                    FastPauli::X_r(cArr5, numQubits, 1);
                    break;
                case UINT16_C(4):  //  0  0  1  0 : X1X2
                    FastPauli::X_r(cArr5, numQubits, 1);
                    FastPauli::X_r(cArr5, numQubits, 2);
                    break;
                case UINT16_C(12):  //  0  0  1  1 : Y3X4
                    FastPauli::Y_r(cArr5, numQubits, 3);
                    FastPauli::X_r(cArr5, numQubits, 4);
                    break;
                case UINT16_C(2):  //  0  1  0  0 : Z2
                    FastPauli::Z_r(cArr5, numQubits, 2);
                    break;
                case UINT16_C(10):  //  0  1  0  1 : Z3X4
                    FastPauli::Z_r(cArr5, numQubits, 3);
                    FastPauli::X_r(cArr5, numQubits, 4);
                    break;
                case UINT16_C(6):  //  0  1  1  0 : X0Y4
                    FastPauli::X_r(cArr5, numQubits, 0);
                    FastPauli::Y_r(cArr5, numQubits, 4);
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
                case UINT16_C(5):  //  1  0  1  0 : X3X4
                    FastPauli::X_r(cArr5, numQubits, 3);
                    FastPauli::X_r(cArr5, numQubits, 4);
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
            break;
        case UINT16_C(0):
            switch (i) {           // a1 a2 a3 a4 : Pi where P ∈ {X,Y,Z} is a single-qubit Pauli operator and i is the index of the qubit to correct with i ∈ {0,1,2,3,4};
                case UINT16_C(8):  //  0  0  0  1 : Y2X3
                    FastPauli::Y_r(cArr5, numQubits, 2);
                    FastPauli::X_r(cArr5, numQubits, 3);
                    break;
                case UINT16_C(4):  //  0  0  1  0 : Z2X3
                    FastPauli::Z_r(cArr5, numQubits, 2);
                    FastPauli::X_r(cArr5, numQubits, 3);
                    break;
                case UINT16_C(12):  //  0  0  1  1 : X0Y1
                    FastPauli::X_r(cArr5, numQubits, 0);
                    FastPauli::Y_r(cArr5, numQubits, 1);
                    break;
                case UINT16_C(2):  //  0  1  0  0 : Z2
                    FastPauli::Z_r(cArr5, numQubits, 2);
                    break;
                case UINT16_C(10):  //  0  1  0  1 : X2X3
                    FastPauli::X_r(cArr5, numQubits, 2);
                    FastPauli::X_r(cArr5, numQubits, 3);
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
                case UINT16_C(9):  //  1  0  0  1 : X0X1
                    FastPauli::X_r(cArr5, numQubits, 0);
                    FastPauli::X_r(cArr5, numQubits, 1);
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
            break;
        default:
            break;
    }
}

// void FiveQubitCode::applyCorrection_flag(const uint16_t i, const uint16_t flag, double *cArr5) {
//     constexpr int numQubits = 5;
//     switch (flag) {
//         case UINT16_C(0):
//             switch (i) {           // a1 a2 a3 a4 : Pi where P ∈ {X,Y,Z} is a single-qubit Pauli operator and i is the index of the qubit to correct with i ∈ {0,1,2,3,4};
//                 case UINT16_C(1):  //  0  0  0  1 : X1
//                     FastPauli::X_r(cArr5, numQubits, 1);
//                     break;
//                 case UINT16_C(2):  //  0  0  1  0 : Z4
//                     FastPauli::Z_r(cArr5, numQubits, 4);
//                     break;
//                 case UINT16_C(3):  //  0  0  1  1 : X2
//                     FastPauli::X_r(cArr5, numQubits, 2);
//                     break;
//                 case UINT16_C(4):  //  0  1  0  0 : X3Z4
//                     FastPauli::X_r(cArr5, numQubits, 3);
//                     FastPauli::Z_r(cArr5, numQubits, 4);
//                     break;
//                 case UINT16_C(5):  //  0  1  0  1 : Z0
//                     FastPauli::Z_r(cArr5, numQubits, 0);
//                     break;
//                 case UINT16_C(6):  //  0  1  1  0 : X3
//                     FastPauli::X_r(cArr5, numQubits, 3);
//                     break;
//                 case UINT16_C(7):  //  0  1  1  1 : Y2
//                     FastPauli::Y_r(cArr5, numQubits, 2);
//                     break;
//                 case UINT16_C(8):  //  1  0  0  0 : X0
//                     FastPauli::X_r(cArr5, numQubits, 0);
//                     break;
//                 case UINT16_C(9):  //  1  0  0  1 : Z3
//                     FastPauli::Z_r(cArr5, numQubits, 3);
//                     break;
//                 case UINT16_C(10):  //  1  0  1  0 : Z1
//                     FastPauli::Z_r(cArr5, numQubits, 1);
//                     break;
//                 case UINT16_C(11):  //  1  0  1  1 : Z3Z4
//                     FastPauli::Z_r(cArr5, numQubits, 3);
//                     FastPauli::Z_r(cArr5, numQubits, 4);
//                     break;
//                 case UINT16_C(12):  //  1  1  0  0 : X4
//                     FastPauli::X_r(cArr5, numQubits, 4);
//                     break;
//                 case UINT16_C(13):  //  1  1  0  1 : Y3Z4
//                     FastPauli::Y_r(cArr5, numQubits, 3);
//                     FastPauli::Z_r(cArr5, numQubits, 4);
//                     break;
//                 case UINT16_C(14):  //  1  1  1  0 : Z0Y1
//                     FastPauli::Z_r(cArr5, numQubits, 0);
//                     FastPauli::Y_r(cArr5, numQubits, 1);
//                     break;
//                 case UINT16_C(15):  //  1  1  1  1 : Z0Z1
//                     FastPauli::Z_r(cArr5, numQubits, 0);
//                     FastPauli::Z_r(cArr5, numQubits, 1);
//                     break;
//                 default:
//                     break;
//             }
//             break;
//         case UINT16_C(1):
//             switch (i) {           // a1 a2 a3 a4 : Pi where P ∈ {X,Y,Z} is a single-qubit Pauli operator and i is the index of the qubit to correct with i ∈ {0,1,2,3,4};
//                 case UINT16_C(1):  //  0  0  0  1 : X1
//                     FastPauli::X_r(cArr5, numQubits, 1);
//                     break;
//                 case UINT16_C(2):  //  0  0  1  0 : Z4
//                     FastPauli::Z_r(cArr5, numQubits, 4);
//                     break;
//                 case UINT16_C(3):  //  0  0  1  1 : X2
//                     FastPauli::X_r(cArr5, numQubits, 2);
//                     break;
//                 case UINT16_C(4):  //  0  1  0  0 : X3Z4
//                     FastPauli::X_r(cArr5, numQubits, 3);
//                     FastPauli::Z_r(cArr5, numQubits, 4);
//                     break;
//                 case UINT16_C(5):  //  0  1  0  1 : Z0
//                     FastPauli::Z_r(cArr5, numQubits, 0);
//                     break;
//                 case UINT16_C(6):  //  0  1  1  0 : X3
//                     FastPauli::X_r(cArr5, numQubits, 3);
//                     break;
//                 case UINT16_C(7):  //  0  1  1  1 : Y2
//                     FastPauli::Y_r(cArr5, numQubits, 2);
//                     break;
//                 case UINT16_C(8):  //  1  0  0  0 : X0
//                     FastPauli::X_r(cArr5, numQubits, 0);
//                     break;
//                 case UINT16_C(9):  //  1  0  0  1 : Z3
//                     FastPauli::Z_r(cArr5, numQubits, 3);
//                     break;
//                 case UINT16_C(10):  //  1  0  1  0 : Z1
//                     FastPauli::Z_r(cArr5, numQubits, 1);
//                     break;
//                 case UINT16_C(11):  //  1  0  1  1 : Z3Z4
//                     FastPauli::Z_r(cArr5, numQubits, 3);
//                     FastPauli::Z_r(cArr5, numQubits, 4);
//                     break;
//                 case UINT16_C(12):  //  1  1  0  0 : X0Z2
//                     FastPauli::X_r(cArr5, numQubits, 0);
//                     FastPauli::Z_r(cArr5, numQubits, 2);
//                     break;
//                 case UINT16_C(13):  //  1  1  0  1 : Y3Z4
//                     FastPauli::Y_r(cArr5, numQubits, 3);
//                     FastPauli::Z_r(cArr5, numQubits, 4);
//                     break;
//                 case UINT16_C(14):  //  1  1  1  0 : Y4
//                     FastPauli::Y_r(cArr5, numQubits, 4);
//                     break;
//                 case UINT16_C(15):  //  1  1  1  1 : Z0Z1
//                     FastPauli::Z_r(cArr5, numQubits, 0);
//                     FastPauli::Z_r(cArr5, numQubits, 1);
//                     break;
//                 default:
//                     break;
//             }
//             break;
//         case UINT16_C(2):
//             switch (i) {           // a1 a2 a3 a4 : Pi where P ∈ {X,Y,Z} is a single-qubit Pauli operator and i is the index of the qubit to correct with i ∈ {0,1,2,3,4};
//                 case UINT16_C(1):  //  0  0  0  1 : X1
//                     FastPauli::X_r(cArr5, numQubits, 1);
//                     break;
//                 case UINT16_C(2):  //  0  0  1  0 : X1X2
//                     FastPauli::X_r(cArr5, numQubits, 1);
//                     FastPauli::X_r(cArr5, numQubits, 2);
//                     break;
//                 case UINT16_C(3):  //  0  0  1  1 : Y3X4
//                     FastPauli::Y_r(cArr5, numQubits, 3);
//                     FastPauli::X_r(cArr5, numQubits, 4);
//                     break;
//                 case UINT16_C(4):  //  0  1  0  0 : Z2
//                     FastPauli::Z_r(cArr5, numQubits, 2);
//                     break;
//                 case UINT16_C(5):  //  0  1  0  1 : Z3X4
//                     FastPauli::Z_r(cArr5, numQubits, 3);
//                     FastPauli::X_r(cArr5, numQubits, 4);
//                     break;
//                 case UINT16_C(6):  //  0  1  1  0 : X0Y4
//                     FastPauli::X_r(cArr5, numQubits, 0);
//                     FastPauli::Y_r(cArr5, numQubits, 4);
//                     break;
//                 case UINT16_C(7):  //  0  1  1  1 : Y2
//                     FastPauli::Y_r(cArr5, numQubits, 2);
//                     break;
//                 case UINT16_C(8):  //  1  0  0  0 : X0
//                     FastPauli::X_r(cArr5, numQubits, 0);
//                     break;
//                 case UINT16_C(9):  //  1  0  0  1 : Z3
//                     FastPauli::Z_r(cArr5, numQubits, 3);
//                     break;
//                 case UINT16_C(10):  //  1  0  1  0 : X3X4
//                     FastPauli::X_r(cArr5, numQubits, 3);
//                     FastPauli::X_r(cArr5, numQubits, 4);
//                     break;
//                 case UINT16_C(11):  //  1  0  1  1 : Y1
//                     FastPauli::Y_r(cArr5, numQubits, 1);
//                     break;
//                 case UINT16_C(12):  //  1  1  0  0 : X4
//                     FastPauli::X_r(cArr5, numQubits, 4);
//                     break;
//                 case UINT16_C(13):  //  1  1  0  1 : Y0
//                     FastPauli::Y_r(cArr5, numQubits, 0);
//                     break;
//                 case UINT16_C(14):  //  1  1  1  0 : Y4
//                     FastPauli::Y_r(cArr5, numQubits, 4);
//                     break;
//                 case UINT16_C(15):  //  1  1  1  1 : Y3
//                     FastPauli::Y_r(cArr5, numQubits, 3);
//                     break;
//                 default:
//                     break;
//             }
//             break;
//         case UINT16_C(3):
//             switch (i) {           // a1 a2 a3 a4 : Pi where P ∈ {X,Y,Z} is a single-qubit Pauli operator and i is the index of the qubit to correct with i ∈ {0,1,2,3,4};
//                 case UINT16_C(1):  //  0  0  0  1 : Y2X3
//                     FastPauli::Y_r(cArr5, numQubits, 2);
//                     FastPauli::X_r(cArr5, numQubits, 3);
//                     break;
//                 case UINT16_C(2):  //  0  0  1  0 : Z2X3
//                     FastPauli::Z_r(cArr5, numQubits, 2);
//                     FastPauli::X_r(cArr5, numQubits, 3);
//                     break;
//                 case UINT16_C(3):  //  0  0  1  1 : X0Y1
//                     FastPauli::X_r(cArr5, numQubits, 0);
//                     FastPauli::Y_r(cArr5, numQubits, 1);
//                     break;
//                 case UINT16_C(4):  //  0  1  0  0 : Z2
//                     FastPauli::Z_r(cArr5, numQubits, 2);
//                     break;
//                 case UINT16_C(5):  //  0  1  0  1 : X2X3
//                     FastPauli::X_r(cArr5, numQubits, 2);
//                     FastPauli::X_r(cArr5, numQubits, 3);
//                     break;
//                 case UINT16_C(6):  //  0  1  1  0 : X3
//                     FastPauli::X_r(cArr5, numQubits, 3);
//                     break;
//                 case UINT16_C(7):  //  0  1  1  1 : Y2
//                     FastPauli::Y_r(cArr5, numQubits, 2);
//                     break;
//                 case UINT16_C(8):  //  1  0  0  0 : X0
//                     FastPauli::X_r(cArr5, numQubits, 0);
//                     break;
//                 case UINT16_C(9):  //  1  0  0  1 : X0X1
//                     FastPauli::X_r(cArr5, numQubits, 0);
//                     FastPauli::X_r(cArr5, numQubits, 1);
//                     break;
//                 case UINT16_C(10):  //  1  0  1  0 : Z1
//                     FastPauli::Z_r(cArr5, numQubits, 1);
//                     break;
//                 case UINT16_C(11):  //  1  0  1  1 : Y1
//                     FastPauli::Y_r(cArr5, numQubits, 1);
//                     break;
//                 case UINT16_C(12):  //  1  1  0  0 : X4
//                     FastPauli::X_r(cArr5, numQubits, 4);
//                     break;
//                 case UINT16_C(13):  //  1  1  0  1 : Y0
//                     FastPauli::Y_r(cArr5, numQubits, 0);
//                     break;
//                 case UINT16_C(14):  //  1  1  1  0 : Y4
//                     FastPauli::Y_r(cArr5, numQubits, 4);
//                     break;
//                 case UINT16_C(15):  //  1  1  1  1 : Y3
//                     FastPauli::Y_r(cArr5, numQubits, 3);
//                     break;
//                 default:
//                     break;
//             }
//             break;
//         default:
//             break;
//     }
// }

FastPauli::CArr<5> FiveQubitCode::applyCorrectionArrays_flag(const double *cArr9, std::array<std::array<double, 16>, 16> const &ancillaKets, const uint16_t flag) {
    constexpr int numAncillaQubits = 4;
    constexpr int numQubits = 5;
    constexpr uint64_t cArrLen = FastPauli::CArrLen<numQubits>;
    std::array<double, cArrLen> dm;
    dm.fill(0.0);
    std::array<double, cArrLen> tempDm;

    for (uint16_t a = 0; a < 16; a++) {
        FiveQubitCode::applyAncillaKet(tempDm.data(), numAncillaQubits, numQubits, ancillaKets[a].data(), cArr9);
        // Apply Pauli correction to tempDm.
        FiveQubitCode::applyCorrection_flag(a, flag, tempDm.data());

        // Store corrected tempDm in dm.
        for (uint64_t k = 0; k < cArrLen; k++) {
            dm[k] += tempDm[k];
        }
    }

    // Return the corrected 5-qubit density matrix (C-array form).
    return dm;
}