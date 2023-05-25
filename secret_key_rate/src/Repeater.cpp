#include "Repeater.hpp"

std::vector<std::pair<double, int>> Repeater::getL0AndntotArray(const double Ltot, const double L0_min,
                                                                const double L0_max, const int numPerArray,
                                                                const bool verbose) {
    // constexpr double Ltot = 10000;  // In units of kilemetres.
    const int ntot_max     = static_cast<int>(std::floor(
        Ltot /
        1e3)); // Maximum total number of links in the repeater network chain. Computation done in units of kilometres.
    const int L0_max_floor = static_cast<int>(std::floor(L0_max / 1e3)); // In units of kilometres.
    int ntot_min           = 1; // Minimum total number of links in the repeater network chain.

    assert(L0_max_floor >= 1 && ntot_max > ntot_min);

    // Set the minimum total number of links in the repeater network chain.
    for (int ntot = 1; ntot < ntot_max; ntot++) {
        if (std::floor(Ltot / ntot / 1e3) == L0_max_floor) {
            ntot_min = ntot;
            if (verbose) {
                std::cout << "ntot_min: " << ntot << ", for Ltot: " << Ltot << "km and L0_max: " << L0_max << "km"
                          << std::endl;
            }
            break;
        }
    }
    assert(ntot_max > ntot_min);

    std::vector<double> L0arr_tmp; // In units of kilometres.
    std::vector<int> ntot_tmp;
    for (int n = ntot_min; n < ntot_max; n++) {
        L0arr_tmp.push_back(Ltot / n / 1e3);
        ntot_tmp.push_back(n);
    }
    std::reverse(L0arr_tmp.begin(), L0arr_tmp.end());
    std::reverse(ntot_tmp.begin(), ntot_tmp.end());

    int L0_min_floor = std::max(1, static_cast<int>(std::floor(L0_min / 1e3)));
    std::vector<std::vector<double>> L0arr(L0_max_floor - L0_min_floor + 1);
    std::vector<std::vector<double>> ntotarr(L0_max_floor - L0_min_floor + 1);

    for (size_t i = 0; i < L0arr_tmp.size(); i++) {
        for (size_t j = 1; j <= L0arr.size(); j++) {
            if (std::floor(L0arr_tmp[i]) == j + L0_min_floor - 1) {
                L0arr[j - 1].push_back(L0arr_tmp[i]);
                ntotarr[j - 1].push_back(ntot_tmp[i]);
            }
        }
    }

    std::vector<std::pair<double, int>> L0AndntotArray;
    for (size_t i = 0; i < L0arr.size(); i++) {
        for (size_t j = 0; j < L0arr[i].size();
             j += std::max(1, static_cast<int>(std::floor(static_cast<double>(L0arr[i].size()) / numPerArray)))) {
            L0AndntotArray.push_back(std::make_pair(L0arr[i][j] * 1e3, static_cast<int>(ntotarr[i][j])));
        }
    }

    return L0AndntotArray;
}

double Repeater::h(const double x) {
    assert(x >= 0 && x <= 1);
    if (x == 0 || x == 1) {
        return 0;
    } else {
        return -x * std::log2(x) - (1 - x) * std::log2(1 - x);
    }
}

double Repeater::f(const double error) {
    const double QBER = 2.0 / 3.0 * error;
    return std::max(1 - Repeater::h(QBER) - QBER - (1 - QBER) * Repeater::h((1 - 3.0 * QBER / 2.0) / (1 - QBER)), 0.0);
}

double Repeater::eta(const double L0, const double Latt) { return std::exp(-L0 / Latt); }

double Repeater::mu(const double L0, const double Latt, const double eta_d) {
    return 1 - Repeater::eta(L0, Latt) * eta_d;
}

double Repeater::eta_e(const double L0, const double Latt, const double eta_d, std::vector<uint16_t> const &treeVec) {
    assert(treeVec.size() == 3);
    const double mu_val = Repeater::mu(L0, Latt, eta_d);
    const double tmp0   = std::pow(1 - mu_val, 1 + treeVec[2]);
    const double tmp1   = mu_val * std::pow((1 - tmp0), treeVec[1]);
    return std::pow(1 - std::pow(mu_val, 1 + treeVec[2]), treeVec[1]) *
           (std::pow(1 - tmp1, treeVec[0]) - std::pow(mu_val - tmp1, treeVec[0]));
}

double Repeater::tauTree(const double tau_ph, const double tau_ss, std::vector<uint16_t> const &treeVec) {
    return treeVec[0] * (tau_ss * (3 + treeVec[1]) + tau_ph * (100 + treeVec[1] * (1 + treeVec[2])));
}

double Repeater::pTrans(const double eta_e_val, const double n, const int c, const int i) {
    assert(n >= 0 && c >= 0 && i >= 0 && c - i >= 0);
    return std::pow(eta_e_val, 5 * n * (c - i)) *
           std::pow(5 * std::pow(eta_e_val, 4 * n) * (1 - std::pow(eta_e_val, n)), i);
}

std::tuple<double, double, double, double, double, double>
Repeater::SKR(BinaryData::ECFlagInfo const &ecFlagInfo, BinaryData::ECErasureInfo const &ecErasureInfo,
              const double Ltot, const int c, const int cMax, const double L0, const double Latt, const int ntot,
              const double eta_d, const double tau_ph, const double tau_ss, const double tau_tele,
              const double tau_meas, const double tauTree, std::vector<uint16_t> const &treeVec, const double kappa,
              const bool verbose) {
    assert(c >= 0 && cMax >= 0 && cMax <= c);
    int which_division_case;
    double result    = 0;
    double pTransSum = 0;
    double cost;
    double logicalError = 0;
    double f_val        = 0;
    double pTrans;
    // double p;  // Ratio of the 1-erasure error occurrences to the total number of type II stations.

    const int n          = ntot / c;
    const double n_float = static_cast<double>(ntot) / c;                               // TESTING.
    const int np         = static_cast<int>(std::floor(static_cast<double>(ntot) / c)); // n'
    const int npp        = ntot - np * (c - 1);                                         // n''
    const int ntotModc   = ntot % c;

    const double tau_nf1 = tau_ss + 3 * tau_tele + tau_meas;
    const double tau_nf2 = 4 * tau_tele + tau_meas;
    const double tau_nf  = tau_nf1 + 3 * tau_nf2; // Total.
    const double tau_f1  = 2 * tau_ss + tau_nf1;
    const double tau_f2  = 2 * tau_tele + tau_nf2;
    // const double tau_loss = tau_nf + 3 * tau_meas;
    // const int worst_n = std::max(n, std::max(npp, np));
    // const double T_f = ecFlagInfo.traces[0][worst_n] * (tau_f + tau_nf) +
    //                    ecFlagInfo.traces[1][worst_n] * (2 * tau_f + tau_nf) +
    //                    ecFlagInfo.traces[2][worst_n] * (3 * tau_f + tau_nf) +
    //                    ecFlagInfo.traces[3][worst_n] * (4 * tau_f + tau_nf) +
    //                    ecFlagInfo.traces[4][worst_n] * (4 * tau_f);

    const double eta_e_val = Repeater::eta_e(L0, Latt, eta_d, treeVec);
    if (verbose) {
        std::cout << "(ntot / c) = (" << ntot << " / " << c << ") = n: " << n << ", n':" << np << ", n''" << npp
                  << std::endl;
    }

    const int iMax = std::min(c, cMax);
    for (int i = 0; i <= iMax; i++) {
        if (ntotModc == 0) {
            logicalError = 1 - std::pow(1 - ecErasureInfo.errs[n], i) * ecFlagInfo.getRecursiveFidelity(n, c);
            pTrans       = Repeater::pTrans(eta_e_val, n, c, i);
            if (verbose) {
                which_division_case = 0;
            }
        } else if ((c - 1 >= i) && (ntotModc > 0)) {
            logicalError = 1 - std::pow(1 - ecErasureInfo.errs[np], i) *
                                   ecFlagInfo.getRecursiveFidelity(np, c - 1 - i) *
                                   ecFlagInfo.getRecursiveFidelity(npp, 1);
            pTrans = std::pow(eta_e_val, 5 * npp) * std::pow(eta_e_val, 5 * np * (c - 1 - i)) *
                     std::pow(5 * std::pow(eta_e_val, 4 * np) * (1 - std::pow(eta_e_val, np)), i);
            if (verbose) {
                which_division_case = 1;
            }
        } else {
            logicalError = 1 - std::pow(1 - ecErasureInfo.errs[np], i - 1) * (1 - ecErasureInfo.errs[npp]);
            pTrans       = std::pow(5 * std::pow(eta_e_val, 4 * np) * (1 - std::pow(eta_e_val, np)), c - 1) * 5 *
                     std::pow(eta_e_val, 4 * npp) * (1 - std::pow(eta_e_val, npp));
            if (verbose) {
                which_division_case = 2;
            }
        }

        // p = (static_cast<double>(i) / c);
        f_val = Repeater::f(logicalError);
        // totProcessingTime = tauTree + (1 - p) * T_f + p * tau_loss;
        if ((!std::isnan(f_val)) && (f_val > 0)) {
            if (verbose) {
                switch (which_division_case) {
                case 0:
                    std::cout << "\x1b[0;91mLogical error: " << logicalError << ", for c = " << c << " and i = " << i
                              << "\x1b[0m" << std::endl;
                    std::cout << "More info:" << std::endl;
                    std::cout << "std::pow(1 - ecErasureInfo.errs[n], i): " << std::pow(1 - ecErasureInfo.errs[n], i)
                              << std::endl;
                    std::cout << "ecFlagInfo.getRecursiveFidelity(n, c): " << ecFlagInfo.getRecursiveFidelity(n, c)
                              << std::endl;
                    break;
                case 1:
                    std::cout << "\x1b[0;92mLogical error: " << logicalError << ", for c = " << c << " and i = " << i
                              << "\x1b[0m" << std::endl;
                    std::cout << "More info:" << std::endl;
                    std::cout << "std::pow(1 - ecErasureInfo.errs[np], i): " << std::pow(1 - ecErasureInfo.errs[np], i)
                              << std::endl;
                    std::cout << " ecFlagInfo.getRecursiveFidelity(np, c - 1 - i) : "
                              << ecFlagInfo.getRecursiveFidelity(np, c - 1 - i) << std::endl;
                    std::cout << " ecFlagInfo.getRecursiveFidelity(npp, 1) : "
                              << ecFlagInfo.getRecursiveFidelity(npp, 1) << std::endl;
                    break;
                case 2:
                    std::cout << "\x1b[0;93mLogical error: " << logicalError << ", for c = " << c << " and i = " << i
                              << "\x1b[0m" << std::endl;
                default:
                    break;
                }
            }
            pTransSum += pTrans;
            result += boost::math::binomial_coefficient<double>(c, i) * pTrans * f_val; // The SKR.
        } else {
            break;
        }
    }
    const double totProcessingTime = tauTree + 14 * tau_ss + 26 * tau_tele + 8 * tau_meas;
    result /= totProcessingTime;

    const double numberType1Nodes = (n_float - 1) * c;
    const double numberType2Nodes = c;

    if (result <= 0) {
        cost = INFINITY;
    } else {
        cost = Latt / (result * tau_ph * Ltot) * (numberType1Nodes + kappa * numberType2Nodes);
        // cost = 1.0 / result;
    }

    return std::make_tuple(cost, result, logicalError, pTransSum, eta_e_val, f_val);
}

std::tuple<double, double, double, double>
Repeater::SKRHomogeneous(const double Ltot, const double L0, const double Latt, const int ntot, const double eta_d,
                         const double tau_ph, const double tau_ss, const double tauTree,
                         std::vector<uint16_t> const &treeVec, const double errReencoding) {
    double result          = 0;        // The SKR.
    double cost            = INFINITY; // Cost.
    const double eta_e_val = Repeater::eta_e(L0, Latt, eta_d, treeVec);
    const double pTrans    = std::pow(eta_e_val, ntot);
    const double eTrans    = 1 - std::pow(1 - errReencoding, ntot);
    const double f_val     = Repeater::f(eTrans);
    if ((!std::isnan(f_val)) && (f_val > 0)) {
        result = pTrans * f_val / tauTree; // The SKR.
    }
    if (result > 0) {
        cost = ntot * Latt / (result * tau_ph * Ltot);
        // cost = 1.0 / result;
    }
    return std::make_tuple(cost, result, pTrans, eta_e_val);
}

std::vector<double> Repeater::getTauTreeArray(const double tau_ph, const double tau_ss,
                                              std::vector<std::vector<uint16_t>> const &treeVecs) {
    std::vector<double> tauTreeArray(treeVecs.size());
    for (size_t i = 0; i < treeVecs.size(); i++) {
        tauTreeArray[i] = Repeater::tauTree(tau_ph, tau_ss, treeVecs[i]);
    }
    return tauTreeArray;
}
