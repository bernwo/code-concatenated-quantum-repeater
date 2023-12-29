#include <boost/program_options.hpp>
#include <chrono>  // https://stackoverflow.com/questions/22387586/measuring-execution-time-of-a-function-in-c
#include <cmath>
#include <cstdint>
#include <cstdlib>  // Provides EXIT_SUCCESS
#include <filesystem>
#ifdef __APPLE__
    #ifndef FMT_HEADER_ONLY
        #define FMT_HEADER_ONLY
    #endif
    #include <fmt/core.h>  // Use `homebrew` and run `brew install fmt` on your terminal to install fmt. We have to use `fmt` instead because `std::format` is not supported by Apple's Clang.
using fmt::format;
#else
    #include <format>
using std::format;
#endif
#include <iomanip>  // https://stackoverflow.com/questions/554063/how-do-i-print-a-double-value-with-full-precision-using-cout
#include <iostream>
#include <limits>
#include <tuple>
#include <vector>

#include "AuxMath.hpp"
#include "Repeater.hpp"
#include "Tree.hpp"

int main(int argc, const char *argv[]) {
    std::cout << "    \x1b[3;97mWritten by Wo Kah Jen.\x1b[0m" << std::endl;
#ifdef _WIN64
    std::cout << "    \x1b[3;97mProgram is running on a Windows x64 machine.\x1b[0m\n"
              << std::endl;
#elif __APPLE__
    std::cout << "    \x1b[3;97mProgram is running on an Apple machine.\x1b[0m\n"
              << std::endl;
#elif __linux__
    std::cout << "    \x1b[3;97mProgram is running on a Linux machine.\x1b[0m\n"
              << std::endl;
#else
    std::cout << "    \x1b[3;97mProgram is running on a machine that is neither Windows x64, MacOS, nor Linux.\x1b[0m\n"
              << std::endl;
#endif

    // Variables to declare immediately.
    const std::string data_path = "./binary_data/";
    const std::string filename = format("{}latest.dat", data_path);

    constexpr double Latt_val = 20e3;      // Attenuation length.
    constexpr double eta_d_val = 0.95;     // Overall efficiency.
    constexpr double tau_ph_val = 1e-9;    // Single photon emission time.
    constexpr double tau_ss_val = 100e-9;  // Spin-spin gate time.
    constexpr double L_min = 100e3;        // Minimum total distance in units of metre.
    constexpr double L_max = 10000e3;      // Maximum total distance in units of metre.
    constexpr double min_L0 = 1e3;         // The minimum inter-repeater distance to check.
    constexpr double max_L0 = 5e3;         // The maximum inter-repeater distance to check.
    std::vector<std::pair<double, int>>
        L0AndntotArray;  // The array of std::pair which stores the possible L0 values and ntot values, where ntot is the
                         // maximum total number of links in the repeater network.

    // Set the variables for reading in pre-computed error rates.
    constexpr int reencodingErrorFactor = 3;
    constexpr int teleportedErrorFactor = 3;
    constexpr int len_er = 5;
    constexpr std::array<double, len_er> er_vals{1e-4, 3e-4, 5e-4, 1e-3, 2e-3};

    // Variables to declare with `boost::program_options`.
    std::string SKRMode;  // Whether to run the SKR calculation for the hybrid repeater protocol or the homogeneous
                          // repeater protocol. Can be either "hybrid" or "homogeneous".
    bool delete_data;
    int L_points;      // The number of points between L_min and L_max, inclusive.
    int numTrees;      // The number of trees to be included from the tree vectors in the calculation of the secret key rate.
    size_t err_index;  // The index of the std::array `er_vals` to perform secret key rate calculation on.
    double kappa;
    boost::program_options::options_description desc("\x1b[4;97mCommand line options\x1b[0m");

    desc.add_options()("help,h", "Produce help message\n")(
        "L_points,l", boost::program_options::value<int>(&L_points)->default_value(21),
        format("Number of data points between {:.2e} and {:.2e} kilometres, inclusive\n", L_min / 1e3, L_max / 1e3)
            .c_str())(
        "numTrees,n", boost::program_options::value<int>(&numTrees)->default_value(10),
        "The number of trees to be included from the tree vectors in the calculation of the secret key rate\n")(
        "kappa,k", boost::program_options::value<double>(&kappa)->default_value(1.0),
        "The relative cost between type i and type ii nodes\n")(
        "err_index,e", boost::program_options::value<size_t>(&err_index)->default_value(0),
        format("The index of the `er_vals = [{:.0e},{:.0e},{:.0e},{:.0e},{:.0e}]` to perform secret key rate "
               "calculation on\n",
               er_vals[0], er_vals[1], er_vals[2], er_vals[3], er_vals[4])
            .c_str())("delete_data,d", boost::program_options::bool_switch(&delete_data)->default_value(false),
                      format("Delete the file at {} if it exists.\n", filename).c_str())(
        "mode,m", boost::program_options::value<std::string>(&SKRMode)->default_value("hybrid"),
        "Whether to run the SKR calculation for the hybrid repeater protocol or the homogeneous repeater protocol. Can "
        "be either \"hybrid\" or \"homogeneous\".\n")(
        "verbose,v", boost::program_options::bool_switch()->default_value(false),
        "Turn verbose mode on\n")("ignoreErasures,i", boost::program_options::bool_switch()->default_value(false),
                                  "If true, then ignores erasure errors on the [[5,1,3]] code level\n");
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
    boost::program_options::notify(vm);

    if (vm.count("help") || argc == 1) {
        std::cout << desc << std::endl;
        return EXIT_SUCCESS;
    }
    const bool verbose = vm["verbose"].as<bool>();  // Program is more verbose if this is set to true.
    const bool ignoreErasures = vm["ignoreErasures"].as<bool>();
    if (vm.count("L_points")) {
        L_points = vm["L_points"].as<int>();
    }
    if (vm.count("numTrees")) {
        numTrees = vm["numTrees"].as<int>();
    }
    kappa = vm["kappa"].as<double>();
    if (vm.count("err_index")) {
        err_index = vm["err_index"].as<size_t>();
        if (err_index >= er_vals.size()) {
            std::cerr << "\x1b[0;91mError: The input err_index is too high!\x1b[0m\n";
            return EXIT_FAILURE;
        } else {
            std::cout << format("Reencoding error chosen: {:.0e}", er_vals[err_index]) << std::endl;
        }
    }

    std::cout << "  kappa = " << kappa << std::endl;

    if (vm["delete_data"].as<bool>()) {
        // Reference: https://stackoverflow.com/questions/13799670/c-cin-only-boolean0-1
        std::cout << format(
                         "\x1b[1;93mWarning\x1b[0m\x1b[0;93m: \x1b[0m\x1b[3;93mAre you sure you want to delete "
                         "\x1b[0m\x1b[4;93m{}\x1b[0m\x1b[0;93m? Please enter 0 or 1 only. 0 for keeping the "
                         "file and 1 for deleting the file. Your input:\x1b[0m",
                         filename)
                  << std::endl;
        while (!(std::cin >> delete_data)) {
            std::cout << format(
                             "\x1b[1;93mWarning\x1b[0m\x1b[0;93m: \x1b[0m\x1b[3;93mAre you sure you want to "
                             "delete \x1b[0m\x1b[4;93m{}\x1b[0m\x1b[0;93m? Please enter 0 or 1 only. 0 for "
                             "keeping the file and 1 for deleting the file. Your input:\x1b[0m",
                             filename)
                      << std::endl;
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
        if (delete_data) {
            // Reference: https://en.cppreference.com/w/cpp/filesystem/remove
            if (std::filesystem::remove(filename)) {
                std::cout << "  File deleted!" << std::endl;
            } else {
                std::cout << "  File doesn't exist to begin with so nothing to delete!" << std::endl;
            }
        } else {
            std::cout << "  File kept!" << std::endl;
        }
    }

    const auto t1 =
        std::chrono::high_resolution_clock::now();  // Start keeping track of elapsed time of running this program.

    // Initialise variables dependent on the `boost::program_options` input values.
    const std::vector<double> L_vals = AuxMath::logspace(
        std::log10(L_min), std::log10(L_max), L_points);  // The array of total distances of length `L_points`.

    // Read the pre-computed error rates.
    BinaryData::ECFlagInfo *ecFlagInfos;
    BinaryData::ECErasureInfo *ecErasureInfos;

    if (SKRMode == "hybrid") {
        ecFlagInfos = new BinaryData::ECFlagInfo[len_er];
        ecErasureInfos = new BinaryData::ECErasureInfo[len_er];
        if (std::filesystem::exists("../effective_error_probability/binary_data/") && err_index < er_vals.size()) {
            if (verbose) {
                std::cout << "\x1b[0;92mPath to the binary data exists!\x1b[0m" << std::endl;
            }
            ecFlagInfos[0].read(
                format("../effective_error_probability/binary_data/flag_err{:.0e}_reencodingFactor{}_teleFactor{}/"
                       "ecInfo0_2023-03-05 20-26-59.dat",
                       er_vals[0], reencodingErrorFactor, teleportedErrorFactor),
                format("../effective_error_probability/binary_data/flag_err{:.0e}_reencodingFactor{}_teleFactor{}/"
                       "ecInfo1_2023-03-05 20-26-59.dat",
                       er_vals[0], reencodingErrorFactor, teleportedErrorFactor));
            ecFlagInfos[1].read(
                format("../effective_error_probability/binary_data/flag_err{:.0e}_reencodingFactor{}_teleFactor{}/"
                       "ecInfo0_2023-03-05 20-17-13.dat",
                       er_vals[1], reencodingErrorFactor, teleportedErrorFactor),
                format("../effective_error_probability/binary_data/flag_err{:.0e}_reencodingFactor{}_teleFactor{}/"
                       "ecInfo1_2023-03-05 20-17-13.dat",
                       er_vals[1], reencodingErrorFactor, teleportedErrorFactor));
            ecFlagInfos[2].read(
                format("../effective_error_probability/binary_data/flag_err{:.0e}_reencodingFactor{}_teleFactor{}/"
                       "ecInfo0_2023-03-05 20-06-53.dat",
                       er_vals[2], reencodingErrorFactor, teleportedErrorFactor),
                format("../effective_error_probability/binary_data/flag_err{:.0e}_reencodingFactor{}_teleFactor{}/"
                       "ecInfo1_2023-03-05 20-06-53.dat",
                       er_vals[2], reencodingErrorFactor, teleportedErrorFactor));
            ecFlagInfos[3].read(
                format("../effective_error_probability/binary_data/flag_err{:.0e}_reencodingFactor{}_teleFactor{}/"
                       "ecInfo0_2023-03-05 19-57-33.dat",
                       er_vals[3], reencodingErrorFactor, teleportedErrorFactor),
                format("../effective_error_probability/binary_data/flag_err{:.0e}_reencodingFactor{}_teleFactor{}/"
                       "ecInfo1_2023-03-05 19-57-33.dat",
                       er_vals[3], reencodingErrorFactor, teleportedErrorFactor));
            ecFlagInfos[4].read(
                format("../effective_error_probability/binary_data/flag_err{:.0e}_reencodingFactor{}_teleFactor{}/"
                       "ecInfo0_2023-03-05 19-33-57.dat",
                       er_vals[4], reencodingErrorFactor, teleportedErrorFactor),
                format("../effective_error_probability/binary_data/flag_err{:.0e}_reencodingFactor{}_teleFactor{}/"
                       "ecInfo1_2023-03-05 19-33-57.dat",
                       er_vals[4], reencodingErrorFactor, teleportedErrorFactor));

            ecErasureInfos[0].read(
                format("../effective_error_probability/binary_data/1erasure_err{:.0e}_reencodingFactor{}_teleFactor{}/"
                       "ecInfo0_2023-03-05 20-33-04.dat",
                       er_vals[0], reencodingErrorFactor, teleportedErrorFactor));
            ecErasureInfos[1].read(
                format("../effective_error_probability/binary_data/1erasure_err{:.0e}_reencodingFactor{}_teleFactor{}/"
                       "ecInfo0_2023-03-05 20-32-09.dat",
                       er_vals[1], reencodingErrorFactor, teleportedErrorFactor));
            ecErasureInfos[2].read(
                format("../effective_error_probability/binary_data/1erasure_err{:.0e}_reencodingFactor{}_teleFactor{}/"
                       "ecInfo0_2023-03-05 20-30-52.dat",
                       er_vals[2], reencodingErrorFactor, teleportedErrorFactor));
            ecErasureInfos[3].read(
                format("../effective_error_probability/binary_data/1erasure_err{:.0e}_reencodingFactor{}_teleFactor{}/"
                       "ecInfo0_2023-03-05 20-29-28.dat",
                       er_vals[3], reencodingErrorFactor, teleportedErrorFactor));
            ecErasureInfos[4].read(
                format("../effective_error_probability/binary_data/1erasure_err{:.0e}_reencodingFactor{}_teleFactor{}/"
                       "ecInfo0_2023-03-05 20-28-37.dat",
                       er_vals[4], reencodingErrorFactor, teleportedErrorFactor));

            std::cout << "\x1b[0;92mBinary data has been read successfully\x1b[0m" << std::endl;
        } else {
            std::cerr << "\x1b[0;91mPath to the binary data does not exist!\x1b[0m\n";
            return EXIT_FAILURE;
        }
    } else {
        std::cout << "\x1b[0;93mSkipping reading pre-computed error rates\x1b[0m" << std::endl;
    }

    // Compute the all the tree vectors that adhere to the max number of photons constraint.
    constexpr uint16_t NMax = 300;
    if (verbose) {
        std::cout << "Max number of photons in the tree vector: " << NMax << std::endl;
    }
    Tree tree = Tree(NMax, 2, 2, 2, verbose);
    std::vector<std::vector<uint16_t>> treevec = tree.getVec();
    std::vector<std::vector<uint16_t>> treevec_sorted;

    // Print the number of photons, transmission rate, and generation time for each tree.
    if (verbose) {
        std::cout << std::setw(4) << "\x1b[1;4;93mNum photons,\x1b[0m " << std::setw(3)
                  << "\x1b[1;4;94mTree vectors,\x1b[0m " << std::setw(std::numeric_limits<double>::digits10 - 2) << " "
                  << "\x1b[1;4;95meta_e,\x1b[0m "
                  << "  "
                  << "\x1b[1;4;96mtau_tree,\x1b[0m" << std::endl;
        for (size_t i = 0; i < treevec.size(); i++) {
            std::cout << "\x1b[0;93m" << std::setw(11) << Tree::getPhotonNumber(treevec[i]) << ", \x1b[0;94m"
                      << std::setw(4) << treevec[i][0] << ", " << std::setw(2) << treevec[i][1] << ", " << std::setw(2)
                      << treevec[i][2] << ", "
                      << "\x1b[0;95m" << std::setprecision(std::numeric_limits<double>::digits10)
                      << std::setw(std::numeric_limits<double>::digits10 + 3)
                      << Repeater::eta_e(min_L0, Latt_val, eta_d_val, treevec[i]) << ", "
                      << "\x1b[0;96m" << std::setw(10) << std::setprecision(4)
                      << Repeater::tauTree(1e-9, 100e-9, treevec[i]) << ", "
                      << "\x1b[0m" << std::endl;
        }
        for (auto const &el : L0AndntotArray) {
            std::cout << std::setw(10) << std::setprecision(7) << el.first << ", " << std::setw(10) << el.second
                      << std::endl;
        }
    }

    std::vector<double> tauTreeArray;
    std::tuple<double, double, double, double, double, double> pair_tmp;
    std::tuple<double, double, double, double> tuple_homogeneous_tmp;
    double cost_tmp;
    double SKR_tmp;
    double max_SKR = std::numeric_limits<double>::max();
    double min_cost = std::numeric_limits<double>::max();
    int max_c;
    size_t argmax_L0andntot;
    size_t argmax_tree;
    std::vector<uint16_t> optim_treeVec;
    double optim_pTransSum = std::numeric_limits<double>::quiet_NaN();
    double optim_logicalError = std::numeric_limits<double>::quiet_NaN();
    double optim_eta_e = std::numeric_limits<double>::quiet_NaN();
    int optim_ntot = std::numeric_limits<int>::quiet_NaN();
    double optim_f = std::numeric_limits<double>::quiet_NaN();
    bool toPerformCalculation = true;

    // Begin main loop.
    Repeater::RepeaterInfo repeaterInfo;
    repeaterInfo.resize(L_vals.size(), 3);              // The std::vector's inside repeaterInfo have length of L_vals.size() and the
                                                        // tree vectors have length of 3, i.e., \Vec{t}={b_0,b_1,b_2}.
    repeaterInfo.reencodingError = er_vals[err_index];  // Denoted as ϵ_r. It is the reencoding error rate.
    repeaterInfo.reencodingErrorFactor =
        reencodingErrorFactor;  // reencodingErrorFactor is such that ϵ_r = reencodingErrorFactor × ϵ_0, where ϵ_0 is the
                                // depolarising error of each single qubit in the generated tree. ϵ_0 is also the error
                                // rate of (non-teleported) spin-spin gates.
    repeaterInfo.teleportedErrorFactor =
        teleportedErrorFactor;             // teleportedErrorFactor is such that ϵ_tele = teleportedErrorFactor × ϵ_0, where ϵ_tele
                                           // is the error rate of the teleported spin-spin gates.
    repeaterInfo.Latt_val = Latt_val;      // Attenuation length.
    repeaterInfo.eta_d_val = eta_d_val;    // Overall efficiency.
    repeaterInfo.tau_ph_val = tau_ph_val;  // Single photon emission time.
    repeaterInfo.tau_ss_val = tau_ss_val;  // Spin-spin gate time.

    if (repeaterInfo.doSKRCalculation(filename)) {
        std::cout << "     \x1b[1;92m>\x1b[0m \x1b[1;4;3;94mBegin main loop.\x1b[0m" << std::endl;
        for (size_t i = 0; i < L_vals.size(); i++) {
            L0AndntotArray = Repeater::getL0AndntotArray(L_vals[i], min_L0, max_L0, 200,
                                                         false);  // At maximum 200 points between min_L0 and max_L0.
            for (size_t arg_L0andntot = 0; arg_L0andntot < L0AndntotArray.size(); arg_L0andntot++) {
                // Get sorted tree vectors according to the transmission rate η_e and truncate this list of trees to
                // `numTrees` trees only.
                treevec_sorted = Tree::getSortedTreeVecs(L0AndntotArray[arg_L0andntot].first, Latt_val, eta_d_val,
                                                         treevec, numTrees);

                // Get the std::vector of tree generation time for each trees in `treevec_sorted`.
                tauTreeArray = Repeater::getTauTreeArray(tau_ph_val, tau_ss_val, treevec_sorted);
                if (SKRMode == "hybrid") {
                    // for (int c = L0AndntotArray[arg_L0andntot].second; c <= L0AndntotArray[arg_L0andntot].second;
                    // c++) {
                    for (int c = 1; c <= L0AndntotArray[arg_L0andntot].second / 2;
                         c++) {  // This is where you control the number of type II nodes in the network
                        // for (int c = 1; c <= L0AndntotArray[arg_L0andntot].second; c++) {
                        for (size_t arg_tree = 0; arg_tree < treevec_sorted.size(); arg_tree++) {
                            if (toPerformCalculation) {
                                pair_tmp = Repeater::SKR(
                                    ecFlagInfos[err_index], ecErasureInfos[err_index], L_vals[i], c,
                                    ignoreErasures
                                        ? 0
                                        : c,  // This is zero if want to exclude erasure error correction. Else put `c`.
                                    L0AndntotArray[arg_L0andntot].first, Latt_val, L0AndntotArray[arg_L0andntot].second,
                                    eta_d_val, tau_ph_val, tau_ss_val, 10 * tau_ss_val, 10 * tau_ss_val,
                                    tauTreeArray[arg_tree], treevec_sorted[arg_tree], kappa, false);
                                cost_tmp = std::get<0>(pair_tmp);
                                SKR_tmp = std::get<1>(pair_tmp);
                            }
                            if (toPerformCalculation && !std::isinf(cost_tmp) && !std::isnan(cost_tmp) &&
                                cost_tmp < min_cost) {
                                argmax_tree = arg_tree;
                                optim_treeVec = treevec_sorted[argmax_tree];
                                argmax_L0andntot = arg_L0andntot;
                                max_c = c;
                                max_SKR = SKR_tmp;
                                min_cost = cost_tmp;
                                optim_logicalError = std::get<2>(pair_tmp);
                                optim_pTransSum = std::get<3>(pair_tmp);
                                optim_eta_e = std::get<4>(pair_tmp);
                                optim_ntot = L0AndntotArray[argmax_L0andntot].second;
                                optim_f = std::get<5>(pair_tmp);
                            }
                        }
                    }
                } else {
                    for (size_t arg_tree = 0; arg_tree < treevec_sorted.size(); arg_tree++) {
                        if (toPerformCalculation) {
                            tuple_homogeneous_tmp = Repeater::SKRHomogeneous(
                                L_vals[i], L0AndntotArray[arg_L0andntot].first, Latt_val,
                                L0AndntotArray[arg_L0andntot].second, eta_d_val, tau_ph_val, tau_ss_val,
                                tauTreeArray[arg_tree], treevec_sorted[arg_tree], er_vals[err_index]);
                            cost_tmp = std::get<0>(tuple_homogeneous_tmp);
                            SKR_tmp = std::get<1>(tuple_homogeneous_tmp);
                        }
                        if (toPerformCalculation && !std::isinf(cost_tmp) && !std::isnan(cost_tmp) &&
                            cost_tmp < min_cost) {
                            argmax_tree = arg_tree;
                            optim_treeVec = treevec_sorted[argmax_tree];
                            argmax_L0andntot = arg_L0andntot;
                            max_SKR = SKR_tmp;
                            min_cost = cost_tmp;
                            optim_pTransSum = std::get<2>(tuple_homogeneous_tmp);
                            optim_eta_e = std::get<3>(tuple_homogeneous_tmp);
                            optim_ntot = L0AndntotArray[argmax_L0andntot].second;
                        }
                    }
                }
            }
            // Save data into struct.
            if (max_SKR >= 1) {
                toPerformCalculation = true;
                repeaterInfo.SKR[i] = max_SKR;
                repeaterInfo.cost[i] = min_cost;
                repeaterInfo.L0[i] = L0AndntotArray[argmax_L0andntot].first;
                repeaterInfo.Ltot[i] = L_vals[i];
                repeaterInfo.c[i] = max_c;
                repeaterInfo.treeVecs[i] = optim_treeVec;
                repeaterInfo.pTransSum[i] = optim_pTransSum;
                repeaterInfo.logicalError[i] = optim_logicalError;
                repeaterInfo.eta_e[i] = optim_eta_e;
                repeaterInfo.ntot[i] = optim_ntot;
                repeaterInfo.f_val[i] = optim_f;
                // Get ready for next iteration.
                SKR_tmp = 0;
                max_SKR = 0;
                cost_tmp = std::numeric_limits<double>::max();
                min_cost = std::numeric_limits<double>::max();
                optim_pTransSum = 0;
                optim_logicalError = 1;
                optim_eta_e = 0;
                optim_ntot = 0;
                optim_f = 0;
            } else {
                // Set values corresponding to ~0 secret key rate.
                toPerformCalculation = false;
                repeaterInfo.SKR[i] = 0;
                repeaterInfo.cost[i] = std::numeric_limits<double>::quiet_NaN();
                repeaterInfo.L0[i] = std::numeric_limits<double>::quiet_NaN();
                repeaterInfo.Ltot[i] = L_vals[i];
                repeaterInfo.c[i] = std::numeric_limits<int>::quiet_NaN();
                repeaterInfo.treeVecs[i].resize(3, std::numeric_limits<uint16_t>::quiet_NaN());
                repeaterInfo.pTransSum[i] = std::numeric_limits<double>::quiet_NaN();
                repeaterInfo.logicalError[i] = std::numeric_limits<double>::quiet_NaN();
                repeaterInfo.eta_e[i] = std::numeric_limits<double>::quiet_NaN();
                repeaterInfo.ntot[i] = std::numeric_limits<int>::quiet_NaN();
                repeaterInfo.f_val[i] = std::numeric_limits<double>::quiet_NaN();
            }
            std::cout
                << format(
                       "    \x1b[1;92m>\x1b[0m \x1b[1;3;97mSKR: {:10.3e}, \x1b[1;3;91mcost: {:10.3e}, \x1b[1;3;92mL0: "
                       "{:8.3f}, \x1b[1;3;93mntot: {:5d}, \x1b[1;3;94mc: {:5d}, \x1b[1;3;95mtree[{:2d}]: "
                       "[{:2d},{:2d},{:2d}], \x1b[1;3;96mfor L = {:9.3f} km  \x1b[0m\x1b[1;92m  ->  {:6.2f}%\x1b[0m",
                       repeaterInfo.SKR[i], repeaterInfo.cost[i], repeaterInfo.L0[i], repeaterInfo.ntot[i],
                       repeaterInfo.c[i], argmax_tree, repeaterInfo.treeVecs[i][0], repeaterInfo.treeVecs[i][1],
                       repeaterInfo.treeVecs[i][2], L_vals[i] / 1e3, static_cast<double>(i + 1) / L_vals.size() * 100)
                << std::endl;
        }
        // Write data to disk.
        std::filesystem::create_directory(data_path);
        repeaterInfo.save(filename);
    } else {
        std::cout
            << format(
                   "     \x1b[1;92m>\x1b[0m \x1b[1;4;3;94mMain loop not executed. Data already exist for {:.0e}\x1b[0m",
                   er_vals[err_index])
            << std::endl;
    }

    // Free memory
    if (SKRMode == "hybrid") {
        delete[] ecFlagInfos;
        delete[] ecErasureInfos;
    }

    return EXIT_SUCCESS;
}
