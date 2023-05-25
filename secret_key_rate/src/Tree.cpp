// Copyright © Kah Jen Wo.
// C++ class for generating all possible tree vectors with fixed tree depth d=2 and constrained number of photons ≤
// NMax.
#include "Tree.hpp"

#include "Repeater.hpp"

Tree::Tree(const uint16_t NMax, const uint16_t aMin, const uint16_t bMin, const uint16_t cMin, const int verbosity) {
    this->verbose = verbosity >= 1;
    this->NMax    = NMax;
    this->aMin    = aMin;
    this->bMin    = bMin;
    this->cMin    = cMin;
    this->aMax    = std::floor(NMax / (1 + (1 + this->cMin) * this->bMin));
    this->cMax    = std::floor((NMax - this->aMin * (1 + this->bMin)) / (this->aMin * this->bMin));
    Tree::calcVecSize_();
    Tree::calcVec_();
}

void Tree::calcVecSize_() {
    // Algorithm self-derived. Original.
    this->vecSize = 0;
    int tmp;
    for (uint16_t a = this->aMin; a <= this->aMax; a++) {
        for (uint16_t c = this->cMin; c <= this->cMax; c++) {
            tmp = std::floor((this->NMax - a) / (a * (1 + c)));
            // std::cout << "tmp: " << tmp << "\n";
            this->vecSize += std::max(tmp - this->bMin + 1, 0);
            // if (tmp > this->bMin){
            //     this->vecSize += tmp-1;
            // }
        }
    }
}

uint32_t Tree::getSize() { return this->vecSize; }

void Tree::calcVec_() {
    const int d = 2; // Tree depth as defined in https://doi.org/10.1103/PhysRevX.10.021071.
    // uint16_t aMax;
    uint16_t a = this->aMin;
    uint16_t b = this->bMin;
    uint16_t c = this->cMin;
    // uint16_t aMaxTmp = std::floor(this->NMax/(1+(1+c)*b)); // Not used.
    uint16_t bMaxTmp = std::floor((this->NMax - a) / (a * (1 + c)));
    uint16_t cMaxTmp = std::floor((this->NMax - a * (1 + b)) / (a * b));

    this->treeVec.resize(this->vecSize, std::vector<uint16_t>(d + 1, 0));
    for (uint32_t i = 0; i < this->vecSize; i++) {
        // std::cout << "aMaxTmp: " << aMaxTmp << ", bMaxTmp: " << bMaxTmp << ", cMaxTmp: " << cMaxTmp << "\n";
        assert((static_cast<void>("Total number of photons exceeded NMax!"), a * (1 + b + b * c) <= this->NMax));
        // std::cout << a*(1+b+b*c) << "\n";
        treeVec[i][0] = a;
        treeVec[i][1] = b;
        treeVec[i][2] = c;
        if (c < cMaxTmp) {
            c++;
        } else {
            c = this->cMin;
            if (b < bMaxTmp) {
                b++;
            } else {
                b = this->bMin;
                if (a < this->aMax) {
                    a++;
                } else if (i != this->vecSize - 1) {
                    std::cout
                        << "\033[1;31mThis branch of the if-statement should never happen because the math is worked "
                           "meticulously.\nIf it happens, you probably need to find the bug and squash it.\033[0m\n";
                }
                bMaxTmp = ((this->NMax - a) / (a * (1 + c)));
            }
            cMaxTmp = ((this->NMax - a * (1 + b)) / (a * b));
        }
        // aMaxTmp = std::floor(this->NMax/(1+(1+c)*b));
    }
    if (a == this->aMax && b == this->bMin && c == this->cMin && this->verbose) {
        std::cout << "\033[1;32mTree vectors correctly generated. vecSize: " << this->vecSize << "\033[0m\n";
    }
}

std::vector<std::vector<uint16_t>> Tree::getVec() { return this->treeVec; }

void Tree::printVec() {
    std::cout << "\033[0;34m";
    for (auto const &i : this->treeVec) {
        for (auto const &j : i) {
            std::cout << std::setw(3) << j << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\033[0m\n";
}

uint16_t Tree::getPhotonNumber(std::vector<uint16_t> const &treeVec) {
    uint16_t tmp;
    uint16_t photonNumber = UINT16_C(0);
    for (size_t i = 0; i < treeVec.size(); i++) {
        tmp = UINT16_C(1);
        for (size_t j = 0; j <= i; j++) {
            tmp *= treeVec[j];
        }
        photonNumber += tmp;
    }
    return photonNumber;
}

std::vector<std::vector<uint16_t>> Tree::getSortedTreeVecs(const double L0, const double Latt, const double eta_d,
                                                           std::vector<std::vector<uint16_t>> const &treeVecs,
                                                           int length) {
    if (length < 0 || length > treeVecs.size()) {
        length = treeVecs.size();
    }

    std::vector<double> eta_eArray(treeVecs.size());
    for (size_t i = 0; i < treeVecs.size(); i++) {
        eta_eArray[i] = Repeater::eta_e(L0, Latt, eta_d, treeVecs[i]);
    }

    // Sort in descending order.
    std::vector<int> indices(treeVecs.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&](int A, int B) -> bool { return eta_eArray[A] > eta_eArray[B]; });

    std::vector<std::vector<uint16_t>> sortedTreeVecs;
    sortedTreeVecs.resize(length, std::vector<uint16_t>(treeVecs[0].size(), 0));
    for (int i = 0; i < length; i++) {
        for (size_t j = 0; j < treeVecs[0].size(); j++) {
            sortedTreeVecs[i][j] = treeVecs[indices[i]][j];
        }
    }
    return sortedTreeVecs;
}