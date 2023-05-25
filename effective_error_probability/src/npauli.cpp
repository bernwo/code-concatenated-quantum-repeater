#include "npauli.hpp"

#include "fpauli.hpp"

void NoisyPauli::SingleQubitDepolarise_r(double *cArray, const int numQubits, const int targetQubit, const double err) {
    int i;
    const uint64_t cArrLen = FastPauli::getCArrLen(numQubits);
    double *x = new double[cArrLen];
    double *y = new double[cArrLen];
    double *z = new double[cArrLen];
    for (i = 0; i < cArrLen; i++) {
        x[i] = cArray[i];
        y[i] = cArray[i];
        z[i] = cArray[i];
    }
    FastPauli::X_r(x, numQubits, targetQubit);
    FastPauli::Y_r(y, numQubits, targetQubit);
    FastPauli::Z_r(z, numQubits, targetQubit);
    for (i = 0; i < cArrLen; i++) {
        cArray[i] = (1 - err) * cArray[i] + (err / 3.0) * (x[i] + y[i] + z[i]);
    }
    delete[] x;
    delete[] y;
    delete[] z;
}

void NoisyPauli::NCX_r(double *cArray, const int numQubits, const int controlQubit, const int targetQubit, const double err) {
    int i;
    const int cArrLen = (UINT64_C(1) << (numQubits - 1)) + (UINT64_C(1) << (2 * numQubits - 1));  // numQubits cannot be too high, otherwise overflow and cArrLen becomes negative!
    FastPauli::CX_r(cArray, numQubits, controlQubit, targetQubit);

    double *i1x2 = new double[cArrLen];
    double *i1y2 = new double[cArrLen];
    double *i1z2 = new double[cArrLen];
    double *x1x2 = new double[cArrLen];
    double *x1y2 = new double[cArrLen];
    double *x1z2 = new double[cArrLen];
    double *x1i2 = new double[cArrLen];
    double *y1x2 = new double[cArrLen];
    double *y1y2 = new double[cArrLen];
    double *y1z2 = new double[cArrLen];
    double *y1i2 = new double[cArrLen];
    double *z1x2 = new double[cArrLen];
    double *z1y2 = new double[cArrLen];
    double *z1z2 = new double[cArrLen];
    double *z1i2 = new double[cArrLen];

    for (i = 0; i < cArrLen; i++) {
        i1x2[i] = cArray[i];
        i1y2[i] = cArray[i];
        i1z2[i] = cArray[i];
        x1i2[i] = cArray[i];
        y1i2[i] = cArray[i];
        z1i2[i] = cArray[i];
    }

    FastPauli::X_r(i1x2, numQubits, targetQubit);
    for (i = 0; i < cArrLen; i++) {
        x1x2[i] = i1x2[i];
        y1x2[i] = i1x2[i];
        z1x2[i] = i1x2[i];
    }
    FastPauli::X_r(x1x2, numQubits, controlQubit);
    FastPauli::Y_r(y1x2, numQubits, controlQubit);
    FastPauli::Z_r(z1x2, numQubits, controlQubit);

    FastPauli::Y_r(i1y2, numQubits, targetQubit);
    for (i = 0; i < cArrLen; i++) {
        x1y2[i] = i1y2[i];
        y1y2[i] = i1y2[i];
        z1y2[i] = i1y2[i];
    }
    FastPauli::X_r(x1y2, numQubits, controlQubit);
    FastPauli::Y_r(y1y2, numQubits, controlQubit);
    FastPauli::Z_r(z1y2, numQubits, controlQubit);

    FastPauli::Z_r(i1z2, numQubits, targetQubit);
    for (i = 0; i < cArrLen; i++) {
        x1z2[i] = i1z2[i];
        y1z2[i] = i1z2[i];
        z1z2[i] = i1z2[i];
    }
    FastPauli::X_r(x1z2, numQubits, controlQubit);
    FastPauli::Y_r(y1z2, numQubits, controlQubit);
    FastPauli::Z_r(z1z2, numQubits, controlQubit);

    for (i = 0; i < cArrLen; i++) {
        cArray[i] = (1 - err) * cArray[i] + (err / 15.0) * (i1x2[i] +
                                                            i1y2[i] +
                                                            i1z2[i] +
                                                            x1x2[i] +
                                                            x1y2[i] +
                                                            x1z2[i] +
                                                            x1i2[i] +
                                                            y1x2[i] +
                                                            y1y2[i] +
                                                            y1z2[i] +
                                                            y1i2[i] +
                                                            z1x2[i] +
                                                            z1y2[i] +
                                                            z1z2[i] +
                                                            z1i2[i]);
    }

    delete[] i1x2;
    delete[] i1y2;
    delete[] i1z2;
    delete[] x1x2;
    delete[] x1y2;
    delete[] x1z2;
    delete[] x1i2;
    delete[] y1x2;
    delete[] y1y2;
    delete[] y1z2;
    delete[] y1i2;
    delete[] z1x2;
    delete[] z1y2;
    delete[] z1z2;
    delete[] z1i2;
}

void NoisyPauli::NCZ_r(double *cArray, const int numQubits, const int controlQubit, const int targetQubit, const double err) {
    int i;
    const int cArrLen = (UINT64_C(1) << (numQubits - 1)) + (UINT64_C(1) << (2 * numQubits - 1));  // numQubits cannot be too high, otherwise overflow and cArrLen becomes negative!
    FastPauli::CZ_r(cArray, numQubits, controlQubit, targetQubit);

    double *i1x2 = new double[cArrLen];
    double *i1y2 = new double[cArrLen];
    double *i1z2 = new double[cArrLen];
    double *x1x2 = new double[cArrLen];
    double *x1y2 = new double[cArrLen];
    double *x1z2 = new double[cArrLen];
    double *x1i2 = new double[cArrLen];
    double *y1x2 = new double[cArrLen];
    double *y1y2 = new double[cArrLen];
    double *y1z2 = new double[cArrLen];
    double *y1i2 = new double[cArrLen];
    double *z1x2 = new double[cArrLen];
    double *z1y2 = new double[cArrLen];
    double *z1z2 = new double[cArrLen];
    double *z1i2 = new double[cArrLen];

    for (i = 0; i < cArrLen; i++) {
        i1x2[i] = cArray[i];
        i1y2[i] = cArray[i];
        i1z2[i] = cArray[i];
        x1i2[i] = cArray[i];
        y1i2[i] = cArray[i];
        z1i2[i] = cArray[i];
    }

    FastPauli::X_r(i1x2, numQubits, targetQubit);
    for (i = 0; i < cArrLen; i++) {
        x1x2[i] = i1x2[i];
        y1x2[i] = i1x2[i];
        z1x2[i] = i1x2[i];
    }
    FastPauli::X_r(x1x2, numQubits, controlQubit);
    FastPauli::Y_r(y1x2, numQubits, controlQubit);
    FastPauli::Z_r(z1x2, numQubits, controlQubit);

    FastPauli::Y_r(i1y2, numQubits, targetQubit);
    for (i = 0; i < cArrLen; i++) {
        x1y2[i] = i1y2[i];
        y1y2[i] = i1y2[i];
        z1y2[i] = i1y2[i];
    }
    FastPauli::X_r(x1y2, numQubits, controlQubit);
    FastPauli::Y_r(y1y2, numQubits, controlQubit);
    FastPauli::Z_r(z1y2, numQubits, controlQubit);

    FastPauli::Z_r(i1z2, numQubits, targetQubit);
    for (i = 0; i < cArrLen; i++) {
        x1z2[i] = i1z2[i];
        y1z2[i] = i1z2[i];
        z1z2[i] = i1z2[i];
    }
    FastPauli::X_r(x1z2, numQubits, controlQubit);
    FastPauli::Y_r(y1z2, numQubits, controlQubit);
    FastPauli::Z_r(z1z2, numQubits, controlQubit);

    for (i = 0; i < cArrLen; i++) {
        cArray[i] = (1 - err) * cArray[i] + (err / 15.0) * (i1x2[i] +
                                                            i1y2[i] +
                                                            i1z2[i] +
                                                            x1x2[i] +
                                                            x1y2[i] +
                                                            x1z2[i] +
                                                            x1i2[i] +
                                                            y1x2[i] +
                                                            y1y2[i] +
                                                            y1z2[i] +
                                                            y1i2[i] +
                                                            z1x2[i] +
                                                            z1y2[i] +
                                                            z1z2[i] +
                                                            z1i2[i]);
    }

    delete[] i1x2;
    delete[] i1y2;
    delete[] i1z2;
    delete[] x1x2;
    delete[] x1y2;
    delete[] x1z2;
    delete[] x1i2;
    delete[] y1x2;
    delete[] y1y2;
    delete[] y1z2;
    delete[] y1i2;
    delete[] z1x2;
    delete[] z1y2;
    delete[] z1z2;
    delete[] z1i2;
}