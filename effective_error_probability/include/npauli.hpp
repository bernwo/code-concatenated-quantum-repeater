#ifndef _NPAULIHPP_
#define _NPAULIHPP_

class NoisyPauli {
   private:
   public:
    static void SingleQubitDepolarise_r(double *cArray, const int numQubits, const int targetQubit, const double err);        // Single-qubit depolarising channel.
    static void NCX_r(double *cArray, const int numQubits, const int controlQubit, const int targetQubit, const double err);  // Two-qubit depolarising channel.
    static void NCZ_r(double *cArray, const int numQubits, const int controlQubit, const int targetQubit, const double err);  // Two-qubit depolarising channel.
};

#endif