#ifndef FW_TDX_H
#define FW_TDX_H

#include <array>

// transform from bounded parameter space Xp to optimization procedure space Xo
std::array<double,4> Xo(const double* x, double zmin);

// transform from optimization procedure space Xo to bounded parameter space Xp
std::array<double,4> Xp(const double* x, double zmin);

extern "C" {
void tdxder(int* m, int* n, const double* x, double* fvec, double* fjac, int* ldjac, int* iflag);

void tdx3(int* m, int* n, const double* x, double* fvec, double* fjac, int* ldjac, int* iflag);
}

#endif // FW_TDX_H
