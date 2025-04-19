#ifndef CALC_COVAR_H
#define CALC_COVAR_H

double* cov_matrix1(int m, int n, double* resw, double* w, const double* Gw);

void get_errellipse(const double* co,
	       	float &dxer, float &dyer, float &dzer, float &dter,
		float &l1, float &l2, float &theta);

#endif // CALC_COVAR_H
