#ifndef CALC_COVAR_H
#define CALC_COVAR_H

#ifdef __cplusplus
extern "C" {
#endif

void cov_matrix2(double* cov, const int m, const int n,
	       	const double resw[], const double Gw[], const double w[],
	       	const double reading_err, const double model_err);

#ifdef __cplusplus
}
#endif
#endif // CALC_COVAR_H
