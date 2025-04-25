#ifndef RT1DLAYERS_H
#define RT1DLAYERS_H

#ifdef __cplusplus
extern"C" {
#endif

void raytr1D(const double s_x, const double s_z,
	       	const int nl, const float d[], const float v[],
	       	double* tt, double* raypar, int* type_of_wave);
		 
#ifdef __cplusplus
}
#endif
#endif // RT1DLAYERS_H
