#ifndef RT1DLAYERS_H
#define RT1DLAYERS_H

#ifdef __cplusplus
extern"C" {
#endif

void raytr1D(double s_x, double s_z,
	       	int nl, float* d, float* v,
	       	double tt, double raypar, int type_of_wave);
		 
#ifdef __cplusplus
}
#endif
#endif // RT1DLAYERS_H
