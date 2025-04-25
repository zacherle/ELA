#ifndef WGSJTSK_H
#define WGSJTSK_H

#ifdef __cplusplus
extern "C" {
#endif

void wgs2jtsk (double phi, double lambda, double *x, double *y ); 

void jtsk2wgs (double X, double Y, double *phi, double *lambda );

#ifdef __cplusplus
}
#endif

#endif // WGSJTSK_H
