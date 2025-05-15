#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <array>

#include "fw_tdx.h"


//http://devernay.free.fr/hacks/cminpack/lmder_.html
/*
void lmder1_(
	void (*fcn) (int *m, int *n, double *x, double *fvec, double *fjac,
	    	int *ldfjac, int *iflag),
        int *m, int * n, double *x, double *fvec, double *fjac, int *ldfjac, 
        double *tol, int *info, int *iwa, double *wa, int *lwa); 

void lmder_ (
	void (*fcn)( int *m, int *n, double *x, double *fvec, double *fjac,
	    	int *ldfjac, int *iflag),
        int *m, int *n, double *x, double *fvec, double *fjac, int *ldfjac, 
        double *ftol, double *xtol, double *gtol, int *maxfev, double *diag, int *mode, 
        double *factor, int *nprint, int *info, 
        int *nfev, int *njev, int *ipvt, 
        double *qtf, double *wa1, double *wa2, double *wa3, double *wa4 ); 

void fcn(int* m, int* n, double *x, double *fvec, double *fjac,
        int *ldfjac, int *iflag)
*/
extern "C" {
// Define a function pointer type corresponding to the Fortran procedure()
// signature used by lmder1, tdx3, and tdxder.
typedef void (*FcnType)(int*, int*, const double*, double*, double*, int*, int*);

void lmder1_(FcnType fcn,
        int *m, int * n, double *x, double *fvec, double *fjac, int *ldfjac, 
        double *tol, int *info, int *iwa, double *wa, int *lwa); 
} // extern "C"

//void lmder1(FcnType fcn, int* m, int n, std::array<double, 4>& x, std::vector<double>& Fvec,
//            std::vector<std::vector<double>>& Fjac, int Ldfjac, double Tol, int& Info,
//            std::vector<int>& Ipvt, std::vector<double>& Wa, int Lwa);


void loc_hypo_lm(std::array<double, 4>& hypo, int m, double* fvec, 
                 double* fjac, bool fix_depth, int& info) {
    
    int n = 4;
    std::array<double, 4> x;
    int ipvt[n];
    double tol = std::sqrt(std::numeric_limits<double>::epsilon());
    double tol100 = 100.0 * tol;
    // tol = 1.0e-6; // 1.0e-6 is the default value in the original code 
    int fvec_size = m;
    int fjac_size = m;
    if (fix_depth) {
        int wa_size = 5*n+m;
        std::vector<double> wa;
        wa.resize(wa_size);
        x = hypo;
        lmder1_(tdx3, &m, &n, (double*) x.data(), fvec, fjac, &fjac_size, &tol100,
               &info, ipvt, (double*) wa.data(), &wa_size);
    } else {
        int wa_size = m*n + 5*n + m;
        std::vector<double> wa;
        wa.resize(wa_size);
        x = Xo(hypo.data(), 0.0); // relief = 0.0
        lmder1_(tdxder, &m, &n, (double*) x.data(), fvec, fjac, &fjac_size, &tol,
               &info, ipvt, (double*) wa.data(), &wa_size);
        x = Xp(x.data(), 0.0);
    }
    hypo = x;
}

