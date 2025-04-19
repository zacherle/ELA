#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <cstdlib>  // for std::abs

#include "gather.h"

extern "C" {



    // Returns |a| with the sign of b.
    inline double sign(double a, double b) {
        return (b >= 0 ? std::fabs(a) : -std::fabs(a));
    }

    // transform from bounded parameter space Xp to optimization procedure space Xo
    std::array<double,4> Xo(const double (& x)[4], double zmin) {
        //double Xo[4];
        double o;
        // reflective transformation zp -> zo
        o = std::sqrt(std::pow((x[2] - zmin + 1.0), 2.0) - 1.0);
        return { x[0] - 1100.0, x[1] - 700.0, o, x[3] };
    }

    // transform from optimization procedure space Xo to bounded parameter space Xp
    std::array<double,4> Xp(const double (& x)[4], double zmin) {
        double p;
        // reflective transformation zo -> zp
        p = x[2];
        p = zmin - 1.0 + std::sqrt(p * p + 1.0);
        return { x[0] + 1100.0, x[1] + 700.0, p, x[3] };
    }
/* alternative with static array
    // transform from optimization procedure space Xo to bounded parameter space Xp
    double (&Xp(const double (& x)[4], double zmin))[4] {
        double p;
        // reflective transformation zo -> zp
        p = x[2];
        p = zmin - 1.0 + std::sqrt(p * p + 1.0);
        static double xp[4] = { x[0] + 1100.0, x[1] + 700.0, p, x[3] };
        return xp;
    }
*/
    // function dpdo(z_o)
    // the derivative of the z-coordinate transformation 
    // df(y)/dx = df(y)/dy * dy(x)/dx
    // x = z_o, y = z_p(z_o) = sqrt(z_o^2+1)-1
    double dpdo(double z_o) {
        double r;
        r = z_o / std::sqrt(1.0 + z_o * z_o);
        // prevention of zero
        r = r + sign(0.001, z_o);
        return r;
    }


void tdxder(int& m, int& n, const double (& x)[4], double* fvec, double* fjac, int& ldjac, int& iflag) {
    // Local variables
    //double xp_[4];
    //std::vector<double> w(m, 0.0);
    double w[m];
    int mres;
    // transform from optimization x to bounded parametric xp space
    // xp_ = Xp(x, relief8(x(1),x(2)))
    // double (&xp_)[4] = Xp( x, 0.0);
    auto xp_=Xp(x, 0.0);
    if (iflag == 1) {
	mres=gather.collect(xp_.data());
        // fvec = fvec * w (elementwise multiplication)
        mres=gather.get_res_w(xp_.data(), m, fvec, w);
        for (int j = 0; j < m; j++) {
            fvec[j] = fvec[j] * w[j];
        }
    } else if (iflag == 2) {
        mres=gather.get_G_w(m, fjac, w);

        // transform dt/dz from bounded xp -> optimization xo space
        // dpdo = dz_p(z_o)/dz_o 
        for (int j = 0; j < m; j++) {
            fjac[j+2*ldjac] *= dpdo(x[2]);
        }

        // Multiply each column elementwise by w
        for (int j = 0; j < m; j++) {
            float wt = w[j];
            fjac[j+0*ldjac] *= wt;
	    fjac[j+1*ldjac] *= wt;
	    fjac[j+2*ldjac] *= wt;
	    fjac[j+3*ldjac] *= wt;
	}
    } else {
        // print xp_
        std::cout << "xp_ = ";
        for (int i = 0; i < 4; i++) {
            std::cout << xp_[i] << " ";
        }
        std::cout << std::endl;
    }
}


// subroutine tdx3 will be called by the optimization procedure lmder
// of standard fortran library MINPACK
// tdx3 is a function that computes the residuals and the Jacobian
// of the function f(x) = 0
// x = (x[0], x[1], x[2], x[3]), where x[2] is the depth fixed
// fvec = (fvec[0], fvec[1], fvec[2], fvec[3])
// fjac = (fjac[0][0], fjac[0][1], fjac[0][2], fjac[0][3])
// ldfjac = leading dimension of fjac
// iflag = 1: compute fvec
// iflag = 2: compute fjac
// m = number of residuals
// n = number of parameters

void tdx3(int& m, int& n, const double (& x)[4], double* fvec, double* fjac, int& ldjac, int& iflag) {
    //std::vector<double> w(m, 0.0);
    double w[m];
    int mres;

    if (iflag == 1) {
	mres=gather.collect(x);
	// fvec = fvec * w (elementwise multiplication)
        mres=gather.get_res_w(x, m, fvec, w);
        for (int j = 0; j < m; j++) {
            fvec[j] = fvec[j] * w[j];
        }
    } else if (iflag == 2) {
        mres=gather.get_G_w(m, fjac, w);
        for (int j = 0; j < m; j++) {
            float wt = w[j];
            fjac[j+0*ldjac] *= wt;
	    fjac[j+1*ldjac] *= wt;
	    fjac[j+2*ldjac] *= wt;   //fix depth
	    fjac[j+3*ldjac] *= wt;
        }
    } else {
        // print x
        std::cout << "x = ";
        for (size_t i = 0; i < 4; i++) {
            std::cout << x[i] << " ";
        }
        std::cout << std::endl;
    }
}
} // extern "C"
