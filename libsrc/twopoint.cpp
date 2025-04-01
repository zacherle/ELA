#include <cmath>
#include <vector>
#include <iostream>
#include <memory>

class cLayers {
public:
    int n_layers;
    std::vector<double> z;
    std::vector<double> v;
};

std::shared_ptr<cLayers> get_layers(const std::string& ph);

double get1Dvel(double depth, int nl, const std::vector<double>& d, const std::vector<double>& v) {
    double vel;
    const double vair = 0.300;

    if (depth < d[0]) {
        return vair;
    }

    vel = vair;
    for (int i = nl - 1; i >= 0; --i) {
        if (depth >= d[i]) break;
    }
    if (i >= 0) {
        vel = v[i];
    } else {
        vel = vair;
    }

    return vel;
}

void two_point_td(const std::vector<double>& srcxyz, const std::vector<double>& obsxyz, const std::string& ph, std::vector<double>& td, double& toa) {
    std::shared_ptr<cLayers> lptr = get_layers(ph);
    int nl;
    double toas;
    double toac;
    double raypar;
    int type_of_wave;
    double v_hypo;
    double delta;
    double tt;

    if (!lptr) {
        td[0] = -1.0;
        td[1] = 0.0;
        td[2] = 0.0;
        td[3] = 0.0;
        return;
    }
    nl = lptr->n_layers;
    delta = std::sqrt(std::pow(obsxyz[0] - srcxyz[0], 2) + std::pow(obsxyz[1] - srcxyz[1], 2));
    if (delta < 0.0001) {
        delta = 0.0001;
    }

    // Call to raytr1D function should be implemented
    // raytr1D(delta, srcxyz[2], nl, lptr->z, lptr->v, tt, raypar, type_of_wave);

    v_hypo = get1Dvel(srcxyz[2], nl, lptr->z, lptr->v);

    td[0] = tt;
    td[1] = -raypar * (obsxyz[0] - srcxyz[0]) / delta;
    td[2] = -raypar * (obsxyz[1] - srcxyz[1]) / delta;
    td[3] = 0.0;
    toas = raypar * v_hypo;
    toac = std::sqrt(1.0 - toas * toas);
    if (std::isnan(toac)) {
        toas = 1.0;
        toac = 0.0;
    }
    if (type_of_wave == 1) {
        toac = -toac;
    }
    td[3] = toac / v_hypo;

    toa = std::acos(-toac) * (180.0 / M_PI); // Convert radians to degrees
}
