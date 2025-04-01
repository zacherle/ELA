#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

namespace iso_fortran_env {
    using real64 = double;
}

struct Arrival {
    double X;
    double Y;
    double Z;
    std::string phase;
    double trec;
    double wt;
};

struct Parrrs {
    int n_arr;
    std::vector<Arrival> arr;
};

void two_point_td(const std::string& ph, const std::vector<double>& srcxyz, const std::vector<double>& obsxyz, std::vector<double>& td, double& toa);

void collect_obs_simul(const std::vector<double>& srcxyz, int m, std::vector<std::vector<iso_fortran_env::real64>>& gather, int& ngather, Parrrs& parrs) {
    int narr = parrs.n_arr;
    if (narr > m) {
        narr = m;
    }
    int j = 0;
    std::vector<double> obsxyz(3);
    std::string ph;
    std::vector<double> td(4);
    double toa;

    for (int i = 0; i < narr; ++i) {
        obsxyz[0] = parrs.arr[i].X;
        obsxyz[1] = parrs.arr[i].Y;
        obsxyz[2] = parrs.arr[i].Z;
        ph = parrs.arr[i].phase;
        two_point_td(ph, srcxyz, obsxyz, td, toa);
        
        j++;
        gather[0][j - 1] = td[0]; // travel time
        gather[1][j - 1] = td[1]; // dt/dx
        gather[2][j - 1] = td[2]; // dt/dy
        gather[3][j - 1] = td[3]; // dt/dz
        
        if (td[0] > 0.0) {
            gather[4][j - 1] = parrs.arr[i].trec;
            gather[5][j - 1] = parrs.arr[i].wt;
            gather[6][j - 1] = toa;
        } else {
            gather[4][j - 1] = parrs.arr[i].trec;
            gather[5][j - 1] = 0.0;
            gather[6][j - 1] = toa;
        }
    }
    ngather = j;
}
