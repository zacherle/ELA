#include <cmath>
#include <vector>
#include <iostream>
#include <memory>
#include <map>
#include <algorithm>

#include "rt1Dlayers.h" // contains the definition of raytr1D
/*
struct TLayers{
    int n_layers;
    std::vector<float> z;
    std::vector<float> v;
};
*/
#include "layers.h" // contains the definition of get_layers



// vector std::vector<float> is sorted in ascending order
// use std::upper_bound and std::distance and find the index of element in the vector x
// corresponding to the float number 'value'
// This function returns the index of the largest element in the vector x that is less than or equal to value.
// If all elements in x are less than value, it returns the last index of the vector.
// If all elements in x are greater than value, it returns -1.
// If the vector is empty, it returns -1.
// If the vector contains only one element, it returns 0 if that element is less than or equal to value, otherwise it returns -1.
// If the vector contains two elements, it returns 0 if the first element is less than or equal to value and the second element is greater than value, otherwise it returns -1.
int find_index(const std::vector<float>& x, float value) {
    if (x.empty()) {
	return -1; // Return -1 for empty vector
    }
    if (x.size() == 1) {
	return (x[0] <= value) ? 0 : -1; // Return 0 or -1 for single element vector
    }

    auto it = std::upper_bound(x.begin(), x.end(), value);  //if (!(val<*it))
    int index = std::distance(x.begin(), it) - 1;
    return (index >= 0 && index < static_cast<int>(x.size())) ? index : -1;
}


double get1Dvel(double depth, int nl, const std::vector<float>& d, const std::vector<float>& v) {
    double vel;
    const double vair = 0.300;

    if (depth < d[0]) {
        return vair;
    }

    vel = vair;
    //int i;
    //for (i = nl - 1; i >= 0; --i) {
    //    if (depth >= d[i]) break;
    //}
    
    int i = find_index(d, depth);
    
    if (i >= 0) {
        vel = v[i];
    } else {
        vel = vair;
    }

    return vel;
}

void two_point_td(const double (&srcxyz)[3], const double (&obsxyz)[3],
	       	const std::string ph, double (&td)[4], double& toa) {
// Function to calculate the travel time and ray parameter for a 1D layered medium
// using the raytr1D function
    std::shared_ptr<TLayers> lptr = getLayersPointer(ph);
    int nl;
    double toas;
    double toac;
    double raypar;
    double v_hypo;
    double delta;
    double tt=-1.0;
    int type_of_wave=-1;

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

    raytr1D(delta, srcxyz[2], nl, lptr->z.data(), lptr->v.data(),
		    &tt, &raypar, &type_of_wave);

    v_hypo = get1Dvel(srcxyz[2], nl, lptr->z, lptr->v);

    td[0] = tt;
    td[1] = -raypar * (obsxyz[0] - srcxyz[0]) / delta;
    td[2] = -raypar * (obsxyz[1] - srcxyz[1]) / delta;
    td[3] = 0.0;
    toas = raypar * v_hypo;
    toac = std::sqrt(1.0 - toas*toas);
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


// main function for testing find_index and get1Dvel
// one argument is the 'value' to be searched in the vector
// test data vector x and vector v is hardcoded
// print the index of the value in the vector x
// and x[index] and v[index] to the console
// return 0
/*
int main(int argc, char* argv[]) {
    std::vector<float> x = {0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<float> v = {5.0, 6.0, 7.0, 8.0, 9.0};
    float value = 2.5;
    if (argc > 1) {
	value = std::stof(argv[1]);
    }
    int index = find_index(x, value);
    if (index >= 0) {
	std::cout << "Index: " << index << ", x[index]: " << x[index] << ", v[index]: " << v[index] << std::endl;
    } else {
	std::cout << "Value not found in the vector." << std::endl;
    }

    double vel= get1Dvel(value, x.size(), x, v);

    std::cout << "Velocity at depth " << value << ": " << vel << std::endl;

    return 0;
}
*/

