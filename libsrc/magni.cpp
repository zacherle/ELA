#include <string>
#include <cmath>
#include <vector>
#include <algorithm>

#include "arrivals.h"

float localmagnitude(const double hypo[4], const TArrival &a)
{
// computes the local magnitude of a seismic event	
// formula Scherbaum , Stoll ,1983

const float min_dist = 5.0;
float hypodist;
float mag=-9.9;

if (a.phase.length() == 0) {
    return -9.9;
}
if (a.freq <= 0.0) {
    return -9.9;
}
if (a.amp <= 0.0) {
    return -9.9;
}

//define the list of phases
std::vector<std::string> phlist{"S", "Sg", "Sn", "Lg"};

//check if r.phase ocurs in the list of phases
if (std::find(std::begin(phlist), std::end(phlist), a.phase) != std::end(phlist)){
   hypodist = sqrt((a.X - hypo[0]) * (a.X - hypo[0])
	         + (a.Y - hypo[1]) * (a.Y - hypo[1])
	         + (a.Z - hypo[2]) * (a.Z - hypo[2]));
   if (hypodist > min_dist) {
       mag = log10(a.amp / 6.283 / a.freq * 2.8e6 / 0.6325) + 0.1 + 1.4 * log10(hypodist);
   }
}

return mag;
}
