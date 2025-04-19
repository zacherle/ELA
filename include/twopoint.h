#ifndef TWOPNT_H
#define TWOPNT_H

#include <string>

void two_point_td(const double (&srcxyz)[3], const double (&obsxyz)[3],
	       	const std::string ph, double (&td)[4], double& toa);

#endif // TWOPNT_H 
