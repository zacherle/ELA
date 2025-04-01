
double get1Dvel(double depth, int nl, const std::vector<double>& d, const std::vector<double>& v);

void two_point_td(const std::vector<double>& srcxyz, const std::vector<double>& obsxyz, const std::string& ph, std::vector<double>& td, double& toa);
