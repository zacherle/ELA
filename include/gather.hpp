
struct Parrrs {
    int n_arr;
    std::vector<Arrival> arr;
};

void two_point_td(const std::string& ph, const std::vector<double>& srcxyz, const std::vector<double>& obsxyz, std::vector<double>& td, double& toa);

void collect_obs_simul(const std::vector<double>& srcxyz, int m, std::vector<std::vector<double>>& gather, int& ngather, Parrrs& parrs);
