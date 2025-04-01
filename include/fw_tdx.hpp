
// fw_tdx namespace containing the translated functions and subroutines
namespace fw_tdx {

    // function Xo(x,zmin)
    // transform from bounded parameter space Xp to optimization procedure space Xo
    std::array<double,4> Xo(const std::array<double,4>& x, double zmin);

    // function Xp(x,zmin)
    // transform from optimization procedure space Xo to bounded parameter space Xp
    std::array<double,4> Xp(const std::array<double,4>& x, double zmin);

    // function dpdo(z_o)
    // the derivative of the z-coordinate transformation 
    // df(y)/dx = df(y)/dy * dy(x)/dx
    // x = z_o, y = z_p(z_o) = sqrt(z_o^2+1)-1
    double dpdo(double z_o);

    // subroutine tdxder(m, n, x, fvec, fjac, ldfjac, iflag)
    void tdxder(int m, int n, const std::vector<double>& x, std::vector<double>& fvec, std::vector<std::vector<double>>& fjac, int ldfjac, int iflag);

    // subroutine tdx3(m, n, x, fvec, fjac, ldfjac, iflag)
    void tdx3(int m, int n, const std::vector<double>& x, std::vector<double>& fvec, std::vector<std::vector<double>>& fjac, int ldfjac, int iflag);

    // subroutine get_topt(m, n, x)
    void get_topt(int m, int n, std::vector<double>& x);

    // subroutine get_res_w(m, n, x, fvec, w)
    void get_res_w(int m, int n, const std::array<double,4>& x, std::vector<double>& fvec, std::vector<double>& w);

    // subroutine get_G_w(m, n, x, fjac, w)
    void get_G_w(int m, int n, const std::array<double,4>& x, std::vector<std::vector<double>>& fjac, std::vector<double>& w);

} // end namespace fw_tdx

