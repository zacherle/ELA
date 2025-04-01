#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <ctime>
#include <algorithm>

// ---------------- Dependencies and Type Definitions ----------------

// Namespace version (from module version)
namespace version {
    const std::string pversion = "ELA v1.0";
}

// Namespace layers (from module layers)
namespace layers {
    const double model_error = 0.123;
    const double reading_error = 0.456;
    const std::string name_model = "DefaultModel";
}

// Namespace arrivals: definitions of cRecordArr and cFileArr
namespace arrivals {

    struct cRecordArr {
        float trec;
        float amp;
        float freq;
        float wt;
        std::string phase;
        int id;
        float X;
        float Y;
        float Z;
    };

    struct cFileArr {
        int n_arr;
        std::vector<cRecordArr> arr;
        // Method equivalent to Fortran's todatetime
        void todatetime(float time_val, int &year, int &month, int &day,
                         int &hour, int &minute, int &isec, int &msec) const {
            // Interpret time_val as time in seconds since epoch (dummy implementation)
            std::time_t t = static_cast<std::time_t>(time_val);
            std::tm *ptm = std::localtime(&t);
            if(ptm) {
                year = ptm->tm_year + 1900;
                month = ptm->tm_mon + 1;
                day = ptm->tm_mday;
                hour = ptm->tm_hour;
                minute = ptm->tm_min;
                isec = ptm->tm_sec;
                msec = 0; // No millisecond information available
            } else {
                year = month = day = hour = minute = isec = msec = 0;
            }
        }
    };

} // end namespace arrivals

// Namespace hypfile: definitions of cRecordHyp and cFileHyp
namespace hypfile {

    struct cRecordHyp {
        std::string sta;
        std::string phase;
        float wt;
    };

    struct cFileHyp {
        std::string filename;
        std::vector<cRecordHyp> recs;
        // Function to mimic hyp%rec(a%id) (Fortran 1-indexed)
        cRecordHyp rec(int id) const {
            if(id > 0 && id <= static_cast<int>(recs.size()))
                return recs[id-1];
            else
                return cRecordHyp{};
        }
    };

} // end namespace hypfile

// Constant conversion factors (from module const_raddeg)
const double RAD2DEG = 57.295779513082320876; // 180/pi
const double DEG2RAD = 0.01745329251994329577;

// Forward declarations for external functions (from INTERFACE block)

// maxgap: real(C_FLOAT) function maxgap(naz, az) BIND(C)
extern "C" float maxgap(int naz, const float* az)
{
    if(naz <= 0) return 0.0f;
    std::vector<float> arr(az, az + naz);
    std::sort(arr.begin(), arr.end());
    float max_gap = 0.0f;
    for (int i = 1; i < naz; i++) {
         float gap = arr[i] - arr[i-1];
         if(gap > max_gap) max_gap = gap;
    }
    float circular_gap = 360.0f - arr.back() + arr.front();
    if(circular_gap > max_gap) max_gap = circular_gap;
    return max_gap;
}

// localmagnitude: function localmagnitude(hypo, r)
float localmagnitude(const float hypo[3], const arrivals::cRecordArr &r)
{
    // Dummy implementation: simply return the amplitude from record r.
    return r.amp;
}

// Function to simulate time breakdown similar to Fortran ltime
void ltime(int stime, int ita[9])
{
    std::time_t t = stime;
    std::tm *ptm = std::localtime(&t);
    if(ptm) {
        // Fortran ltime: ita(6)=year, ita(5)=month-1, ita(4)=day, ita(3)=hour, ita(2)=minute, ita(1)=second
        ita[5] = ptm->tm_year + 1900;
        ita[4] = ptm->tm_mon;  // months from 0 to 11
        ita[3] = ptm->tm_mday;
        ita[2] = ptm->tm_hour;
        ita[1] = ptm->tm_min;
        ita[0] = ptm->tm_sec;
        // Remaining elements set to zero
        for (int i = 6; i < 9; i++) {
            ita[i] = 0;
        }
    } else {
        for(int i=0; i<9; i++)
            ita[i] = 0;
    }
}

// Function XY2FL: Converts local coordinates to geographic (dummy implementation)
void XY2FL(double X, double Y, double &fi, double &rla)
{
    // Dummy conversion: simply scale the coordinates to get angles
    fi = X / 100000.0;
    rla = Y / 100000.0;
}

// Function get_errellipse: calculates error ellipse parameters (dummy implementation)
void get_errellipse(const double co[4][4], double &dxer, double &dyer, double &dzer,
                    double &dter, double &l1, double &l2, double &theta)
{
    // Dummy implementation using some of the co matrix elements
    dxer = (std::fabs(co[0][0]) < 0.001 ? 0.0 : co[0][0]);
    dyer = (std::fabs(co[1][1]) < 0.001 ? 0.0 : co[1][1]);
    dzer = (std::fabs(co[2][2]) < 0.001 ? 0.0 : co[2][2]);
    dter = (std::fabs(co[3][3]) < 0.001 ? 0.0 : co[3][3]);
    l1 = dxer + 1.0;
    l2 = dyer + 1.0;
    theta = 45.0;
}

// Function collect_obs_simul: collects simulated observations (dummy implementation)
namespace gather_module {
    void collect_obs_simul(const float hypo[4], int n_arr, std::vector<std::vector<double>> &g, int &ngather)
    {
        // Allocate a 7 x n_arr array (g has 7 rows)
        g.resize(7);
        for(auto &row : g)
            row.resize(n_arr, 0.0);
        // Dummy fill for simulation purposes.
        for (int i = 0; i < n_arr; i++) {
            // Fortran: g(1,i) -> row 1, so index 0 in C++
            g[0][i] = 0.5;  // calibration offset (dummy value)
            // g(5,i) -> row 5 -> index 4 and g(6,i) -> row 6 -> index 5, and g(7,i) -> row 7 -> index 6
            g[4][i] = hypo[3] + 0.2; // simulate calculated time offset
            g[5][i] = 0.05; // weight (dummy)
            g[6][i] = 0.7;  // atoa (dummy)
            // The other rows remain as 0.0
        }
        ngather = n_arr;
    }
} // end namespace gather_module

// ---------------- Module output_list (translated) ----------------

namespace output_list {

    // Global variable: character(LEN=255), save :: name_output_hy3 = ' '
    std::string name_output_hy3 = std::string(255, ' ');

    // Internal function mconvergence (translated from Fortran contains)
    double mconvergence(double X, double Y)
    {
        // X,Y in kilometers
        return 0.008257 * Y + 2.373 * Y / X;
    }

    // Subroutine write_hy3 translated from Fortran
    // Note: The parameter 'lu' is translated to an output std::ostream&
    void write_hy3(std::ostream &lu, const float hypo[4], const float startpoint[4],
                   int marr, const arrivals::cFileArr &arrs, const hypfile::cFileHyp &hyp,
                   const std::string &modfn, const double co[4][4], int info)
    {
        // Local variable declarations
        arrivals::cRecordArr a;
        hypfile::cRecordHyp r;

        std::string whole_date;
        int ita[9];
        int stime, i, j;
        int isec, msec;
        double dxer, dyer, dzer, dter, l1, l2, theta;
        bool ee_nan[4] = {false, false, false, false};
        double az_theta;
        double rmsres;

        std::vector<float> az(marr, 0.0f);
        float gap;

        double fi, rla;
        double meridian_con;

        int year_ref, month_ref, day_ref, hour_ref, minute_ref;
        double t_ref;
        int year_orig, month_orig, day_orig, hour_orig, minute_orig;

        std::vector<std::vector<double>> g; // 7 x marr matrix
        int ngather;
        double aobs, acal, ares, depi, dhypo, aaz, amag, dx, dy, dz, atoa;
        int nmag;
        double smag, smag2;

        // Call to arrs%todatetime for origin time using hypo(4)
        arrs.todatetime(hypo[3], year_orig, month_orig, day_orig, hour_orig, minute_orig, isec, msec);
        year_orig = year_orig % 100;
        {
            // Formatting whole_date with format: (2(i2.2,"-"), i2.2, 2x, 2(i2.2,":"), i2.2,".", i3.3)
            std::ostringstream oss;
            oss << std::setw(2) << std::setfill('0') << year_orig << "-"
                << std::setw(2) << std::setfill('0') << month_orig << "-"
                << std::setw(2) << std::setfill('0') << day_orig << "  "
                << std::setw(2) << std::setfill('0') << hour_orig << ":"
                << std::setw(2) << std::setfill('0') << minute_orig << ":"
                << std::setw(2) << std::setfill('0') << isec << "."
                << std::setw(3) << std::setfill('0') << msec;
            whole_date = oss.str();
        }

        t_ref = 0.0;
        int dummy1, dummy2, dummy3, dummy4, dummy5, dummy6;
        arrs.todatetime(static_cast<float>(t_ref), year_ref, month_ref, day_ref, hour_ref, minute_ref, dummy1, dummy2);
        year_ref = year_ref % 100;
        // --------------------------------------------------------------------

        // Write program and model parameters
        lu << "program       :ELA, " << version::pversion << std::endl;
        // The following line was commented in Fortran:
        // write (lu,'("model         :",a)') modfn(1:lnblnk(modfn))
        lu << "model         :" << layers::name_model << std::endl;
        lu << "model error   :" << std::fixed << std::setprecision(3) << layers::model_error << " s" << std::endl;
        lu << "reading error :" << std::fixed << std::setprecision(3) << layers::reading_error << " s" << std::endl;

        // Get system time
        stime = static_cast<int>(std::time(nullptr));
        ltime(stime, ita);
        // Using Fortran order: mod(ita(6),100), ita(5)+1, ita(4), ita(3), ita(2), ita(1)
        int create_year = ita[5] % 100;
        int create_month = ita[4] + 1;
        int create_day = ita[3];
        int create_hour = ita[2];
        int create_minute = ita[1];
        int create_sec = ita[0];
        lu << "create time   :" 
           << std::setw(2) << std::setfill('0') << create_year << "-"
           << std::setw(2) << std::setfill('0') << create_month << "-"
           << std::setw(2) << std::setfill('0') << create_day << " "
           << std::setw(2) << std::setfill('0') << create_hour << ":"
           << std::setw(2) << std::setfill('0') << create_minute << ":"
           << std::setw(2) << std::setfill('0') << create_sec << std::endl;
        // Write event information
        {
            std::string filename = hyp.filename;
            // Remove trailing blanks (simulate lnblnk)
            filename.erase(filename.find_last_not_of(" ") + 1);
            lu << "event         :" << filename << std::endl;
        }
        // Write startpoint
        lu << "start(x,y,z,t):(" 
           << std::fixed << std::setprecision(2)
           << startpoint[0] << ","
           << startpoint[1] << ","
           << startpoint[2] << ","
           << startpoint[3] << ")" << std::endl;
        // Write reference time
        lu << "reference time:" 
           << std::setw(2) << std::setfill('0') << year_ref << "-"
           << std::setw(2) << std::setfill('0') << month_ref << "-"
           << std::setw(2) << std::setfill('0') << day_ref << " "
           << std::setw(2) << std::setfill('0') << hour_ref << ":"
           << std::setw(2) << std::setfill('0') << minute_ref << std::endl;

        // header for station data
        lu << "---------------------------------------------------------------" << std::endl;
        lu << " sta     |obs. t.|cal. t.|res. |amplitude|freq|w| epi |hypo |azm|ain|xmag" << std::endl;
        lu << "         |  [s]  |  [s]  | [s] |  [m/s]  |[Hz]| |[km] |[km] |[o]|[o]|    " << std::endl;
        lu << "---------------------------------------------------------------" << std::endl;

        // Calculate meridian convergence
        meridian_con = mconvergence(hypo[0], hypo[1]);

        // station data: call collect_obs_simul
        gather_module::collect_obs_simul(hypo, arrs.n_arr, g, ngather);

        nmag = 0;
        smag = 0.0;
        smag2 = 0.0;
        j = 0;
        // Loop over stations (Fortran: do i=1,arrs%n_arr; here using 0-indexed)
        for (i = 0; i < arrs.n_arr; i++) {
            a = arrs.arr[i];
            r = hyp.rec(a.id);
            aobs = static_cast<double>(a.trec);
            dx = static_cast<double>(hypo[0]) - a.X;
            dy = static_cast<double>(hypo[1]) - a.Y;
            dz = static_cast<double>(hypo[2]) - a.Z;
            dhypo = std::sqrt(dx*dx + dy*dy + dz*dz);
            depi = std::sqrt(dx*dx + dy*dy);
            aaz = std::fmod(720.0 + std::atan2(dy, dx) * RAD2DEG - 180.0 - meridian_con, 360.0);
            if(aaz < 0) aaz += 360.0;
            if (g[0][i] > 0.0) {
                acal = g[0][i] + hypo[3];
                ares = aobs - acal;
                if(a.wt > 0.0) {
                    j = j + 1;
                    if(j <= static_cast<int>(az.size()))
                        az[j-1] = static_cast<float>(aaz);
                }
                {
                    float hypo3[3] = {hypo[0], hypo[1], hypo[2]};
                    amag = localmagnitude(hypo3, a);
                }
                if (amag > -9.9) {
                    nmag = nmag + 1;
                    smag += amag;
                    smag2 += amag * amag;
                }
                atoa = g[6][i];
                // Write formatted line using format 100:
                // Format: (a5, ' ', a3, '|', f7.2, '|', f7.2, '|', f5.3, '|', 1pe9.2, '|', 0pf4.1, '|', i1, '|', 0pf5.1, '|', 0pf5.1, '|', i3, '|', i3, '|')
                lu << std::setw(5) << std::left << r.sta << " " 
                   << std::setw(3) << std::left << r.phase << "|"
                   << std::right << std::fixed << std::setprecision(2) << std::setw(7) << aobs << "|"
                   << std::fixed << std::setprecision(2) << std::setw(7) << acal << "|"
                   << std::fixed << std::setprecision(3) << std::setw(5) << ares << "|"
                   << std::scientific << std::setprecision(2) << std::setw(9) << a.amp << "|"
                   << std::fixed << std::setprecision(1) << std::setw(4) << a.freq << "|"
                   << std::setw(1) << static_cast<int>(std::lround(r.wt)) << "|"
                   << std::fixed << std::setprecision(1) << std::setw(5) << depi << "|"
                   << std::fixed << std::setprecision(1) << std::setw(5) << dhypo << "|"
                   << std::setw(3) << static_cast<int>(std::lround(aaz)) << "|"
                   << std::setw(3) << static_cast<int>(std::lround(atoa)) << "|";
                if (!a.phase.empty() && (a.phase[0]=='S' || a.phase[0]=='s' || a.phase[0]=='L')) {
                    lu << std::fixed << std::setprecision(1) << std::setw(4) << amag << std::endl;
                } else {
                    lu << std::endl;
                }
            } else {
                // Else block for t_cal < 0, using format 101:
                // Format: (a5,' ',a3,'|', f7.2, '|',7X, '|',5X, '|', 1pe9.2, '|', 0pf4.1, '|', i1, '|', 0pf5.1, '|', 0pf5.1, '|', i3, '|',3X,'|')
                lu << std::setw(5) << std::left << r.sta << " " 
                   << std::setw(3) << std::left << r.phase << "|"
                   << std::right << std::fixed << std::setprecision(2) << std::setw(7) << aobs << "|"
                   << std::string(7, ' ') << "|"
                   << std::string(5, ' ') << "|"
                   << std::scientific << std::setprecision(2) << std::setw(9) << a.amp << "|"
                   << std::fixed << std::setprecision(1) << std::setw(4) << a.freq << "|"
                   << std::setw(1) << static_cast<int>(std::lround(r.wt)) << "|"
                   << std::fixed << std::setprecision(1) << std::setw(5) << depi << "|"
                   << std::fixed << std::setprecision(1) << std::setw(5) << dhypo << "|"
                   << std::setw(3) << static_cast<int>(std::lround(aaz)) << "|"
                   << std::string(3, ' ') << "|" << std::endl;
            }
        } // end do loop

        gap = maxgap(j, az.data());
        // Calculate rmsres = sqrt(sum(((g(5,:)-g(1,:)-hypo(4))*g(6,:))**2) / sum(g(6,:)))
        double sum_num = 0.0;
        double sum_w = 0.0;
        for (i = 0; i < arrs.n_arr; i++) {
            double diff = g[4][i] - g[0][i] - hypo[3]; // g(5,:) - g(1,:) - hypo(4)
            double prod = diff * g[5][i]; // multiplied by g(6,:) which is row 6 -> index 5
            sum_num += prod * prod;
            sum_w += g[5][i];
        }
        if(sum_w != 0.0)
            rmsres = std::sqrt(sum_num / sum_w);
        else
            rmsres = 0.0;

        // error ellipse: call get_errellipse(co, dxer, dyer, dzer, dter, l1, l2, theta)
        get_errellipse(co, dxer, dyer, dzer, dter, l1, l2, theta);
        ee_nan[0] = true;
        ee_nan[1] = true;
        ee_nan[2] = true;
        ee_nan[3] = true;
        if (dxer > 0 && dxer < 100) ee_nan[0] = false;
        if (dyer > 0 && dyer < 100) ee_nan[1] = false;
        if (dzer > 0 && dzer < 50) ee_nan[2] = false;
        if (dter > 0 && dter < 5) ee_nan[3] = false;

        // coordinates local --> Krovak
        // coordinates Krovak --> geographic
        az_theta = theta - 180.0 - meridian_con;
        az_theta = std::fmod(360.0 + az_theta, 360.0);
        if (az_theta < 0) az_theta += 360.0;

        // XY2FL: convert coordinates to geographic. Note: Fortran calls XY2FL (hypo(2)*1000, hypo(1)*1000,...)
        XY2FL(hypo[1] * 1000.0, hypo[0] * 1000.0, fi, rla);

        lu << std::endl << "hypocenter data:" << std::endl
           << "----------------" << std::endl;
        lu << "origin time          t:  " << std::setw(22) << whole_date << "$" << std::endl;
        if(ee_nan[3]) {
            lu << "                     +-    NaN" << std::endl;
        } else {
            lu << "                     +- " << std::fixed << std::setprecision(3) << dter << std::endl;
        }
        lu << "x-coordinate         x:  " << std::fixed << std::setprecision(2) << hypo[0] << "$" << std::endl;
        if(ee_nan[0]) {
            lu << "                     +-    NaN   km$" << std::endl;
        } else {
            lu << "                     +- " << std::fixed << std::setprecision(2) << dxer << "   km$" << std::endl;
        }
        lu << std::setw(7) << " " << "(fi:" << std::fixed << std::setprecision(6) << fi << " deg)" << std::endl;
        lu << "y-coordinate         y:  " << std::fixed << std::setprecision(2) << hypo[1] << "$" << std::endl;
        if(ee_nan[1]) {
            lu << "                     +-    NaN   km$" << std::endl;
        } else {
            lu << "                     +- " << std::fixed << std::setprecision(2) << dyer << "   km$" << std::endl;
        }
        lu << std::setw(3) << " " << "(lambda:" << std::fixed << std::setprecision(6) << rla << " deg)" << std::endl;
        lu << "depth                z:  " << std::fixed << std::setprecision(2) << hypo[2] << "$" << std::endl;
        if(ee_nan[2]) {
            lu << "                     +-    NaN   km" << std::endl;
        } else {
            lu << "                     +- " << std::fixed << std::setprecision(2) << dzer << "   km" << std::endl;
        }
        if(nmag != 0) {
            lu << "magnitude           ml:" << std::fixed << std::setprecision(2)
               << (smag / nmag) << " +- " << std::fixed << std::setprecision(2)
               << std::sqrt((smag2 - smag * smag / nmag) / nmag) << std::endl;
        } else {
            lu << "magnitude           ml: NaN" << std::endl;
        }
        lu << "rms of time residuals :" << std::setw(8) << std::fixed << std::setprecision(2)
           << rmsres << std::setw(8) << "s" << std::endl;
        lu << "angular gap           :" << std::setw(11) << static_cast<int>(std::lround(gap))
           << std::setw(8) << "deg" << std::endl;
        lu << "info                  :" << std::setw(11) << info << std::endl;

        if(ee_nan[0] || ee_nan[1]) {
            lu << "error ellipse axis l1 :" << std::setw(10) << " NaN" << std::setw(8) << "km" << std::endl;
            lu << "              axis l2 :" << std::setw(10) << " NaN" << std::setw(8) << "km" << std::endl;
        } else {
            lu << "error ellipse axis l1 :" << std::setw(8) << std::fixed << std::setprecision(2) << l1
               << std::setw(8) << "km" << std::endl;
            lu << "              axis l2 :" << std::setw(8) << std::fixed << std::setprecision(2) << l2
               << std::setw(8) << "km" << std::endl;
        }

        if(ee_nan[0] || ee_nan[1]) {
            lu << "              theta   :  NaN deg (to grid)$" << std::endl;
            lu << "    (azimuth:  NaN deg)" << std::endl;
        } else {
            lu << "              theta   :  " << std::fixed << std::setprecision(1) << theta
               << " deg (to grid)$" << std::endl;
            lu << "    (azimuth: " << std::fixed << std::setprecision(1) << az_theta << " deg)" << std::endl;
        }
        return;
    } // end write_hy3

} // end namespace output_list

// ---------------- Main Function (for testing purposes) ----------------

int main()
{
    // Dummy input parameters matching Fortran subroutine signature

    // Output stream (lu) simulated by std::cout
    std::ostream &lu = std::cout;

    // Hypocenter and startpoint arrays (4 elements each)
    float hypo[4] = {34.05f, -118.25f, 10.0f, 0.5f};       // x, y, z, t (dummy values)
    float startpoint[4] = {33.95f, -118.35f, 5.0f, 0.4f};

    // marr: number of stations for azimuth array
    int marr = 10;

    // Create a dummy arrivals::cFileArr object with some station records
    arrivals::cFileArr arrs;
    arrs.n_arr = 3;
    arrs.arr.resize(arrs.n_arr);
    // Fill dummy station records
    arrs.arr[0] = {1.23f, 2.5f, 5.6f, 1.0f, "Sphase", 1, 34.00f, -118.20f, 12.0f};
    arrs.arr[1] = {2.34f, 3.5f, 6.7f, 1.0f, "Lphase", 2, 34.10f, -118.30f, 15.0f};
    arrs.arr[2] = {3.45f, 4.5f, 7.8f, 1.0f, "Pphase", 3, 34.05f, -118.25f, 8.0f};

    // Create a dummy hypfile::cFileHyp object with dummy records and a filename
    hypfile::cFileHyp hyp;
    hyp.filename = "dummy_event.dat   "; // trailing blanks simulated
    hyp.recs.resize(3);
    hyp.recs[0] = {"STA1", "S", 1.0f};
    hyp.recs[1] = {"STA2", "L", 1.0f};
    hyp.recs[2] = {"STA3", "P", 1.0f};

    // modfn: model file name (dummy)
    std::string modfn = "model_file.mod";

    // co: 4x4 covariance matrix (dummy values)
    double co[4][4] = {
        {1.2, 0,   0,   0},
        {0,   2.3, 0,   0},
        {0,   0,   3.4, 0},
        {0,   0,   0,   4.5}
    };

    // info: integer information (dummy)
    int info = 99;

    // Call the translated subroutine write_hy3
    output_list::write_hy3(lu, hypo, startpoint, marr, arrs, hyp, modfn, co, info);

    return 0;
}
