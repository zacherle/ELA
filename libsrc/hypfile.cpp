#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace HypFile {

// Global file streams corresponding to Fortran units
std::ifstream iu; // For input file unit
std::ofstream ou; // For output file unit

//---------------------------------------------------------------------
// class cRecordHyp
//---------------------------------------------------------------------
class cRecordHyp {
public:
    int id = 0;
    std::string sta = "";
    std::string phase = "";
    int ichan = 0;
    int year = 0;
    int month = 0;
    int day = 0;
    int hour = 0;
    int minute = 0;
    int isec = 0;
    int msec = 0;
    int micros = 0;
    float wt = 0;
    float amp_d = 0;
    int pol = 0;
    float amp_v = 0;
    float period = 0;

    // Reads a record from the global input stream (iu)
    // Corresponds to Fortran: function recordhyp_read(r) result(res)
    int rrec() {
        int res = 0;
        // Create a buffer of length 100 characters
        std::string cbuf;
        // Read one line from the input stream using the entire line (similar to Fortran read(iu,'(a)'))
        if (!std::getline(iu, cbuf)) {
            if (iu.eof()) {
                return -1; // res < 0: end-of-file
            } else {
                return 1;  // res > 0: error occurred
            }
        }
        // Use a string stream to parse the buffer
        std::istringstream iss(cbuf);
        // Read tokens from the input stream.
        // The tokens appear in the following order:
        // sta, phase, ichan, year, month, day, hour, minute, isec, msec, micros, wt, amp_d, pol, amp_v, period
        if (!(iss >> sta >> phase
                  >> ichan >> year >> month >> day
                  >> hour >> minute >> isec >> msec
                  >> micros >> wt >> amp_d >> pol
                  >> amp_v >> period)) {
            return 2; // Some error in conversion (non-zero error code)
        }
        return res;
    }

    // Writes a record to the global output stream (ou) using the first format.
    // Corresponds to Fortran: subroutine recordhyp_write(r)
    void wrec() const {
        // Compute mod(year,100)
        int mod_year = year % 100;
        // Buffer for formatted output
        char buffer[256];
        // Format string corresponding to Fortran format:
        // 101   format(a5,1x,a,1x,i2,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2,
        //        1x,i3.3,1x,i3.3,1x,i1,1x,1pe12.6,2x,i1,1x,1pe12.6,1x,1pe12.6)
        // Field-by-field conversion:
        //   sta      -> a5       -> %-5s
        //   1x       -> " "
        //   phase    -> a        -> %s
        //   1x       -> " "
        //   ichan    -> i2       -> %2d
        //   1x       -> " "
        //   mod(year,100) -> i2.2 -> %02d
        //   1x       -> " "
        //   month    -> i2.2   -> %02d
        //   1x       -> " "
        //   day      -> i2.2   -> %02d
        //   1x       -> " "
        //   hour     -> i2.2   -> %02d
        //   1x       -> " "
        //   minute   -> i2.2   -> %02d
        //   1x       -> " "
        //   isec     -> i2.2   -> %02d
        //   1x       -> " "
        //   msec     -> i3.3   -> %03d
        //   1x       -> " "
        //   micros   -> i3.3   -> %03d
        //   1x       -> " "
        //   int(wt)  -> i1     -> %1d
        //   1x       -> " "
        //   amp_d    -> 1pe12.6 -> %12.6e
        //   2x       -> "  "
        //   pol      -> i1     -> %1d
        //   1x       -> " "
        //   amp_v    -> 1pe12.6 -> %12.6e
        //   1x       -> " "
        //   period   -> 1pe12.6 -> %12.6e
        std::snprintf(buffer, sizeof(buffer), "%-5s %s %2d %02d %02d %02d %02d %02d %02d %03d %03d %1d %12.6e  %1d %12.6e %12.6e",
                      sta.substr(0,5).c_str(), phase.c_str(), ichan, mod_year, month, day,
                      hour, minute, isec, msec, micros, int(wt), amp_d, pol, amp_v, period);
        ou << buffer << std::endl;
    }

    // Writes a record to the global output stream (ou) using the second format.
    // Corresponds to Fortran: subroutine recordhyp_write_hyr(r)
    void wrec_hyr() const {
        int mod_year = year % 100;
        char buffer[256];
        // For recordhyp_write_hyr, the format string is:
        // 101   format(a5,1x,a,1x,i2,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2,
        //        1x,i3.3,1x,i3.3,1x,f5.0,1x,1pe12.6,2x,i1,1x,1pe12.6,1x,1pe12.6)
        // The difference is that the field for wt is printed with f5.0 instead of int(wt)
        std::snprintf(buffer, sizeof(buffer), "%-5s %s %2d %02d %02d %02d %02d %02d %02d %03d %03d %5.0f %12.6e  %1d %12.6e %12.6e",
                      sta.substr(0,5).c_str(), phase.c_str(), ichan, mod_year, month, day,
                      hour, minute, isec, msec, micros, wt, amp_d, pol, amp_v, period);
        ou << buffer << std::endl;
    }

    // Gets the full year from the stored year field.
    // Corresponds to Fortran: function recordhyp_getyear(r) result(year)
    int getYear() const {
        int mod_year = year % 100;
        int full_year = mod_year;
        if (mod_year > 70) {
            full_year = mod_year + 1900;
        } else {
            full_year = mod_year + 2000;
        }
        return full_year;
    }
};

//---------------------------------------------------------------------
// class cFileHyp
//---------------------------------------------------------------------
class cFileHyp {
public:
    int n_rec = 0;
    bool hyr = false;
    std::vector<cRecordHyp> rec; // Allocatable array of records
    std::string filename = "";

    // Loads the hyp file.
    // Corresponds to Fortran: function filehyp_load(hyp,hypname) result(rmsg)
    int load(const std::string & hypname) {
        int rmsg = 0;
        int ios = 0; // iostat equivalent
        // For iomsg we use a character buffer of size 256
        char iom[256] = {0};
        std::string cbuffer;
        cRecordHyp r;
        int i = 0;

        // Store trimmed filename (Fortran trim is simulated by std::string's erase of trailing spaces)
        filename = hypname;
        // Open the file for reading (FORMATTED, STATUS='OLD')
        // Using a new file stream for input (iu)
        iu.close();
        iu.clear();
        iu.open(filename);
        if (!iu.is_open()) {
            std::cerr << "iostat = " << 1 << std::endl;
            std::snprintf(iom, sizeof(iom), "Error opening file: %s", filename.c_str());
            std::cerr << "iomsg: " << iom << std::endl;
            rmsg = 1;
            std::exit(rmsg);
        }

        // Set hyr flag based on file extension.
        hyr = false;
        if (filename.size() >= 4) {
            std::string ext = filename.substr(filename.size() - 4);
            if (ext == ".hyr" || ext == ".HYR") {
                hyr = true;
            }
        }

        n_rec = 0;
        // First pass: count the number of records by reading lines.
        while (std::getline(iu, cbuffer)) {
            // !res == 0  ok
            n_rec++;
            // cycle (continue reading)
        }
        if (iu.bad()) {
            iu.close();
            rmsg = 2;
            std::exit(rmsg);
        }
        // Rewind the file.
        iu.clear();
        iu.seekg(0, std::ios::beg);

        // Allocate the record vector with the counted number of records.
        rec.resize(n_rec);

        // Second pass: read each record.
        for (i = 0; i < n_rec; i++) {
            reswitch:
            {
                int res = r.rrec();
                if (res == 0) { // res == 0 ok
                    r.id = i + 1;
                    rec[i] = r;
                    // cycle -> continue the loop
                    continue;
                }
                if (res < 0) { // res == -1 end of file
                    break;    // exit loop
                }
                // If res > 0 then error occurred.
                iu.close();
                rmsg = res;
                std::exit(rmsg);
            }
        }
        iu.close();
        rmsg = 0;
        return rmsg;
    }

    // Dumps the hyp file.
    // Corresponds to Fortran: function filehyp_dump(hyp,hypname) result(rmsg)
    int dump(const std::string & hypname) {
        int rmsg = 0;
        int ios = 0;
        char iom[256] = {0};
        cRecordHyp r;
        int i = 0;

        // Open the file for writing (ACTION='write',STATUS='replace')
        ou.close();
        ou.clear();
        ou.open(hypname, std::ios::out | std::ios::trunc);
        if (!ou.is_open()) {
            std::cerr << "iostat = " << 1 << std::endl;
            std::snprintf(iom, sizeof(iom), "Error opening file for write: %s", hypname.c_str());
            std::cerr << "iomsg: " << iom << std::endl;
            rmsg = 1;
            std::exit(rmsg);
        }
        // Write each record according to the hyr flag.
        for (i = 0; i < n_rec; i++) {
            r = rec[i];
            if (hyr) {
                r.wrec_hyr();  // call r%wrec_hyr()
            } else {
                r.wrec();      // call r%wrec()
            }
        }
        ou.close();
        rmsg = 0;
        return rmsg;
    }
};

// Create a global instance of cFileHyp corresponding to Fortran: type(cFileHyp),public :: hyp
cFileHyp hyp;

} // end namespace HypFile

//---------------------------------------------------------------------
// The following main program is commented out as in the original Fortran code.
// It demonstrates how to use the HypFile module.
/*
#include "HypFile.hpp" // If module were in separate file

int main() {
    int rmsg;
    // Create instances (global hyp already exists in the module)
    // For Arrivals, the Fortran code had type(cFileArr) arrs; here we do not implement it.
    // rmsg = HypFile::hyp.load("test.hyp");
    // rmsg = HypFile::hyp.dump("testo.hyp");
    // call arrs.fill(hyp);
    return 0;
}
*/
//---------------------------------------------------------------------
 
// End of complete translation.
 
int main() {
    // This main can be used for quick testing.
    // Uncomment and modify the file names as needed.
    // Note: The Arrivals module is not implemented.
    // Example usage:
    // int rmsg = HypFile::hyp.load("test.hyp");
    // rmsg = HypFile::hyp.dump("testo.hyp");
    return 0;
}
 
