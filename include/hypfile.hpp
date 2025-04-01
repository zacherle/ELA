
namespace HypFile {

//---------------------------------------------------------------------
// class cRecordHyp
//---------------------------------------------------------------------
class cRecordHyp {
public:
    int id;
    std::string sta;
    std::string phase;
    int ichan;
    int year;
    int month;
    int day;
    int hour;
    int minute;
    int isec;
    int msec;
    int micros;
    float wt;
    float amp_d;
    int pol;
    float amp_v;
    float period;
    // Method: getYear()
    int getYear() const {
        return year;
    }

    // Reads a record from the global input stream (iu)
    int rrec();
    // Writes a record to the global output stream (ou) using the first format.
    void wrec() const; 
    // Writes a record to the global output stream (ou) using the second format.
    void wrec_hyr() const;
};

//---------------------------------------------------------------------
// class cFileHyp
//---------------------------------------------------------------------
class cFileHyp {
public:
    int n_rec;
    bool hyr;
    std::vector<cRecordHyp> rec; // Allocatable array of records
    std::string filename = "";
    
    cFileHyp() : n_rec(0), hyr(false) {}
    // Loads the hyp file.
    int load(const std::string & hypname);
    // Dumps the hyp file.
    int dump(const std::string & hypname);
};

// global instance of cFileHyp corresponding to Fortran: type(cFileHyp),public :: hyp
extern cFileHyp hyp;

} // end namespace HypFile

