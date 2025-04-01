//---------------------------------------------------------------------
// Module Arrivals
//---------------------------------------------------------------------
// All code below corresponds to the Fortran module Arrivals

// type, public :: cRecordArr
class cRecordArr {
public:
    int id = 0;
    std::string sta = "";
    std::string phase = "";
    double trec = 0.0;
    double wt = 0.0;
    double amp = 0.0;
    double freq = 0.0;
    int key = 0;
    double X = 0.0;
    double Y = 0.0;
    double Z = 0.0;
    // No procedures defined in Fortran for cRecordArr
};

// Forward declaration of class cFileArr for member procedure pointers if needed
class cFileArr;

// type, public :: cFileArr
class cFileArr {
public:
    int n_arr = 0;
    datetime reftime;
    double sumwt = 0.0;
    double sumwt2 = 0.0;
    double maxwt = 0.0;
    std::vector<cRecordArr> arr;
    
    // procedure :: fill => cfilearr_fill
    void fill(const cFileHyp &hyp);
    // procedure :: toDateTime => cfilearr_todatetime
    void toDateTime(double rtime, int &year, int &month, int &day, int &hour, int &minute, int &isec, int &msec);
};

// Global objects analogous to Fortran targets and pointers
extern cFileArr arrs;
extern cFileArr* parrs;


