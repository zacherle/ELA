#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <iomanip>
#include <sstream>

//---------------------------------------------------------------------
// Supporting datetime and timedelta classes (equivalent to datetime_module)
//---------------------------------------------------------------------
class timedelta; // forward declaration

class datetime {
private:
    std::chrono::system_clock::time_point tp;
public:
    // Default constructor
    datetime() {
        tp = std::chrono::system_clock::now();
    }
    // Constructor from year only (e.g., datetime(3000))
    explicit datetime(int year) {
        std::tm t = {};
        t.tm_year = year - 1900;
        t.tm_mon = 0;
        t.tm_mday = 1;
        t.tm_hour = 0;
        t.tm_min = 0;
        t.tm_sec = 0;
        std::time_t time_tt = std::mktime(&t);
        tp = std::chrono::system_clock::from_time_t(time_tt);
    }
    // Constructor from full date and time (seconds and millisecond)
    datetime(int year, int month, int day, int hour, int minute, int second, int msec=0) {
        std::tm t = {};
        t.tm_year = year - 1900;
        t.tm_mon = month - 1;
        t.tm_mday = day;
        t.tm_hour = hour;
        t.tm_min = minute;
        t.tm_sec = second;
        std::time_t time_tt = std::mktime(&t);
        tp = std::chrono::system_clock::from_time_t(time_tt) + std::chrono::milliseconds(msec);
    }
    // Copy constructor
    datetime(const datetime &other) {
        tp = other.tp;
    }
    // Overloaded assignment operator
    datetime& operator=(const datetime &other) {
        if(this != &other) {
            tp = other.tp;
        }
        return *this;
    }
    // Addition of timedelta
    datetime operator+(const timedelta &td) const;
    // Subtraction of two datetime objects gives a timedelta
    timedelta operator-(const datetime &other) const;
    // Returns ISO formatted date string
    std::string isoformat() const {
        auto time_tt = std::chrono::system_clock::to_time_t(tp);
        std::tm *ptm = std::localtime(&time_tt);
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(tp.time_since_epoch()) % 1000;
        std::ostringstream oss;
        oss << std::put_time(ptm, "%Y-%m-%dT%H:%M:%S") << "." << std::setw(3) << std::setfill('0') << ms.count();
        return oss.str();
    }
    // Getters for datetime components
    int getYear() const {
        auto time_tt = std::chrono::system_clock::to_time_t(tp);
        std::tm *ptm = std::localtime(&time_tt);
        return ptm->tm_year + 1900;
    }
    int getMonth() const {
        auto time_tt = std::chrono::system_clock::to_time_t(tp);
        std::tm *ptm = std::localtime(&time_tt);
        return ptm->tm_mon + 1;
    }
    int getDay() const {
        auto time_tt = std::chrono::system_clock::to_time_t(tp);
        std::tm *ptm = std::localtime(&time_tt);
        return ptm->tm_mday;
    }
    int getHour() const {
        auto time_tt = std::chrono::system_clock::to_time_t(tp);
        std::tm *ptm = std::localtime(&time_tt);
        return ptm->tm_hour;
    }
    int getMinute() const {
        auto time_tt = std::chrono::system_clock::to_time_t(tp);
        std::tm *ptm = std::localtime(&time_tt);
        return ptm->tm_min;
    }
    int getSecond() const {
        auto time_tt = std::chrono::system_clock::to_time_t(tp);
        std::tm *ptm = std::localtime(&time_tt);
        return ptm->tm_sec;
    }
    int getMillisecond() const {
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(tp.time_since_epoch()) % 1000;
        return static_cast<int>(ms.count());
    }
    // Comparison operator (less than)
    bool operator<(const datetime &other) const {
        return tp < other.tp;
    }
};

class timedelta {
private:
    // Store duration in milliseconds
    std::chrono::milliseconds duration;
public:
    // Constructor from milliseconds value via parameters (e.g., timedelta(0,0,0,0, rtime_m))
    timedelta(int days, int hours, int minutes, int seconds, int millis) {
        long long total_ms = millis;
        total_ms += static_cast<long long>(seconds) * 1000;
        total_ms += static_cast<long long>(minutes) * 60 * 1000;
        total_ms += static_cast<long long>(hours) * 3600 * 1000;
        total_ms += static_cast<long long>(days) * 24 * 3600 * 1000;
        duration = std::chrono::milliseconds(total_ms);
    }
    // Default constructor
    timedelta() : duration(0) {}
    // Getter for total seconds as double
    double total_seconds() const {
        return duration.count() / 1000.0;
    }
    // Allow addition of two timedeltas
    timedelta operator+(const timedelta &other) const {
        timedelta temp;
        temp.duration = duration + other.duration;
        return temp;
    }
    // Allow subtraction of two timedeltas
    timedelta operator-(const timedelta &other) const {
        timedelta temp;
        temp.duration = duration - other.duration;
        return temp;
    }
    // Expose the internal duration for use in datetime arithmetic
    const std::chrono::milliseconds& getDuration() const {
        return duration;
    }
};

// Definition of datetime + timedelta operator
datetime datetime::operator+(const timedelta &td) const {
    datetime temp(*this);
    temp = datetime(); // reset temp to force assignment below
    temp.tp = this->tp + td.getDuration();
    return temp;
}

// Definition of datetime - datetime operator to yield timedelta
timedelta datetime::operator-(const datetime &other) const {
    timedelta td;
    auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(this->tp - other.tp);
    // Directly assign the computed difference in milliseconds
    td = timedelta(0,0,0,0, static_cast<int>(diff.count()));
    return td;
}

//---------------------------------------------------------------------
// Supporting hypfile definitions
//---------------------------------------------------------------------
class cRecordHyp {
public:
    int id;
    int year;
    int month;
    int day;
    int hour;
    int minute;
    int isec;
    int msec;
    std::string sta;
    std::string phase;
    double amp_v;
    double wt;
    double period;
    // Method: getYear()
    int getYear() const {
        return year;
    }
};

class cFileHyp {
public:
    int n_rec;
    bool hyr; // additional flag from Fortran usage: hyp%hyr
    std::vector<cRecordHyp> rec;
    cFileHyp() : n_rec(0), hyr(false) {}
};

//---------------------------------------------------------------------
// Supporting stations module definitions (from stations)
//---------------------------------------------------------------------
namespace stations {
    int nstat = 0;
    std::vector<std::string> stat_name;
    std::vector<double> xstat;
    std::vector<double> ystat;
    std::vector<double> zstat;
}

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
cFileArr arrs;
cFileArr* parrs = &arrs;

// contains

// subroutine cfilearr_todatetime(this, rtime, year, month, day, hour, minute, isec, msec)
void cFileArr::toDateTime(double rtime, int &year, int &month, int &day, int &hour, int &minute, int &isec, int &msec)
{
    // implicit none
    // class(cFileArr)      :: this
    // real, intent(IN)     :: rtime
    // integer, intent(OUT) :: year,month,day,hour,minute,isec,msec
    datetime pcas;
    int rtime_m;
    
    rtime_m = static_cast<int>(std::lround(rtime * 1000));
    std::cout << rtime_m << std::endl;
    pcas = this->reftime + timedelta(0,0,0,0, rtime_m);
    std::cout << pcas.isoformat() << std::endl;
    year = pcas.getYear();
    month = pcas.getMonth();
    day   = pcas.getDay();
    hour  = pcas.getHour();
    minute= pcas.getMinute();
    isec  = pcas.getSecond();
    msec  = pcas.getMillisecond();
} 

// subroutine cfilearr_fill(this, hyp)
void cFileArr::fill(const cFileHyp &hyp)
{
    // use stations,   only: nstat, stat_name, xstat, ystat, zstat
    // implicit none
    // class(cFileArr)      :: this
    // type (cFileHyp)      :: hyp

    int i;
    int j;
    
    double pwt;
    double psum, psum2, mwt;
    datetime pcas, mcas;
    timedelta dcas(0,0,0,0,0);
    cRecordHyp r;
    cRecordArr a;
    // reftime
    mcas = datetime(3000);
    for(i = 0; i < hyp.n_rec; i++) {
        r = hyp.rec[i];
        pcas = datetime(r.getYear(), r.month, r.day, r.hour, r.minute, r.isec, r.msec);
        if(pcas < mcas) {
            mcas = pcas;
        }
    }
    
    this->reftime = datetime(mcas.getYear(), mcas.getMonth(), mcas.getDay(), mcas.getHour(), mcas.getMinute(), 0, 0);
    std::cout << this->reftime.isoformat() << std::endl;
    
    // if (allocated(this%arr)) deallocate(this%arr)
    // allocate(this%arr(hyp%n_rec))
    this->arr.clear();
    this->arr.resize(hyp.n_rec);
    
    this->n_arr = 0;
    for(i = 0; i < hyp.n_rec; i++) {
        r = hyp.rec[i];
        //     if (r%wt .eq. 4) cycle
        a.id = r.id;
        a.sta = r.sta;
        a.phase = r.phase;
        a.amp = r.amp_v;
        if(r.period > 0) {
            a.freq = 1.0 / r.period;
        } else {
            a.freq = 99.9;
        }
        
        // compute weights
        pwt = r.wt;
        
        if (hyp.hyr) {
            if (pwt > 0) {
                pwt = pwt / 1000.0;
                pwt = 1.0 / pwt;
            } else if (pwt < 0) {
                pwt = 0.0;
            } else {
                pwt = 0.0;
            }
        } else {
            pwt = (4.0 - pwt) / 4.0;
        }
        a.wt = pwt;
        
        pcas = datetime(r.getYear(), r.month, r.day, r.hour, r.minute, r.isec, r.msec);
        dcas = pcas - this->reftime;
        a.trec = dcas.total_seconds();
        // key
        a.key = 0;
        for(j = 0; j < stations::nstat; j++) {
            if(a.sta == stations::stat_name[j]) {
                a.key = j + 1;
                a.X = stations::xstat[j];
                a.Y = stations::ystat[j];
                a.Z = stations::zstat[j];
            }
        }
        if(a.key == 0) {
            // write (*,'(1x,"Station ",a5," not found.")') this%arr(i)%sta
            std::cout << " Station " << a.sta << " not found." << std::endl;
            continue;
        }
        this->arr[this->n_arr] = a;
        this->n_arr = this->n_arr + 1;
    }   // i loop over n_rec
    
    //   sum weights
    psum = 0.0;
    psum2 = 0.0;
    mwt = 0.0;
    for(i = 0; i < this->n_arr; i++) {
        pwt = this->arr[i].wt;
        psum = psum + pwt;
        psum2 = psum2 + pwt * pwt;
        if(pwt > mwt)
            mwt = pwt;
    }
    this->sumwt = psum;
    this->sumwt2 = psum2;
    this->maxwt = mwt;
    
    return;
} 

// End of module Arrivals
// (No further code in module)

int main() {
    // For demonstration purposes only, we initialize sample data for hyp and stations.
    
    // Initialize stations data
    stations::nstat = 2;
    stations::stat_name.push_back("STA1");
    stations::stat_name.push_back("STA2");
    stations::xstat.push_back(10.0);
    stations::xstat.push_back(20.0);
    stations::ystat.push_back(100.0);
    stations::ystat.push_back(200.0);
    stations::zstat.push_back(1000.0);
    stations::zstat.push_back(2000.0);
    
    // Create a cFileHyp object and populate with sample records
    cFileHyp hyp;
    hyp.n_rec = 3;
    hyp.hyr = true;
    cRecordHyp rec1;
    rec1.id = 1;
    rec1.year = 2020;
    rec1.month = 5;
    rec1.day = 15;
    rec1.hour = 12;
    rec1.minute = 30;
    rec1.isec = 5;
    rec1.msec = 500;
    rec1.sta = "STA1";
    rec1.phase = "P";
    rec1.amp_v = 3.5;
    rec1.wt = 2000;
    rec1.period = 2.0;
    
    cRecordHyp rec2;
    rec2.id = 2;
    rec2.year = 2020;
    rec2.month = 5;
    rec2.day = 15;
    rec2.hour = 12;
    rec2.minute = 31;
    rec2.isec = 10;
    rec2.msec = 250;
    rec2.sta = "STA2";
    rec2.phase = "S";
    rec2.amp_v = 4.5;
    rec2.wt = 1500;
    rec2.period = 1.5;
    
    cRecordHyp rec3;
    rec3.id = 3;
    rec3.year = 2020;
    rec3.month = 5;
    rec3.day = 15;
    rec3.hour = 12;
    rec3.minute = 32;
    rec3.isec = 20;
    rec3.msec = 750;
    rec3.sta = "STA3"; // This station is not in the stations module list.
    rec3.phase = "P";
    rec3.amp_v = 2.5;
    rec3.wt = 2500;
    rec3.period = 0.0;
    
    hyp.rec.push_back(rec1);
    hyp.rec.push_back(rec2);
    hyp.rec.push_back(rec3);
    
    // Call fill on global arrs
    arrs.fill(hyp);
    
    // Demonstrate toDateTime conversion for a sample rtime (e.g., 1.234 seconds)
    int year, month, day, hour, minute, isec, msec;
    arrs.toDateTime(1.234, year, month, day, hour, minute, isec, msec);
    
    // Output the computed weight sums
    std::cout << "n_arr = " << arrs.n_arr << std::endl;
    std::cout << "sumwt = " << arrs.sumwt << std::endl;
    std::cout << "sumwt2 = " << arrs.sumwt2 << std::endl;
    std::cout << "maxwt = " << arrs.maxwt << std::endl;
    
    return 0;
}
