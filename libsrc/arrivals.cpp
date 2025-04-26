#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <map>

#include "hypfile.h"
//#include "datetime.hpp"
#include "stations.h"
/*
struct TArrival {
    int id;
    std::string sta;
    std::string phase;
    double trec;
    double wt;
    double amp;
    double freq;
    double X;
    double Y;
    double Z;
    double delayP;
    double delayS;
};

class CArrivalFile {
public:
    int n_arr = 0;
    double reftime;
    double sumwt = 0.0;
    double sumwt2 = 0.0;
    double maxwt = 0.0;
    std::vector<TArrival> arr;
    
    void fill(const std::vector<TRecordHyp> &hyp, const std::map<std::string,TStation> &stations);
};
*/

#include "arrivals.h"


// function for converting time and date to epoch time
// input integer numbers of: year, month, day, hour, minute, second, millisecond
// output double number of seconds since epoch
// epoch time is the number of seconds since 1970-01-01 00:00:00 UTC
// both input and output are in UTC
double toEpochTime(int year, int month, int day, int hour, int minute, int second, int millisecond) {
    std::tm tm = {};
    tm.tm_year = year - 1900; // tm_year is years since 1900
    tm.tm_mon = month - 1;    // tm_mon is months since January (0-11)
    tm.tm_mday = day;
    tm.tm_hour = hour;
    tm.tm_min = minute;
    tm.tm_sec = second;
    tm.tm_isdst=-1;
    tm.tm_zone = NULL;
    tm.tm_gmtoff = 0;
    
    // Convert to time_t (seconds since epoch)
    //std::time_t time_since_epoch = std::mktime(&tm);
    //std::time_t time_since_epoch = std::timelocal(&tm);
    std::time_t time_since_epoch = timelocal(&tm);
    
    // Add milliseconds
    double epoch_time = static_cast<double>(time_since_epoch) + static_cast<double>(millisecond) / 1000.0;
    
    return epoch_time;
}


CArrivalFile arrs;
CArrivalFile* parrs = &arrs;

void CArrivalFile::fill(const std::vector<TRecordHyp> &hyp, const std::map<std::string,TStation> &stations){

// fill std::vector<TArrival> arr with records computed from all records in hyp joinded with stations on station name
// hyp is a vector of TRecordHyp
// stations is a map of station name to TStation
// arr is a vector of TArrival
// fill only arr.sta, arr.phase, arr.freq, arr.X, arr.Y, arr.Z, arr.delayP, arr.delayS
// arr.freq is computed from hyp.period
// arr.X, arr.Y, arr.Z, arr.delayP, arr.delayS are taken from stations
// arr.wt is computed from hyp.wt

// test hyp.wt content and set flag hyr
// if one of hyp.wt are not in enumerated set of integer numbers {0,1,2,3,4}
// then set hyr flag to true else set it to false
   bool hyr = false; // initialize hyr to false
   for (const auto& r : hyp) {
      if (r.wt < 0 || r.wt > 4) {
	  // If any weight is outside the range [0, 4], set hyr to true
	  hyr = true;
	  break;
      }
      // check if r.wt is an integer
      if (r.wt != static_cast<int>(r.wt)) {
          // if r.wt is not integer number then set hyr to true
          hyr = true;
          break;
      }
   }

   n_arr = 0; // reset number of arrivals
	     
   for (const auto& r : hyp) {
	TArrival a;
	a.id = n_arr+1; // set id to current number of arrivals + 1
	a.sta = r.sta;
	a.phase = r.phase;
	a.trec = 0.0; // initialize to zero
	a.wt = 0.0; // initialize to zero
	a.amp = r.amp_v;

        // convert r.year from YY format to YYYY format
	// r.year is in YY format
	// int yyyy to the year in YYYY format
	int yyyy;
        if (r.year < 70) {
	    yyyy = r.year + 2000; // convert to YYYY format
	} else {
	    yyyy = r.year + 1900; // convert to YYYY format
	}
	// compute epoch time a.trec from fields of TRecordHyp
 
	a.trec = toEpochTime(yyyy, r.month, r.day, r.hour, r.minute, r.isec, r.msec);

	// find station in stations map
	auto it = stations.find(r.sta);
	if (it != stations.end()) {
	    const TStation& station = it->second;
	    a.X = station.Xk;
	    a.Y = station.Yk;
	    a.Z = station.elevation;
	    a.delayP = station.delayP;
	    a.delayS = station.delayS;
	} else {
	    std::cerr << "Station " << r.sta << " not found in stations map." << std::endl;
	    continue; // skip this record if station not found
	}
	
	// compute frequency
	if (r.period > 0) {
	    a.freq = 1.0 / r.period;
	} else {
	    a.freq = 99.9; // default value for frequency
	}
	
	// compute weights
	double pwt = r.wt;
	if (hyr) {
	// hyp.wt is in miliseconds
	    if (pwt > 0) {
		pwt = pwt / 1000.0;
		pwt = 1.0 / pwt;
	    } else if (pwt < 0) {
		pwt = 0.0; // negative weight set to zero
	    } else {
		pwt = 0.0; // zero weight remains zero
	    }
	} else {
	// hyp.wt is acording HYPO71
	    pwt = (4.0 - pwt) / 4.0; // normalize weight
	}
	
	a.wt = pwt; // assign computed weight to arrival
	
	arr.push_back(a); // add arrival to vector
	n_arr++; // increment number of arrivals
   } // end of loop over hyp records
   // now fill the other members of the class
   // reftime is the minimal time of the first arrival rounded to the floor to whole tents
   reftime = arr[0].trec; // initialize
   for (int i = 0; i < n_arr; i++) {
      if (arr[i].trec < reftime) {
      reftime = arr[i].trec;
      }
   }
   // reftime round to the floor to whole tens
   reftime = floor(reftime/10.0)*10.0; // round to the floor of the second

   // sum weights
   sumwt = 0.0;
   sumwt2 = 0.0;
   maxwt = 0.0;
   for(int i = 0; i < n_arr; i++) {
	// arr[i].trec = arr[i].trec - reftime; // remove reference time
      double pwt = arr[i].wt;
      sumwt = sumwt + pwt;
      sumwt2 = sumwt2 + pwt * pwt;
      if(pwt > maxwt) maxwt = pwt;
    }

   for(int i = 0; i < n_arr; i++) {
      arr[i].trec = arr[i].trec - reftime; // remove reference time
     // arr[i].wt = arr[i].wt / sumwt; // normalize weights
   }
    // std::cout << "Number of arrivals: " << n_arr << " / "  << arr.size() << std::endl;
    // std::cout << "reftime = " << reftime << std::endl;
    // std::cout << "sumwt = " << sumwt << std::endl;
    // std::cout << "sumwt2 = " << sumwt2 << std::endl;
    // std::cout << "maxwt = " << maxwt << std::endl;
}

/*    
    // Demonstrate toDateTime conversion for a sample rtime (e.g., 1.234 seconds)
    int year, month, day, hour, minute, isec, msec;
    arrs.toDateTime(1.234, year, month, day, hour, minute, isec, msec);
    
    // Output the computed weight sums
    std::cout << "n_arr = " << arrs.n_arr << std::endl;
    std::cout << "sumwt = " << arrs.sumwt << std::endl;
    std::cout << "sumwt2 = " << arrs.sumwt2 << std::endl;
    std::cout << "maxwt = " << arrs.maxwt << std::endl;
    
    return 0;
*/
