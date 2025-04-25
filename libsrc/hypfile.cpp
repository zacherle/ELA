#include<iostream> 
#include<vector> 
#include<string>
#include<fstream>
//#include<cstdio>
#include <sstream>
#include <iomanip>

/*
struct TRecordHyp {
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
};
*/

#include "hypfile.h"

std::vector<TRecordHyp> hyp;

// load hyp from file
// the file format is:
// id sta phase ichan year month day hour minute isec msec micros wt amp_d pol amp_v period
// for example:
// VRAC  S 99 25 02 11 10 03 22 496 000 0 6.030000e-09  0 7.041980e-07 5.380249e-02

void loadHyp(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
	std::cerr << "Error opening file: " << filename << std::endl;
	return;
    }

    std::string line;
    while (std::getline(file, line)) {
	std::istringstream iss(line);
	TRecordHyp record;
	if (!(iss >> record.sta >> record.phase >> record.ichan
		>> record.year >> record.month >> record.day
	       	>> record.hour >> record.minute >> record.isec >> record.msec >> record.micros
	       	>> record.wt >> record.amp_d >> record.pol >> record.amp_v >> record.period)) {
	    std::cerr << "Error reading line: " << line << std::endl;
	    continue; // Skip to the next line
	}
	hyp.push_back(record);
    }
}

// Print hyp records
// the file format is:
// sta phase ichan year month day hour minute isec msec micros wt amp_d pol amp_v period
// for example:
// VRAC  S 99 25 02 11 10 03 22 496 000 0 6.030000e-09  0 7.041980e-07 5.380249e-02

void printHyp() {
    std::vector<TRecordHyp>::iterator ihyp;
    for (ihyp = hyp.begin(); ihyp != hyp.end(); ++ihyp) {
	std::cout	<< std::setfill(' ') << std::left
		        << std::setw(6) 	<< ihyp->sta
			<< std::setw(5) 	<< ihyp->phase;
        std::cout       << std::setfill('0') << std::right
			<< " "	<< std::setw(2)	<< ihyp->ichan
		        << " "	<< std::setw(2)	<< ihyp->year
			<< " " 	<< std::setw(2)	<< ihyp->month
			<< " "	<< std::setw(2)	<< ihyp->day
			<< " "	<< std::setw(2) << ihyp->hour
			<< " "	<< std::setw(2)	<< ihyp->minute
			<< " "	<< std::setw(2)	<< ihyp->isec
			<< " "	<< std::setw(3)	<< ihyp->msec
			<< " "	<< std::setw(3)	<< ihyp->micros
			<< std::setfill(' ')
			<< std::defaultfloat
			<< std::setw(6)         << ihyp->wt
			<< std::scientific
			<< std::setprecision (5)
			<< std::setw(12)        << ihyp->amp_d
			<< std::setw(3)         << ihyp->pol
			<< std::setw(12)        << ihyp->amp_v
			<< std::setw(12)        << ihyp->period
			<< std::endl;
    }
}

// main() for testing
// load hyp records from file and print them
/*
int main(int argc, char* argv[]) {
	if (argc != 2) {
		std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
		return 1;
	}

	loadHyp(argv[1]);
	printHyp();

	return 0;
}
*/
