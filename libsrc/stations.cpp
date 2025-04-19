#include<iostream> 
#include<map> // for map operations 
#include<string>
#include<fstream>
//#include<cstdio>
#include <sstream>
#include <iomanip>
/*
struct TStation {
	std::string name;
	double Xk;
	double Yk;
	float elevation;
	float delayP = 0.0; // default value
	float delayS = 0.0; // default value
};
*/
#include "stations.h"

std::map<std::string, TStation> stations;

// load stations from file
// the file format is:
// name Xk Yk elevation delayP delayS
// or alternatively:
// name Xk Yk elevation delayP
// the last two values are optional
// for example:
// 'VRAC' 1148.113 597.870 0.470 0.0 0.050
// 'VRAC1' 1148.113 597.870 0.470 
// 'VRAC2' 1148.113 597.870 0.470 0.0

void loadStations(const std::string& filename) {
	std::ifstream file(filename);
	if (!file.is_open()) {
		std::cerr << "Error opening file: " << filename << std::endl;
		return;
	}

	std::string line;
	while (std::getline(file, line)) {
		std::istringstream iss(line);
		TStation station;
		if (!(iss >> station.name >> station.Xk >> station.Yk >> station.elevation)) {
			std::cerr << "Error reading line: " << line << std::endl;
			continue; // Skip to the next line
		}
		if (!(iss >> station.delayP)) {
			station.delayP = 0.0; // default value
		}
		if (!(iss >> station.delayS)) {
			station.delayS = station.delayP; // default value
		}
		// Strip the station name of quotes if present

		if (station.name.front() == '\'' && station.name.back() == '\'') {
			station.name = station.name.substr(1, station.name.size() - 2);
		 	}
		// Check if the station name is already in the map
		if (stations.find(station.name) != stations.end()) {
		 	std::cerr << "Duplicate station name: " << station.name << std::endl;
		 	continue; // Skip to the next line
			}

		// Add the station to the map
		stations[station.name] = station;
	}
	file.close();
}

// print all stations
// the format is:
// the same as in the file
// for example:
// 'VRAC' 1148.113 597.870 0.470 0.0 0.050
// the stations are sorted by name

void printStations() {
   std::map<std::string, TStation>::iterator istations;
   std::cout << std::fixed << std::setprecision(3);
   for (istations = stations.begin(); istations != stations.end(); ++istations) {
	std::cout << "'" << istations->second.name << "' " 
	          << istations->second.Xk << " " 
	          << istations->second.Yk << " " 
	          << istations->second.elevation << " " 
	          << istations->second.delayP << " " 
	          << istations->second.delayS << std::endl;
   }
}

// main() for testing
// load stations from file and print them
// the file name is passed as a command line argument
// for example:
// ./stations stations.txt

int main(int argc, char* argv[]) {
	if (argc != 2) {
		std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
		return 1;
	}

	loadStations(argv[1]);
	printStations();

	return 0;
}

