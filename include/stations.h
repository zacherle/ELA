#ifndef STATIONS_H
#define STATIONS_H
#include <string>

struct TStation {
	std::string name;
	double Xk;
	double Yk;
	float elevation;
	float delayP = 0.0; // default value
	float delayS = 0.0; // default value
};

void loadStations(const std::string& filename);

extern std::map<std::string, TStation> stations;
#endif // STATIONS_H
