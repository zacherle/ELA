#ifndef ARRIVALS_H
#define ARRIVALS_H

#include <vector>
#include <map>
#include "hypfile.h"
#include "stations.h"

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

extern CArrivalFile arrs;
extern CArrivalFile* parrs;

#endif // ARRIVALS_H 
