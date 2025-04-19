#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <ctime>
#include <algorithm>

#include "version.h"
//#include "layers.h"
#include "arrivals.h"
#include "hypfile.h"
#include "const_raddeg.h"
#include "gather.h"
#include "wgsjtsk.h"
#include "hy3file.h"
#include "magni.h"


extern "C" void get_errellipse_(const double* co,
	       	float dxer, float dyer, float dzer, float dter,
	       	float l1, float l2, float theta);

double mconvergence(double X, double Y)
{
	// Calculate meridian convergence
	// X,Y in kilometers
	return 0.008257 * Y + 2.373 * Y / X;
}

float maxgap(std::vector<float>& az)
{
    std::sort(az.begin(), az.end());
    float max_gap = 0.0f;
    for (int i = 1; i < az.size(); i++) {
         float gap = az[i] - az[i-1];
         if(gap > max_gap) max_gap = gap;
    }
    float circular_gap = 360.0f - az.back() + az.front();
    if(circular_gap > max_gap) max_gap = circular_gap;
    return max_gap;
}

void output_hy3(struct hy3_file hy3, int marr,
	       	const double hypo[4], const float startpoint[4],
		const CArrivalFile &arrs, std::vector<TRecordHyp> &hyp, Gather &gather,
		const std::string &modfn, const double* co, int info)
{ 
/// Output the results of the inversion to the hy3file struct


    modfn.copy(hy3.model, 40);
    //hy3.model_error = model_error;
    //hy3.reading_error = reading_error;
    //hy3.event =
    hy3.start_x = startpoint[0];
    hy3.start_y = startpoint[1];
    hy3.start_z = startpoint[2];
    hy3.start_t = startpoint[3];
    //hy3.fix[0] = 0;


    // origin time epoch
    double origin_et=arrs.reftime+hypo[3];
    // the reference time is the origin time round to the nearst minute
    // reference time epoch
    double ref_et = std::floor(origin_et/60.0) * 60.0;
    // set the origin time (relative to the reference time)
    hy3.origin_time[0] = static_cast<float>(origin_et-ref_et);

    // convert ref_time to time_t
    time_t ref_time_t = static_cast<time_t>(ref_et);
    // convert to tm struct
    struct tm* ref_tm;
    ref_tm = std::localtime(&ref_time_t);
    // set the reference time
    hy3.ref_time.tm_year = ref_tm->tm_year;
    hy3.ref_time.tm_mon = ref_tm->tm_mon;
    hy3.ref_time.tm_mday = ref_tm->tm_mday;
    hy3.ref_time.tm_hour = ref_tm->tm_hour;
    hy3.ref_time.tm_min = ref_tm->tm_min;
    hy3.ref_time.tm_sec = ref_tm->tm_sec;       //=0
    hy3.ref_time.tm_isdst = ref_tm->tm_isdst;   //=-1
    hy3.ref_time.tm_zone = ref_tm->tm_zone;     //=NULL
    hy3.ref_time.tm_gmtoff = ref_tm->tm_gmtoff; //=0
 
    // gather.collect(hypo);

    // Calculate meridian convergence
    double meridian_con = mconvergence(hypo[0], hypo[1]);
    
    std::vector<float> az(marr, 0.0f); // azimuths for gap calculation
    int j_az = 0;
    int nmag = 0;
    double smag = 0.0;
    double smag2 = 0.0;
    for (int i = 0; i < arrs.n_arr; i++) {
        TArrival a;
        TRecordHyp r;
	a = arrs.arr[i];
	//r = hyp.rec(a.id);
	
	double dx = a.X - hypo[0];
	double dy = a.Y - hypo[1];
	double dz = a.Z - hypo[2];
	// calculate azimuth
	double aaz = std::atan2(dy, dx) * RAD2DEG;
   	aaz = std::fmod(720.0 + aaz - 180.0 - meridian_con, 360.0);
	if(aaz < 0) aaz += 360.0;
        if(a.wt > 0.0) {
	    az[j_az++] = static_cast<float>(aaz);
	}

	// calculate distance
	double dhypo = std::sqrt(dx * dx + dy * dy + dz * dz);
	double depi = std::sqrt(dx * dx + dy * dy);
       
        // station magnitude	
	double amag = localmagnitude(hypo, a);
        if (amag > -9.9) {
            nmag = nmag + 1;
            smag += amag;
            smag2 += amag * amag;
        }

	double aobs = a.trec + arrs.reftime-ref_et;
	//double aobs = gather.trec[i] + arrs.reftime-ref_et;
        double acal = gather.ttcal[i] + origin_et-ref_et;
        double ares = aobs - acal;
        double atoa = gather.toa[i];

	a.sta.copy(hy3.rec[i].sta,10);
	a.phase.copy(hy3.rec[i].ph,10);
	hy3.rec[i].obs_t = aobs;
	hy3.rec[i].cal_t = acal;
	hy3.rec[i].res = ares;
	hy3.rec[i].amp = a.amp;
	hy3.rec[i].freq = a.freq;
	hy3.rec[i].w = a.wt;
	hy3.rec[i].epi = static_cast<float>(depi);
	hy3.rec[i].hypo = static_cast<float>(dhypo);
	hy3.rec[i].azm = static_cast<float>(aaz);
	hy3.rec[i].ain = static_cast<float>(atoa);
	hy3.rec[i].xmag = static_cast<float>(amag);

    } // end do loop

    // azimuthal gap
    hy3.gap = maxgap(az);

    // local magnitude
    hy3.ml[0] = smag/nmag;
    hy3.ml[1] = std::sqrt((smag2-(smag*smag)/nmag)/nmag);
       
    // rms of residuals 
    double sum_num = 0.0;
    double sum_w = 0.0;
    for (int i = 0; i < arrs.n_arr; i++) {
        double res = gather.trec[i] - gather.ttcal[i] - hypo[3];
        double w_res = res * gather.wt[i];
        sum_num += w_res * w_res;
        sum_w += gather.wt[i];
    }
    if(sum_w != 0.0)
        hy3.rms = std::sqrt(sum_num / sum_w);
    else
        hy3.rms = 0.0;

    // error ellipse
    float dxer, dyer, dzer, dter, l1, l2, theta;
    get_errellipse_(co, dxer, dyer, dzer, dter, l1, l2, theta);
    hy3.x[0] = hypo[0];
    hy3.x[1] = dxer;
    hy3.y[0] = hypo[1];
    hy3.y[1] = dyer;
    hy3.z[0] = hypo[2];
    hy3.z[1] = dzer;
    hy3.origin_time[1] = dter;
    hy3.eel1 = l1; // semi-major axis in km
    hy3.eel2 = l2; // semi-minor axis in km

    double az_theta;
    az_theta = theta - 180.0 - meridian_con;
    az_theta = std::fmod(360.0 + az_theta, 360.0);
    if (az_theta < 0) az_theta += 360.0;
    // azimuth in degrees
    hy3.theta = theta;   // to grid
    hy3.azim = az_theta;

    // coordinates Krovak --> geographic
    double fi, rla;
    jtsk2wgs(hypo[0]*1000.0, hypo[1]*1000.0, &fi, &rla);
    hy3.lon = rla;
    hy3.lat = fi;
    hy3.z[0] = hypo[2];
    hy3.z[1] = dzer;
    return;

} // end output_hy3

