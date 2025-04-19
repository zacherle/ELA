#ifndef _HY3FILE_H
#define _HY3FILE_H

#include <time.h>

#ifdef __cplusplus
extern "C" {
#endif

struct hy3_record
  {
              char sta[10];  // station name
              char ph[10];   // phase name
	      float obs_t;   // observed time relative to the reference time
	      float cal_t;   // calculated time relative to the reference time
	      float res;     // residual time
	      float amp;     // amplitude displacement
	      float freq;    // frequency
	      float w;       // weight
	      float epi;     // epicentral distance
	      float hypo;    // hypocentral distance
	      float azm;     // azimuth
	      float ain;     // take-off angle
	      float xmag;    // magnitude
}; // end of struct hy3_record

struct hy3_file
{
//model         :testdata/ete_3d_a.mod
	      char model[40];
//model error   :0.000 s
	      float model_error;
//reading error :0.016 s
	      float reading_error;
//create time   :17-04-25 20:06:31
//	      struct tm create_time;
//event         :testdata/testP.hyp
	      char event[40];
//start(x,y,z,t):(   0.00,  15.00, 999.00,  21.00)
	      float start_x, start_y, start_z, start_t;
//fixed coordinates:(       ,       , fix Z ,       )
              int fix[4];
//reference time:17-04-25 08:59
	      struct tm ref_time;
//origin time          t:  17-04-25  08:59:09.305 +- 23.270
	      //struct date_time origin_time;
	      float origin_time[2];  // to the reference time
//x-coordinate         x:  1309.90 +- 155.58    km     (fi: 47.664119 deg)
	      float x[2];
	      float lon;
//y-coordinate         y:   780.35 +- 138.51    km (lambda: 14.410148 deg)
	      float y[2];
	      float lat;
//depth                z:    -1.26 +-   0.00    km
	      float z[2];
//magnitude           ml:    -9.90 +-   0.00
	      float ml[2];
//rms of time residuals :          2.54         s
	      float rms;
//angular gap           :           320         deg
	      float gap;
//info of iterations  :            11
	      int info;
//error ellipse axis l1 :        202.51         km
	      float eel1;
//              axis l2 :         48.78         km
	      float eel2;
//              theta   :     41.3 deg (to grid)   (azimuth: 213.4 deg)
	      float theta;
	      float azim;
//
              int nrec;
              struct hy3_record * rec;
}; // end of struct hy3_file


void hy3print(struct hy3_file * hy3);
int hy3load (struct hy3_file *hy3, FILE *fin);

#ifdef __cplusplus
}
#endif
#endif // _HY3FILE_H
