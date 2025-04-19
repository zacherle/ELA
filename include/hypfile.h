#ifndef HYPFILE_H
#define HYPFILE_H

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

extern std::vector<TRecordHyp> hyp;

#endif // HYPFILE_H
