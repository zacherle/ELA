#ifndef PARAM_H
#define PARAM_H

using std::string;
#include <array>

struct TParams {
    double reading_err;
    double model_err;
    string name_model;
    string name_event;
    std::array<bool,4> fix;
    bool init_nea;
    std::array<double,4> startpt;
};

#endif // PARAM_H
