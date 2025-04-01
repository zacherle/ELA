#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>

class cLayers {
public:
    std::string phase = " ";
    int n_layers = 0;
    std::vector<double> z;
    std::vector<double> v;

    int layers_load(std::ifstream &iu);
};

cLayers velmod[5];

double model_error = 0.0;
double reading_error = 0.016;

int layers_load(cLayers &lay, std::ifstream &iu) {
    int rmsg = 0;
    int ios;
    std::string iom;
    char cbuf[40];
    int icbuf;
    int i = 0;
    double z_tmp[500];
    double v_tmp[500];
    int iphase = 0;
    std::string phase;

    iphase = 0;
    i = 0;
    while (true) {
        iu.getline(cbuf, 40);
        ios = iu.fail() ? -1 : 0; // Simulating IOSTAT

        if (ios > 0) { // ios > 0 error
            iu.close();
            std::cout << "Encountered Error:\n" << iom << std::endl;
            rmsg = ios;
            return rmsg;
        }
        if (ios < 0) { // ios = -1 end of file
            if (i > 0) {
                lay.n_layers = i;
                lay.z.resize(i);
                lay.v.resize(i);
                for (int j = 0; j < i; ++j) {
                    lay.z[j] = z_tmp[j];
                    lay.v[j] = v_tmp[j];
                }
            }
            rmsg = ios;
            return rmsg;
        }
        if (ios == 0) { // ios == 0 ok
            icbuf = std::string(cbuf).find('#');
            if (icbuf > 0 && icbuf < 5) continue;
            icbuf = std::string(cbuf).find("wave");
            if (icbuf > 0) { // keyword wave
                phase = std::string(cbuf).substr(5);
                iphase++;
                if (iphase == 1) {
                    lay.phase = phase;
                } else if (iphase > 1) {
                    iu.clear(); // Simulating backspace
                    lay.n_layers = i;
                    lay.z.resize(i);
                    lay.v.resize(i);
                    for (int j = 0; j < i; ++j) {
                        lay.z[j] = z_tmp[j];
                        lay.v[j] = v_tmp[j];
                    }
                    break;
                }
            } else {
                i++;
                sscanf(cbuf, "%lf %lf", &z_tmp[i - 1], &v_tmp[i - 1]);
            }
        }
    }
    rmsg = ios;
    return rmsg;
}

cLayers* get_layers(const std::string &ph) {
    cLayers* lptr = nullptr;
    for (int i = 0; i < 5; ++i) {
        if (velmod[i].n_layers > 0) {
            if (velmod[i].phase.find(ph) == 0) {
                lptr = &velmod[i];
                return lptr;
            }
        }
    }
    return lptr;
}
