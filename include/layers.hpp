class cLayers {
public:
    std::string phase;
    int n_layers;
    std::vector<double> z;
    std::vector<double> v;
    cLayers() {}
    int layers_load(std::ifstream &iu);
};

extern cLayers velmod[5];

extern double model_error = 0.0;
extern double reading_error = 0.016;

// Namespace layers (from module layers)
namespace layers {
    const double model_error = 0.123;
    const double reading_error = 0.456;
    const std::string name_model = "DefaultModel";
}
int layers_load(cLayers &lay, std::ifstream &iu);

cLayers* get_layers(const std::string &ph);
