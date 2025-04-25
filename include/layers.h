#ifndef LAYERS_H
#define LAYERS_H

struct TLayer {
    float z;
    float v;
};

//extern std::map<std::string, std::vector<TLayer>*> layers_pmap;
//std::vector<TLayer>* getLayers(const std::string& wave_name);
extern std::map<std::string, std::shared_ptr<std::vector<TLayer>>> layers_pmap;
std::shared_ptr<std::vector<TLayer>> getLayersVector(const std::string& wave_name);

struct TLayers{
    int n_layers;
    std::vector<float> z;
    std::vector<float> v;
};

void loadLayersPointerMap(const std::string& filename);
void printLayersPointerMap();

extern std::map<std::string, std::shared_ptr<TLayers>> p_layers_map;
std::shared_ptr<TLayers> getLayersPointer(const std::string& ph);
#endif // LAYERS_H
