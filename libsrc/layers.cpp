#include<iostream> 
#include<vector> 
#include<string>
#include<fstream>
//#include<cstdio>
#include <sstream>
#include <iomanip>
#include<map>
#include<memory>
/*
struct TLayer {
    float z;
    float v;
};

struct TLayers{
    int n_layers;
    std::vector<float> z;
    std::vector<float> v;
};
*/


#include "layers.h"

std::map<std::string, std::vector<TLayer>> layers_map;
//std::map<std::string, std::vector<TLayer>*> layers_pmap;
// map like layers_pmap but using std::shared_ptr instead of raw pointer
//std::map<std::string, std::shared_ptr<std::vector<TLayer>>> layers_pmap;
std::map<std::string, std::shared_ptr<TLayers>> p_layers_map;

std::shared_ptr<TLayers> get_layers(const std::string& ph);




std::vector<TLayer> layers;
//std::vector<TLayer>::iterator ilayer;

// load layers from file
// the file contains more layers for different seismic waves
// the first line contains the keyword 'wave' and name of the wave
// the next lines contain the layers
// the first column is the depth (z) and the second column is the velocity (v)
// the depth is in km and the velocity is in km/s
// next 'wave` keyword starts an other layers
// the file format is for example:
// wave P
// 0.0 5.400
// 0.1 5.600
// 0.3 5.700
// 2.4 5.800
// 6.4 6.000
// 15.4 6.500
// 21.4 6.800
// 25.4 7.000
// 35.4 8.140
// wave S
// 0.0 3.214
// 0.1 3.333
// 0.3 3.393
// 2.4 3.452
// 6.4 3.571
// 15.4 3.869
// 21.4 4.048
// 25.4 4.167
// 35.4 4.845

void loadLayersMap(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
	std::cerr << "Error opening file: " << filename << std::endl;
	return;
    }

    std::string line;
    std::string current_wave;
    while (std::getline(file, line)) {
	if (line.find("wave") == 0) {
	    current_wave = line.substr(5); // Extract wave name
	    layers_map[current_wave] = std::vector<TLayer>(); // Initialize vector for this wave
	} else {
	    std::istringstream iss(line);
	    TLayer record;
	    if (!(iss >> record.z >> record.v)) {
		std::cerr << "Error reading line: " << line << std::endl;
		continue; // Skip to the next line
	    }
	    layers_map[current_wave].push_back(record);
	}
    }
}


// Alternative function to load p_layers_map from file
// similar to loadLayersMap but stores a pointer
// to the struct TLayers

void loadLayersPointerMap(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
	std::cerr << "Error opening file: " << filename << std::endl;
	return;
    }

    std::string line;
    std::string wave_names="";
    std::shared_ptr<TLayers> p_layers = nullptr; // Shared pointer to struct
    while (std::getline(file, line)) {
	// Skip empty lines
	// Skip lines that start with '#'
	// Skip lines that start with '!'
        if (line.empty() || line[0] == '#' || line[0] == '!') {
	    continue; // Skip to the next line
	}	
	// Skip lines that start with 'wave' and contain only whitespace
	if (line.find("wave") == 0 && line.find_first_not_of(" \t", 4) == std::string::npos) {
	    p_layers = nullptr;
	    continue; // Skip to the next line
	}
	// Check if the line starts with "wave"
	if (line.find("wave") == 0) {
	    wave_names = line.substr(5); // Extract list of wave names
	    // Split the wave_names string acording to separators (whitespace)
	    // and store parts (seismic phases) in a vector
	    std::vector<std::string> phases;
	    std::istringstream iss(wave_names);
	    std::string phase;
	    while (iss >> phase) {
		    phases.push_back(phase);
	    }
            if (phases.size() > 0) {
		// Create a new struct and vectors for the layers 
	        p_layers = std::make_shared<TLayers>(); // Initialize struct for this wave 
	        p_layers->n_layers=0;
	    }
	    else { 
		std::cerr << "Error: no wave name found" << std::endl;
		p_layers = nullptr;
		continue; // Skip to the next line
	    }
	    for (const auto& phase : phases) {
		    p_layers_map[phase] = p_layers; // Store pointer to vector to the map
	    }
	} else {
	    // If the line does not start with "wave", read the layer data
	    std::istringstream iss(line);
	    TLayer record;
	    if (!(iss >> record.z >> record.v)) {
		std::cerr << "Error reading line: " << line << std::endl;
		continue; // Skip to the next line
	    }
	    // Store record to the pointed struct
	    if (p_layers) { 
	    	p_layers->z.push_back(record.z);
	    	p_layers->v.push_back(record.v);
	    	p_layers->n_layers++;
	    }
	}
    }
}


// Print records of layers_map to standard output
// the file format is the same as the input file:
// wave P
// 0.0 5.400
// 0.1 5.600
// 0.3 5.700
// 2.4 5.800
// ...
// wave S  	
// 0.0 3.214
// 0.1 3.333
// 0.3 3.393
// 2.4 3.452
// ...
//
// the function prints the layers_map to standard output
// the function is used for testing
// the function is used in the main() function
void printLayersMap() {
    for (const auto& pair : layers_map) {
	std::cout << "wave " << pair.first << std::endl;
	for (const auto& layer : pair.second) {
	    std::cout << std::fixed << std::setprecision(1) << layer.z << " " 
		      << std::fixed << std::setprecision(3) << layer.v << std::endl;
	}
    }
}

// Print records of p_layers_map to standard output
// the file format is the same as the input file:
// similar to printLayersMap

void printLayersPointerMap() {
    for (const auto& pair : p_layers_map) {
	std::cout << "wave " << pair.first << std::endl;
	const TLayers& layers= *(pair.second);
	for (size_t i = 0; i < layers.z.size(); ++i) {
	    std::cout << std::fixed << std::setprecision(1) << layers.z[i] << " " 
		      << std::fixed << std::setprecision(3) << layers.v[i] << std::endl;
	}
    }
}
/*
// Get the std::vector<TLayer> for the given wave name
// the function returns a pointer to the vector
std::vector<TLayer>* getLayers(const std::string& wave_name) {
    auto it = layers_pmap.find(wave_name);
    if (it != layers_pmap.end()) {
	return it->second; // Return pointer to the vector
    } else {
	std::cerr << "Wave name not found: " << wave_name << std::endl;
	return nullptr; // Return nullptr if not found
    }
}
*/

// Get the std::vector<TLayer> for the given wave name
// the function returns a pointer to the vector
// the function is similar to getLayers
// but uses std::shared_ptr instead of raw pointer
//std::shared_ptr<std::vector<TLayer>> getLayersPointer(const std::string& wave_name) {
std::shared_ptr<TLayers> getLayersPointer(const std::string& wave_name) {
    auto it = p_layers_map.find(wave_name);
    if (it != p_layers_map.end()) {
	return it->second; // Return pointer to the struct
    } else {
	std::cerr << "Wave name not found: " << wave_name << std::endl;
	return nullptr; // Return nullptr if not found
    }
}

/*
// main() for testing
// load layers_map records from file and print them
int main(int argc, char* argv[]) {
    if (argc != 2) {
	std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
	return 1;
    }

    //loadLayersMap(argv[1]);
    //printLayersMap();

    loadLayersPointerMap(argv[1]);
    printLayersPointerMap();


    // select data for wave name Pn
    std::string wave_name = "Pn";
    std::shared_ptr<TLayers> layers = getLayersPointer(wave_name);
    std::cout << "wave " << wave_name << std::endl;
    if (layers) {
	for (size_t i = 0; i < layers->z.size(); ++i) {
	    std::cout << std::fixed << std::setprecision(1) << layers->z[i] << " " 
		      << std::fixed << std::setprecision(3) << layers->v[i] << std::endl;
       	}
    }
    else {
      std::cerr << "No data found for wave name: " << wave_name << std::endl;
    }
   
    // Clean up dynamically allocated memory im p_layers_map
    // delete the vector pointed by the shared_ptr
    for (auto& pair : p_layers_map) {
	// No need to delete the shared_ptr, it will be automatically deleted when it goes out of scope
	// However, if you want to clear the vector, you can do it like this:
	pair.second->z.clear(); // Clear the vector
	pair.second->v.clear(); // Clear the vector
    }
   
    p_layers_map.clear(); // Clear the map

    return 0;
}
*/
