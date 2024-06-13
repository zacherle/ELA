// https://stackoverflow.com/questions/45466939/parsing-command-line-arguments-with-boost-program-options-c
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <iostream>
#include <fstream>
#include <iterator>
#include <exception>

//using namespace std;
using std::cerr;
using std::cout;
using std::endl;
using std::exception;
//using std::vector;
using std::string;

#include "../include/c_f_elalib.h"

#ifdef __cplusplus
extern"C" {
#endif
void ela_d_(void);
#ifdef __cplusplus
}
#endif


int main(int ac, char** av){

string file_i;
string file_o;
string file_m;
string file_s;
double reading_err;
double model_err;
string config_file;

bool option_missing(false);

try {

    po::options_description desc("Usage: ela [options] input \n" \
		    "Allows options:");
    
    po::options_description generic("Generic options");
    generic.add_options()
	("version,v", "print version string")
        ("help,h", "produce help message")
	("config,c",
	 po::value<string>(&config_file)->default_value("ela.cfg"),
	 "name of a file of a configuration.");
	
    po::options_description cmd("Command line I/O options");
    cmd.add_options()
        ("input,i", po::value<string>(), "input hyp file")
        ("output,o", po::value<string>(), "output hy3 file");

    po::options_description config("Configuration");
    config.add_options()
        ("model,m", po::value<string>(), "input model file")
        ("stations,s", po::value<string>(), "input list of stations")
        ("reading_err", po::value<double>(&reading_err)->default_value(0),
	 "set error of arrival time observationd [s]")
        ("model_err", po::value<double>(&model_err)->default_value(0),
	 "set error of propagation times predictions [s]");

    desc.add(generic).add(cmd).add(config);

    po::options_description config_file_options;
    config_file_options.add(cmd).add(config);

    po::variables_map vm;
//    po::store(po::parse_command_line(ac, av, desc), vm);


    po::positional_options_description p;
    p.add("input", -1);

    po::store(po::command_line_parser(ac, av).
          options(desc).positional(p).run(), vm);
    po::notify(vm);


    if (vm.count("help")) {
        cout << desc << "\n";
        return 0;
      }

    if (vm.count("version")) {
        cout << "ELA v00.1 2024-05" << "\n";
        return 0;
      }

    std::ifstream ifs(config_file.c_str());
    if (!ifs) {
        cout << "can not open config file: " << config_file << "\n";
        return 0;
      }
    else {
	po::store(po::parse_config_file(ifs, config_file_options), vm);
	po::notify(vm);
      }

    if (vm.count("input")) {
	file_i = vm["input"].as<string>();
        cout << "Input file: " << file_i << "\n";
        }
    else {
        cout << "Input file was not set (option -i).\n";
	option_missing=true;
        }

    if (vm.count("output")) {
        file_o = vm["output"].as<string>();
        cout << "Output file: " << file_o << "\n";
        }
    else if (vm.count("input")) {
        file_o = vm["input"].as<string>();
	file_o = file_o.substr(file_o.find_last_of("/")+1,file_o.npos); 
	file_o = "./ela_"+file_o.substr(0, file_o.find_last_of("."))+".hy3"; 
        cout << "Output file: " << file_o << "\n";
        }
    else {
        cout << "Output file was not set (option -o).\n";
	option_missing=true;
        }

    if (vm.count("model")) {
	file_m = vm["model"].as<string>();
        cout << "Velocity model file: " << file_m << "\n";
        }
    else {
        cout << "Velocity model file was not set (option -m).\n";
	option_missing=true;
        }

    if (vm.count("stations")) {
        file_s = vm["stations"].as<string>();
        cout << "Stations list file: " << file_s << "\n";
        }
    else {
        cout << "Stations list was not set (option -s).\n";
	option_missing=true;
        }

    if (vm.count("reading_err")) {
        cout << "Reading error: "
             << vm["reading_err"].as<double>() << "\n";
        }

    if (vm.count("model_err")) {
        cout << "Model error: "
             << vm["model_err"].as<double>() << "\n";
        }

    } //try

    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
    return 1;
    }

    if (option_missing) {return 1;}

string name_model=file_m.substr(file_m.find_last_of("/")+1,file_m.npos); 
setpar_model(model_err,reading_err, (char*)name_model.c_str());
// 1)   read_stations     and fill list of stations
filesta_read((char*)file_s.c_str());
// 2)   read_hyp          and fill cFileHyp::hyp, cFileArr::arrs
filehyp_read((char*)file_i.c_str());
// 3)   read_model        and fill cLayers::velmod(5)
filemod_read((char*)file_m.c_str());

setpar_name_o((char*)file_o.c_str());
ela_d_();

return 0;
}
