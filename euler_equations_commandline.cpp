#include <iostream>
#include <cmath>
#include <complex>
#include <array>
#include <vector>
#include <functional>
#include <algorithm>
#include <iomanip>
#include <getopt.h>
#include <fstream>
#include <map>

double ncells = 1000;
double maxtime = 0.15;
double densityL = 1;
double densityR = 1;
double velocityL = -2;
double velocityR = 2;
double pressureL = 0.4;
double pressureR = 0.4;
std::string outPath = "-";

void PrintHelp()
{
    std::cout <<
        "Parameter settings:\n"
        "   -n, --steps           Set number of steps across domain\n"
        "   -t, --time            Set time to simulate until\n"
        "   -d, --density_left    Set density on left side of domain\n"
        "   -D, --density_right   Set density on right side of domain\n"
        "   -v, --velocity_left   Set velocity on left side of domain\n"
        "   -V, --velocity_right  Set velocity on right side of domain\n"
        "   -p, --pressure_left   Set pressure on left side of domain\n"
        "   -P, --pressure_right  Set pressure on right side of domain\n"
        "   -o, --outfile         File to write solution\n"
        "   -h, --help            Show help\n";
    exit(1);
}

// https://gist.github.com/ashwin/d88184923c7161d368a9
void ProcessArgs(int argc, char** argv)
{
    const char* const short_opts = "n:t:d:D:v:V:p:P:o:h";
    const option long_opts[] = {
        {"steps", required_argument, nullptr, 'n'},
        {"time", required_argument, nullptr, 't'},
        {"density_left", required_argument, nullptr, 'd'},
        {"density_right", required_argument, nullptr, 'D'},
        {"velocity_left", required_argument, nullptr, 'v'},
        {"velocity_right", required_argument, nullptr, 'V'},
        {"pressure_left", required_argument, nullptr, 'p'},
        {"pressure_right", required_argument, nullptr, 'P'},
        {"outfile", required_argument, nullptr, 'o'},
        {"help", no_argument, nullptr, 'h'},
        {nullptr, no_argument, nullptr, 0}
    };

    while (true)
    {
        const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

        if (-1 == opt)
            break;

        switch (opt)
        {
        case 'n':
            ncells = std::stoi(optarg);
            break;
        case 't':
           maxtime = std::stof(optarg);
           break;
        case 'd':
            densityL = std::stof(optarg);
            break;
        case 'D':
            densityR = std::stof(optarg);
            break;
        case 'v':
            velocityL = std::stof(optarg);
            break;
        case 'V':
            velocityR = std::stof(optarg);
            break;
        case 'p':
            pressureL = std::stof(optarg);
            break;
        case 'P':
            pressureR = std::stof(optarg);
            break;
        case 'o':
            outPath = std::string(optarg);
            break;
        case 'h': // -h or --help
        case '?': // Unrecognized option
        default:
            PrintHelp();
            break;
        }
    }
}

double get_energy(double pressure, double density, 
                  double velocity, double gamma = 1.4) {
    double internal_energy = pressure/((gamma - 1)*density);
    return (density * internal_energy) + (0.5*density*velocity*velocity);
}

double get_velocity(std::map<std::string, double>& q) {
    return q["momentum"] / q["density"];
}

double get_internal_energy(std::map<std::string, double>& q) {
    double velocity = get_velocity(q);
    return (q["energy"] / q["density"]) - (0.5 * velocity * velocity);
}

double get_pressure(std::map<std::string, double>& q, double gamma = 1.4) {;
    double internal_energy = get_internal_energy(q);
    return (gamma - 1) * (q["density"] * internal_energy);
}

double delta_x(double ncells) {
  return 1/ncells;
}

double speed_of_sound(std::map<std::string, double>& q, double gamma = 1.4) {
  double pressure = get_pressure(q);
  double c_squared = (gamma*pressure)/q["density"];            
  return std::abs(sqrt(c_squared));
}

double get_amax(std::vector<std::map<std::string, double>>& q) {
  int size = q.size();
  double c;
  std::vector<double> all_a;
  for (int i = 0; i < size; i++) {
    double velocity = get_velocity(q[i]);
    c = speed_of_sound(q[i]);
    all_a.push_back(std::abs(velocity) + c);
  } 
  return *std::max_element(all_a.begin(), all_a.end());
}

double get_dt(double amax, double dx, double CFL = 0.9) {
  return CFL*(dx/amax);
}

double f_density(std::map<std::string, double>& q) {
  double velocity = get_velocity(q);
  return q["density"] * velocity;
}

double f_momentum(std::map<std::string, double>& q) {
    double velocity = get_velocity(q);
    double pressure = get_pressure(q);
    return (q["density"]*velocity*velocity) + pressure;
}

double f_energy(std::map<std::string, double>& q) {
  double velocity = get_velocity(q);
  double pressure = get_pressure(q);
  return (q["energy"] + pressure)*velocity;
}

std::map<std::string, double> get_q_ihalf(std::map<std::string, double>& qi , 
                                          std::map<std::string, double>& q_i1, 
                                          double dx, double dt) {
                          
    // Declare output vector
    std::map<std::string, double> q_out;
    
    q_out["density"] = (0.5 * (qi["density"] + q_i1["density"])) 
               + (0.5 * (dt/dx) * (f_density(qi) - f_density(q_i1)));
    q_out["momentum"] = (0.5 * (qi["momentum"] + q_i1["momentum"])) 
               + (0.5 * (dt/dx) * (f_momentum(qi) - f_momentum(q_i1)));
    q_out["energy"] = (0.5 * (qi["energy"] + q_i1["energy"])) 
               + (0.5 * (dt/dx) * (f_energy(qi) - f_energy(q_i1)));
               
    return q_out;
}

std::map<std::string, double>  get_flux_RI(std::map<std::string, double>& q) {
    // Declare output vector
    std::map<std::string, double> flux_RI;
  
    flux_RI["density"] = f_density(q);
    flux_RI["momentum"] = f_momentum(q);
    flux_RI["energy"] = f_energy(q);
    
    return flux_RI;
}

std::map<std::string, double> get_flux_LF(std::map<std::string, double>& qi , 
                                          std::map<std::string, double>& q_i1, 
                                          double dx, double dt) {
                          
    // Declare output vector
    std::map<std::string, double> flux_LF;
    
    flux_LF["density"] = (0.5 * (f_density(qi) + f_density(q_i1))) 
                 + (0.5 * (dx/dt) * (qi["density"] - q_i1["density"]));
    flux_LF["momentum"] = (0.5 * (f_momentum(qi) + f_momentum(q_i1))) 
                 + (0.5 * (dx/dt) * (qi["momentum"] - q_i1["momentum"]));
    flux_LF["energy"] = (0.5 * (f_energy(qi) + f_energy(q_i1))) 
                 + (0.5 * (dx/dt) * (qi["energy"] - q_i1["energy"]));
    
    return flux_LF;
}

std::map<std::string, double> get_force(std::map<std::string, double>& flux_LF, 
                                        std::map<std::string, double>& flux_RI) {
                          
    // Declare output vector
    std::map<std::string, double> force_flux;
    
    force_flux["density"] = 0.5 * (flux_LF["density"] + flux_RI["density"]);
    force_flux["momentum"] = 0.5 * (flux_LF["momentum"] + flux_RI["momentum"]);
    force_flux["energy"] = 0.5 * (flux_LF["energy"] + flux_RI["energy"]);
    
    return force_flux;
}

void initiaseData(std::vector<std::map<std::string, double>>& q,
                  double densityL, double velocityL, double pressureL,
                  double densityR, double velocityR, double pressureR) {
  
  int size = q.size();
  double energy, density, velocity, pressure, momentum;
  for (int i = 0; i < size; i++) {
    if (i < (size/2)) {
      density = densityL;
      velocity = velocityL;
      pressure = pressureL;
    }
    else {
      density = densityR;
      velocity = velocityR;
      pressure = pressureR;
    }
    energy = get_energy(pressure, density, velocity);
    momentum = density * velocity;
    q[i]["density"] = density;
    q[i]["momentum"] = momentum;
    q[i]["energy"] = energy;    
  }
}
      

int main(int argc, char* argv[])
{
  ProcessArgs(argc, argv);
  std::ofstream outFile;
  if (outPath != "-") {
    outFile.open(outPath, std::ios::out);
  } 
  std::ostream& out = (outPath != "-" ? outFile : std::cout);
    
  double dx = delta_x(ncells);
  std::vector<std::map<std::string, double>> q(ncells);
  std::vector<std::map<std::string, double>> q_next(ncells);
  initiaseData(q, densityL, velocityL, pressureL,
               densityR, velocityR, pressureR);
  double t = 0;
  double amax, dt;
  while (t < maxtime) {
    
    amax = get_amax(q);
    dt = get_dt(amax, dx);
    int size = q.size();
    std::map<std::string, double> q_ihalf, q_n1, flux_RI, flux_LF, force_flux, 
                               q_ihalf_less, flux_RI_less, flux_LF_less, 
                               force_flux_less;
    for (int i = 0; i < size; i++) {

        std::map<std::string, double> q_1less;
        std::map<std::string, double> q_1plus;
        if (i == 0) {
            q_1less = q[i];
        } else {
            q_1less = q[i-1];   
        }
        
        if (i == size - 1) {
            q_1plus = q[size - 1];
        } else {
            q_1plus = q[i+1];   
        }
            
        q_ihalf = get_q_ihalf(q[i], q_1plus, dx, dt);
        flux_RI = get_flux_RI(q_ihalf);
        flux_LF = get_flux_LF(q[i], q_1plus, dx, dt);
        force_flux = get_force(flux_LF, flux_RI);
        
        q_ihalf_less = get_q_ihalf(q_1less, q[i], dx, dt);
        flux_RI_less = get_flux_RI(q_ihalf_less);
        flux_LF_less = get_flux_LF(q_1less, q[i], dx, dt);
        force_flux_less = get_force(flux_LF_less, flux_RI_less);
        
        q_n1["density"] = q[i]["density"] + (dt/dx)*(force_flux_less["density"] - force_flux["density"]);
        q_n1["momentum"] = q[i]["momentum"] + (dt/dx)*(force_flux_less["momentum"] - force_flux["momentum"]);
        q_n1["energy"] = q[i]["energy"] + (dt/dx)*(force_flux_less["energy"] - force_flux["energy"]);

        q_next[i] = q_n1;
    }
    t += dt;
    q = q_next;
    }
    
    int size = q.size();
    for (int i = 0; i < size; i++) {
      out << q[i]["density"] << "\t" 
          << get_velocity(q[i]) << "\t" 
          << get_pressure(q[i]) << "\t" 
          << get_internal_energy(q[i]) << std::endl;
    }
    return 0;
}
