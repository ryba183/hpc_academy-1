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

unsigned int ncells = 1000;
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

double get_velocity(std::array<double, 3>& q) {
    return q[1] / q[0];
}

double get_internal_energy(std::array<double, 3>& q) {
    double velocity = get_velocity(q);
    return (q[2] / q[0]) - (0.5 * velocity * velocity);
}

double get_pressure(std::array<double, 3>& q, double gamma = 1.4) {;
    double internal_energy = get_internal_energy(q);
    return (gamma - 1) * (q[0] * internal_energy);
}

double get_dx(unsigned int ncells, double domain_start = 0, double domain_end = 1) {
  return (domain_end - domain_start)/ncells;
}

double speed_of_sound(std::array<double, 3>& q, double gamma = 1.4) {
  double pressure = get_pressure(q);
  double c_squared = (gamma*pressure)/q[0];            
  return std::abs(sqrt(c_squared));
}

double get_amax(std::vector<std::array<double, 3>>& q) {
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

double f_density(std::array<double, 3>& q) {
  double velocity = get_velocity(q);
  return q[0] * velocity;
}

double f_momentum(std::array<double, 3>& q) {
    double velocity = get_velocity(q);
    double pressure = get_pressure(q);
    return (q[0]*velocity*velocity) + pressure;
}

double f_energy(std::array<double, 3>& q) {
  double velocity = get_velocity(q);
  double pressure = get_pressure(q);
  return (q[2] + pressure)*velocity;
}

std::array<double, 3> get_q_ihalf(std::array<double, 3>& qi , 
                                  std::array<double, 3>& q_i1, 
                                  double dx, double dt) {
                          
    // Declare output vector
    std::array<double, 3> q_out;
    
    q_out[0] = (0.5 * (qi[0] + q_i1[0])) 
               + (0.5 * (dt/dx) * (f_density(qi) - f_density(q_i1)));
    q_out[1] = (0.5 * (qi[1] + q_i1[1])) 
               + (0.5 * (dt/dx) * (f_momentum(qi) - f_momentum(q_i1)));
    q_out[2] = (0.5 * (qi[2] + q_i1[2])) 
               + (0.5 * (dt/dx) * (f_energy(qi) - f_energy(q_i1)));
               
    return q_out;
}

std::array<double, 3>  get_flux_RI(std::array<double, 3>& q) {
    // Declare output vector
    std::array<double, 3> flux_RI;
  
    flux_RI[0] = f_density(q);
    flux_RI[1] = f_momentum(q);
    flux_RI[2] = f_energy(q);
    
    return flux_RI;
}

std::array<double, 3> get_flux_LF(std::array<double, 3>& qi , 
                                          std::array<double, 3>& q_i1, 
                                          double dx, double dt) {
                          
    // Declare output vector
    std::array<double, 3> flux_LF;
    
    flux_LF[0] = (0.5 * (f_density(qi) + f_density(q_i1))) 
                 + (0.5 * (dx/dt) * (qi[0] - q_i1[0]));
    flux_LF[1] = (0.5 * (f_momentum(qi) + f_momentum(q_i1))) 
                 + (0.5 * (dx/dt) * (qi[1] - q_i1[1]));
    flux_LF[2] = (0.5 * (f_energy(qi) + f_energy(q_i1))) 
                 + (0.5 * (dx/dt) * (qi[2] - q_i1[2]));
    
    return flux_LF;
}

std::array<double, 3> get_force(std::array<double, 3>& qi, 
                                        std::array<double, 3>& q_i1,
                                        double dx, double dt) {
                          
    std::array<double, 3>  q_ihalf = get_q_ihalf(qi, q_i1, dx, dt);
    std::array<double, 3>  flux_RI = get_flux_RI(q_ihalf);
    std::array<double, 3>  flux_LF = get_flux_LF(qi, q_i1, dx, dt);
    
    // Declare output vector
    std::array<double, 3> force_flux;
    
    force_flux[0] = 0.5 * (flux_LF[0] + flux_RI[0]);
    force_flux[1] = 0.5 * (flux_LF[1] + flux_RI[1]);
    force_flux[2] = 0.5 * (flux_LF[2] + flux_RI[2]);
    
    return force_flux;
}

std::array<double, 3> get_q_n1(std::array<double, 3>& qi, 
                                       std::array<double, 3>& q_iplus1,
                                       std::array<double, 3>& q_iless1,
                                       double dx, double dt) {
    
    std::array<double, 3> force_plus_half = get_force(qi, q_iplus1, dx, dt);
    std::array<double, 3> force_less_half = get_force(q_iless1, qi, dx, dt);

    // Declare output vector
    std::array<double, 3> q_nplus1;
    
    q_nplus1[0] = qi[0] + (dt/dx)*(force_less_half[0] - force_plus_half[0]);
    q_nplus1[1] = qi[1] + (dt/dx)*(force_less_half[1] - force_plus_half[1]);
    q_nplus1[2] = qi[2] + (dt/dx)*(force_less_half[2] - force_plus_half[2]);
    
    return q_nplus1;
}

std::array<std::array<double, 3>, 2> get_q_isurround(std::vector<std::array<double, 3>>& q, int i) {
    // Returns q[i-1] and q[i+1] accounting for boundary conditions
    
    // Declare output map
    std::array<std::array<double, 3>, 2> q_isurround;

    int size = q.size();
    if (i == 0) {
        q_isurround[0] = q[i];
    } else {
        q_isurround[0] = q[i-1];   
    }
    if (i == size - 1) {
        q_isurround[1] = q[size - 1];
    } else {
        q_isurround[1] = q[i+1];   
    }
    return q_isurround;
}


std::vector<std::array<double, 3>> initialiseData(unsigned int ncells,
    double densityL, double velocityL, double pressureL,
    double densityR, double velocityR, double pressureR) {
    
    // Declare initial data
    std::vector<std::array<double, 3>> q(ncells);
    // std::vector<std::array<double, 3>>
    double energy, density, velocity, pressure, momentum;
    for (unsigned int i = 0; i < ncells; i++) {
        if (i < (ncells/2)) {
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
        q[i][0] = density;
        q[i][1] = momentum;
        q[i][2] = energy;    
    }
    return q;
}
      

int main(int argc, char* argv[]) {
    
    // Extract command line argumenets
    ProcessArgs(argc, argv);
    std::ofstream outFile;
    
    // Set output to file path or stdout (default)
    if (outPath != "-") {
        outFile.open(outPath, std::ios::out);
    } 
    std::ostream& out = (outPath != "-" ? outFile : std::cout);
    
    // Initialise data
    std::vector<std::array<double, 3>> q = initialiseData(ncells, densityL, velocityL, pressureL, densityR, velocityR, pressureR);
    
    // Compute delta x across domain space
    double dx = get_dx(ncells);
    
    double t = 0;
    while (t < maxtime) {
        
        // Declare vector for storing q(n+1) for all i
        std::vector<std::array<double, 3>> q_next(ncells);
        
        // Compute dt for q[i]
        double amax = get_amax(q);
        double dt = get_dt(amax, dx);
        
        for (unsigned int i = 0; i < q.size(); i++) {
            // Return q(n, i+1) and q(n, i-1)
            // May use std::array<std::array<double, 3>, 2> instead of map to store q_surrond for increased speed!
            std::array<std::array<double, 3>, 2> q_surround = get_q_isurround(q, i);
            // Compute q(i, n+1)
            q_next[i] = get_q_n1(q[i], q_surround[1], q_surround[0], dx, dt);
        }
        t += dt;
        q = q_next;
    }
    
    int size = q.size();
    for (int i = 0; i < size; i++) {
        out << q[i][0] << "\t" 
            << get_velocity(q[i]) << "\t" 
            << get_pressure(q[i]) << "\t" 
            << get_internal_energy(q[i]) << std::endl;
    }
    return 0;
}
