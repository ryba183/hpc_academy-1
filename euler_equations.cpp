#include <iostream>
#include <cmath>
#include <complex>
#include <array>
#include <vector>
#include <functional>
#include <algorithm>

double get_energy(double pressure, double density, 
                  double velocity, double gamma = 1.4) {
    double internal_energy = pressure/((gamma - 1)*density);
    return (density * internal_energy) + (0.5*density*velocity*velocity);
}

double get_velocity(std::array<double, 3>& q) {
    double density = q[0];
    double momentum = q[1];
    return momentum / density;
}

double get_internal_energy(std::array<double, 3>& q) {
    double density = q[0];
    double energy = q[2]; 
    double velocity = get_velocity(q);
    return (energy / density) - (0.5 * velocity * velocity);
}

double get_pressure(std::array<double, 3>& q, double gamma = 1.4) {
    double density = q[0];
    double internal_energy = get_internal_energy(q);
    return (gamma - 1) * (density * internal_energy);
}

double delta_x(double ncells) {
  return 1/ncells;
}

double speed_of_sound(std::array<double, 3>& q, double gamma = 1.4) {
  double density = q[0];
  double pressure = get_pressure(q);
  double c_squared = (gamma*pressure)/density;            
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

double f_density(std::array<double, 3> q) {
  double density = q[0];
  double velocity = get_velocity(q);
  return density * velocity;
}

double f_momentum(std::array<double, 3> q) {
    double density = q[0];
    double velocity = get_velocity(q);
    double pressure = get_pressure(q);
    return (density*velocity*velocity) + pressure;
}

double f_energy(std::array<double, 3> q) {
  double energy = q[2];
  double velocity = get_velocity(q);
  double pressure = get_pressure(q);
  return (energy + pressure)*velocity;
}

std::array<double, 3> get_q_ihalf(std::array<double, 3> qi , 
                                  std::array<double, 3> q_i1, 
                                  double dx, double dt) {
                          
    // Declare output vector
    std::array<double, 3> q_out;
    
    // Get variable values for i and i+1
    double density = qi[0];
    double momentum = qi[1];
    double energy = qi[2];
    double density_i1 = q_i1[0];
    double momentum_i1 = q_i1[1];
    double energy_i1 = q_i1[2]; 

    q_out[0] = (0.5 * (density + density_i1)) 
               + (0.5 * (dt/dx) * (f_density(qi) - f_density(q_i1)));
    q_out[1] = (0.5 * (momentum + momentum_i1)) 
               + (0.5 * (dt/dx) * (f_momentum(qi) - f_momentum(q_i1)));
    q_out[2] = (0.5 * (energy + energy_i1)) 
               + (0.5 * (dt/dx) * (f_energy(qi) - f_energy(q_i1)));
               
    return q_out;
}

std::array<double, 3>  get_flux_RI(std::array<double, 3> q) {
    // Declare output vector
    std::array<double, 3> flux_RI;
  
    flux_RI[0] = f_density(q);
    flux_RI[1] = f_momentum(q);
    flux_RI[2] = f_energy(q);
    
    return flux_RI;
}

std::array<double, 3> get_flux_LF(std::array<double, 3> qi , 
                                  std::array<double, 3> q_i1, 
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

std::array<double, 3> get_force(std::array<double, 3> flux_LF, 
                                std::array<double, 3> flux_RI) {
                          
    // Declare output vector
    std::array<double, 3> force_flux;
    
    force_flux[0] = 0.5 * (flux_LF[0] + flux_RI[0]);
    force_flux[1] = 0.5 * (flux_LF[1] + flux_RI[1]);
    force_flux[2] = 0.5 * (flux_LF[2] + flux_RI[2]);
    
    return force_flux;
}

void initiaseData(std::vector<std::array<double, 3>>& q,
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
    q[i][0] = density;
    q[i][1] = momentum;
    q[i][2] = energy;    
  }
}
      

int main(int argc, char* argv[])
{
  double densityL = 1;
  double velocityL = -2;
  double pressureL = 0.4;
  double densityR = 1;
  double velocityR = 2;
  double pressureR = 0.4;
  double time = 0.15;
  double ncells = 1000;
  double dx = delta_x(ncells);
  std::vector<std::array<double, 3>> q(ncells);
  std::vector<std::array<double, 3>> q_next(ncells);
  initiaseData(q, densityL, velocityL, pressureL,
               densityR, velocityR, pressureR);
  double t = 0;
  double amax, dt;
  while (t < time) {
    
    amax = get_amax(q);
    dt = get_dt(amax, dx);
    int size = q.size();
    std::array<double, 3> q_ihalf, q_n1, flux_RI, flux_LF, force_flux, 
                          q_ihalf_less, flux_RI_less, flux_LF_less, 
                          force_flux_less;
    for (int i = 0; i < size; i++) {

        std::array<double, 3> q_1less;
        std::array<double, 3> q_1plus;
        if (i == 0) {
            q_1less = {q[i][0],q[i][1],q[i][2]};
        } else {
            q_1less = {q[i-1][0],q[i-1][1],q[i-1][2]};   
        }
        
        if (i == size - 1) {
            q_1plus = {q[size - 1][0],q[size - 1][1],q[size - 1][2]};
        } else {
            q_1plus = {q[i+1][0],q[i+1][1],q[i+1][2]};   
        }
            
        q_ihalf = get_q_ihalf(q[i], q_1plus, dx, dt);
        flux_RI = get_flux_RI(q_ihalf);
        flux_LF = get_flux_LF(q[i], q_1plus, dx, dt);
        force_flux = get_force(flux_LF, flux_RI);
        
        q_ihalf_less = get_q_ihalf(q_1less, q[i], dx, dt);
        flux_RI_less = get_flux_RI(q_ihalf_less);
        flux_LF_less = get_flux_LF(q_1less, q[i], dx, dt);
        force_flux_less = get_force(flux_LF_less, flux_RI_less);
        
        q_n1[0] = q[i][0] + (dt/dx)*(force_flux_less[0] - force_flux[0]);
        q_n1[1] = q[i][1] + (dt/dx)*(force_flux_less[1] - force_flux[1]);
        q_n1[2] = q[i][2] + (dt/dx)*(force_flux_less[2] - force_flux[2]);

        q_next[i] = q_n1;
    }
    t += dt;
    q = q_next;
    }
    
    int size = q.size();
    for (int i = 0; i < size; i++) {
      std::cout << q[i][0] << "\t" 
                << get_velocity(q[i]) << "\t" 
                << get_pressure(q[i]) << "\t" 
                << get_internal_energy(q[i]) << std::endl;
    }
    return 0;
}
