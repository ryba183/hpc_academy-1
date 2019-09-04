#include <iostream>
#include <cmath>
#include <complex>
#include <array>
#include <vector>
#include <functional>


double get_energy(double pressure, double density, double velocity, double gamma = 1.4) {
    double internal_energy = pressure/((gamma - 1)*density);
    return (density * internal_energy) + (0.5*density*velocity*velocity);
}

double get_pressure(double energy, double density, double velocity, double gamma = 1.4) {
    double internal_energy = (energy - (0.5 * density * velocity * velocity))/density;
    return (gamma - 1) * (density * internal_energy);
}

double delta_x(double ncells) {
  return 1/ncells;
}

double speed_of_sound(double density, double pressure, double gamma = 1.4) {
  double c_squared = (gamma*pressure)/density;
  return sqrt(c_squared);
}

double get_amax(std::vector<std::array<double, 5>>& q) {
  int size = q.size();
  double a, c, density, velocity, pressure, amax;
  for (int i = 0; i < size; i++) {
    density = q[i][0];
    velocity = q[i][1];
    pressure = q[i][3];
    c = speed_of_sound(density, pressure);
    a = std::abs(velocity) + c;
    if (i == 0) {
      amax = a;
    }
    else { 
      double amax = std::max(amax, a);
    }
  }
  return amax;
}

double delta_t(double amax, double dx, double CFL = 0.9) {
  return CFL*(dx/amax);
}

double f_momentum(double density, double velocity, double pressure) {
  return (density*velocity*velocity) + pressure;
}

double get_momentum_ihalf(std::array<double, 5> q_i , 
                        std::array<double, 5> q_i1, 
                        double dx, double dt) {
    double density_i = q_i[0];
    double velocity_i = q_i[1];
    double pressure_i = q_i[3];
    double momentum_i = q_i[4];
    double density_i1 = q_i1[0];
    double velocity_i1 = q_i1[1];
    double pressure_i1 = q_i1[3];
    double momentum_i1 = q_i1[4];
    return 0.5 * (momentum_i + momentum_i1) 
            + (0.5 * (dt/dx) * (f_momentum(density_i, velocity_i, pressure_i) 
                             - f_momentum(density_i1, velocity_i1, pressure_i1)));
}

double get_flux_momentum(std::array<double, 5> q_i , 
                         std::array<double, 5> q_i1, 
                         double dx, double dt,
                         double q_ihalf_momentum,
                         double q_ihalf_density,
                         double q_ihalf_velocity,
                         double q_ihalf_pressure) {
    double density_i = q_i[0];
    double velocity_i = q_i[1];
    double pressure_i = q_i[3];
    double momentum_i = q_i[4];
    double density_i1 = q_i1[0];
    double velocity_i1 = q_i1[1];
    double pressure_i1 = q_i1[3];
    double momentum_i1 = q_i1[4];
    
    double f_RI = f_momentum(q_ihalf_density, q_ihalf_velocity, q_ihalf_pressure);
    double f_LF = 0.5 * (f_momentum(density_i, velocity_i, pressure_i) 
                  + f_momentum(density_i1, velocity_i1, pressure_i1))
                  + (0.5 * (dx/dt) * (momentum_i - momentum_i1));
    double f_FORCE = 0.5 * (f_LF + f_RI);
    
    return f_FORCE;
}

double f_density(double density, double velocity) {
  return density * velocity;
}

double get_density_ihalf(std::array<double, 5> q_i , 
                        std::array<double, 5> q_i1, 
                        double dx, double dt) {
    double density_i = q_i[0];
    double velocity_i = q_i[1];
    double density_i1 = q_i1[0];
    double velocity_i1 = q_i1[1];
    return 0.5 * (density_i + density_i1) 
            + (0.5 * (dt/dx) * (f_density(density_i, velocity_i) 
                             - f_density(density_i1, velocity_i1)));
}

double get_flux_density(std::array<double, 5> q_i , 
                         std::array<double, 5> q_i1, 
                         double dx, double dt,
                         double q_ihalf_density,
                         double q_ihalf_velocity) {
    double density_i = q_i[0];
    double velocity_i = q_i[1];
    double energy_i = q_i[2];
    double pressure_i = q_i[3];
    double density_i1 = q_i1[0];
    double velocity_i1 = q_i1[1];
    double energy_i1 = q_i1[2];
    double pressure_i1 = q_i1[3];

    double f_RI = f_density(q_ihalf_density, q_ihalf_velocity);
    double f_LF = 0.5 * (f_density(density_i,  velocity_i) 
                  + f_density(density_i1,  velocity_i1) )
                  + (0.5 * (dx/dt) * (density_i - density_i1));
    double f_FORCE = 0.5 * (f_LF + f_RI);
    
    return f_FORCE;
}

double f_energy(double energy, double pressure, double velocity) {
  return (energy + pressure)*velocity;
}

double get_energy_ihalf(std::array<double, 5> q_i , 
                        std::array<double, 5> q_i1, 
                        double dx, double dt) {
    double density_i = q_i[0];
    double velocity_i = q_i[1];
    double energy_i = q_i[2];
    double pressure_i = q_i[3];
    double density_i1 = q_i1[0];
    double velocity_i1 = q_i1[1];
    double energy_i1 = q_i1[2];
    double pressure_i1 = q_i1[3];
    return 0.5 * (energy_i + energy_i1) 
            + (0.5 * (dt/dx) * (f_energy(energy_i, pressure_i, velocity_i) 
                             - f_energy(energy_i1, pressure_i1, velocity_i1)));
}

double get_flux_energy(std::array<double, 5> q_i , 
                         std::array<double, 5> q_i1, 
                         double dx, double dt,
                         double q_ihalf_momentum,
                         double q_ihalf_density,
                         double q_ihalf_energy,
                         double q_ihalf_velocity,
                         double q_ihalf_pressure) {
    double density_i = q_i[0];
    double velocity_i = q_i[1];
    double energy_i = q_i[2];
    double pressure_i = q_i[3];
    double density_i1 = q_i1[0];
    double velocity_i1 = q_i1[1];
    double energy_i1 = q_i1[2];
    double pressure_i1 = q_i1[3];
    
    

    double f_RI = f_energy(q_ihalf_energy, q_ihalf_pressure, q_ihalf_velocity);
    double f_LF = 0.5 * (f_energy(energy_i, pressure_i, velocity_i) 
                  + f_energy(energy_i1, pressure_i1, velocity_i1))
                  + (0.5 * (dx/dt) * (energy_i - energy_i1));
    double f_FORCE = 0.5 * (f_LF + f_RI);
    
    return f_FORCE;
}

void initiaseData(std::vector<std::array<double, 5>>& q,
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
    q[i][1] = velocity;
    q[i][2] = energy;
    q[i][3] = pressure;
    q[i][4] = momentum;
  }
}
      

int main(int argc, char* argv[])
{
  double densityL = 1;
  double velocityL = 0;
  double pressureL = 1;
  double densityR = 0.125;
  double velocityR = 0;
  double pressureR = 0.1;
  double time = 0.25;
  double ncells = 100;
  double dx = delta_x(ncells);
  std::vector<std::array<double, 5>> q(ncells);
  std::vector<std::array<double, 5>> q_next(ncells);
  initiaseData(q, densityL, velocityL, pressureL,
               densityR, velocityR, pressureR);
  double t = 0;
  
  double amax, dt;
  while (t < time) {
    amax = get_amax(q);
    dt = delta_t(amax, dx);
    int size = q.size();
    double momentum_ihalf, density_ihalf, energy_ihalf, velocity_ihalf, pressure_ihalf;
    std::array<double, 5> q_i, q_i1;
    for (int i = 0; i < size; i++) {
        
        std::array<double, 5> q_1less;
        if (i == 0) {
            q_1less = {q[i][0],q[i][1],q[i][2],q[i][3],q[i][4]};
        } else {
            q_1less = {q[i-1][0],q[i-1][1],q[i-1][2],q[i-1][3],q[i-1][4]};
        }
            
        
        q_i = q[i];
        q_i1 = q[i+1];
        double momentum_ihalf = get_momentum_ihalf(q_i, q_i1, dx, dt);
        double density_ihalf = get_density_ihalf(q_i, q_i1, dx, dt);
        double energy_ihalf = get_energy_ihalf(q_i, q_i1, dx, dt);
        double velocity_ihalf = momentum_ihalf / density_ihalf;
        double pressure_ihalf = get_pressure(energy_ihalf, density_ihalf, velocity_ihalf);
        
        double flux_energy = get_flux_energy(q_i, q_i1, dx, dt,
            momentum_ihalf, density_ihalf,
            energy_ihalf, velocity_ihalf,
            pressure_ihalf);
                         
        double flux_density = get_flux_density(q_i, q_i1, dx, dt,
            density_ihalf,
            velocity_ihalf);
                         
        double flux_momentum = get_flux_momentum(q_i, q_i1, dx, dt,
            momentum_ihalf, density_ihalf,
            velocity_ihalf, pressure_ihalf);
            
        // Less 1
        double momentum_ihalf_1less = get_momentum_ihalf(q_1less, q_i, dx, dt);
        double density_ihalf_1less = get_density_ihalf(q_1less, q_i, dx, dt);
        double energy_ihalf_1less = get_energy_ihalf(q_1less, q_i, dx, dt);
        double velocity_ihalf_1less = momentum_ihalf / density_ihalf;
        double pressure_ihalf_1less = get_pressure(energy_ihalf, density_ihalf, velocity_ihalf);
        
        double flux_energy_1less = get_flux_energy(q_1less, q_i, dx, dt,
            momentum_ihalf, density_ihalf,
            energy_ihalf, velocity_ihalf,
            pressure_ihalf);
                         
        double flux_density_1less = get_flux_density(q_1less, q_i, dx, dt,
            density_ihalf,
            velocity_ihalf);
                         
        double flux_momentum_1less = get_flux_momentum(q_1less, q_i, dx, dt,
            momentum_ihalf, density_ihalf,
            velocity_ihalf, pressure_ihalf);


        double q_1_energy = (dt/dx) * (flux_energy_1less - flux_energy);
        double q_1_density = (dt/dx) * (flux_density_1less - flux_density);
        double q_1_momentum = (dt/dx) * (flux_momentum_1less - flux_momentum);
        double q_1_velocity = q_1_momentum / q_1_density;
        double q_1_pressure = get_pressure(q_1_energy, q_1_density, q_1_velocity);
        q_next[i][0] = q_1_density;
        q_next[i][1] = q_1_velocity;
        q_next[i][2] = q_1_energy;
        q_next[i][3] = q_1_pressure;
        q_next[i][4] = q_1_momentum;
        
        std::cout << q_1_density << std::endl;
        std::cout << q_1_velocity << std::endl; // returns nan
        std::cout << q_1_energy << std::endl;
        std::cout << q_1_pressure << std::endl; //returns nan
        std::cout << q_1_momentum << std::endl;
        t += dt;
        q = q_next;
    }
    }
    return 0;
}
