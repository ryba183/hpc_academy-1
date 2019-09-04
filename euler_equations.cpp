#include <iostream>
#include <cmath>
#include <complex>
#include <array>
#include <vector>
#include <functional>

double get_epsilon(double density, double pressure, double gamma = 1.4) {
  return pressure/((gamma - 1)*density);
}

double get_energy(double epsilon, double density, double velocity) {
  return (density * epsilon) + (0.5*density*velocity*velocity);
}

double f_mass (double ro, double v, double p) {
  return (ro*v*v) + p;
}

double f_energy (double energy, double p, double v) {
  return (energy + p)*v;
}

double delta_x (double ncells) {
  return 1/ncells;
}

double speed_of_sound(double density, double pressure, double gamma = 1.4) {
  double c_squared = (gamma*pressure)/density;
  return sqrt(c_squared);
}

double get_amax (std::vector<std::array<double, 4>>& q) {
  int size = q.size();
  double a, c, density, velocity, pressure, amax;
  for (int i = 0; i < size; i++) {
    density = q[i][0];
    velocity = q[i][1];
    pressure = q[i][4];
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

double delta_t (double amax, double dx, double CFL = 0.9) {
  return CFL*(dx/amax);
}

double f_density(double density, double velocity) {
  return density * velocity;
}

double q_ihalf_density(double density_i, double density_i1, 
                       double velocity_i, double velocity_i1, 
                       double dx, double dt) {
  return 0.5 * (density_i + density_i1) 
         + (0.5 * (dt/dx) * (f_density(density_i, velocity_i) 
                             - f_density(density_i1, velocity_i1)));
}


// Compute qhalf for pressure and velocity

// Compute all fluxes with qi, qhalf and qi1



double force_flux(std::array<double, 4> q_i , 
                  std::array<double, 4> q_i1, 
                  double dx, double dt) {
  double density_i = q_i[0];
  double velocity_i = q_i[1];
  double energy_i = q_i[2];
  double pressure_i = q_i[3];
  double density_i1 = q_i1[0];
  double velocity_i1 = q_i1[1];
  double energy_i1 = q_i[2];
  double pressure_i1 = q_i[3];

  double density_ihalf = q_ihalf_density(density_i, density_i1,
                                         velocity_i, velocity_i1,
                                         dx, dt);
}

void initiaseData(std::vector<std::array<double, 4>>& q,
                  double densityL, double velocityL, double pressureL,
                  double densityR, double velocityR, double pressureR) {
  
  int size = q.size();
  double energy, epsilon, density, velocity, pressure;
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
    epsilon = get_epsilon(density, pressure);
    energy = get_energy(epsilon, density, velocity);
    q[i][0] = density;
    q[i][1] = velocity;
    q[i][2] = energy;
    q[i][3] = pressure;
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
  std::vector<std::array<double, 4>> q(ncells);
  
  initiaseData(q, densityL, velocityL, pressureL,
               densityR, velocityR, pressureR);
  double t = 0;
  
  double amax, dt;
  while (t < time) {
    amax = get_amax(q);
    dt = delta_t(amax, dx);
    
  }

  std::cout << dx << std::endl;
  std::cout << dt << std::endl;

  return 0;
}
