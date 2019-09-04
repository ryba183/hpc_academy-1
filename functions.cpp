#include <iostream>
#include <cmath>
#include <complex>

int sum(int a, int b) {
  return a + b;
}

// Take a by reference to modify a non-locally
void decrease(int& a, int b) {
  a = a - b;
}

int q(double a, double b, double c, double& x1, double& x2) {
  double disc = (b*b) - (4*c*a);
  if(disc < 0) {
    return 0;
  }
  else { 
    double disc_sqrt = sqrt(disc);
    if(b > 0) {
      x1 = (-b - disc_sqrt)/(2*a);
      x2 = (2*c) / (-b - disc_sqrt);
    }
    else {
      x1 = (-b + disc_sqrt)/(2*a);
      x2 = (2*c) / (-b + disc_sqrt);
    }
  }
  return 2;
}

int main(void) {
  int a = 1;
  int a0 = a;
  int b = 2;
  
  int sum_ab = sum(a, b);
  std::cout << a << " + " << b << " = " << sum_ab << std::endl;
  
  decrease(a, b);
  std::cout << a0 << " - " << b << " = " << a << std::endl;
  
  double A = 1;
  double B = 2;
  double C = 1;
  double x1, x2;
  
  int k = q(A, B, C, x1, x2);
  std::cout << k  << " solution(s) for a = " << A << ", b = " << B << ", c = " << C << std::endl;
  std::cout << "Solutions are : " << x1 << " and " << x2 << std::endl;
  
}
