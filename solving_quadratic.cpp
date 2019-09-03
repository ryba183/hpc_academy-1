#include <iostream>
#include <cmath>
#include <complex>
int main(void)
{
  int a, b, c;
  std::cout << "Enter a:" << std::flush;
  std::cin >> a;
  std::cout << "Enter b:" << std::flush;
  std::cin >> b;
  std::cout << "Enter c:" << std::flush;
  std::cin >> c;
  
  double s1, s2;
  double j = (b*b) - (4*c*a);
  if(j < 0) {
    std::cout << "There are no real solutions" << std::endl;
  }
  else { 
    j = sqrt(j);
    if(b > 0) {
      s1 = (-b - j)/(2*a);
      s2 = (2*c) / (-b -j);
    }
    else {
      s1 = (-b + j)/(2*a);
      s2 = (2*c) / (-b + j);
    }
    std::cout << s1 << std::endl;
    std::cout << s2 << std::endl;
  }
  return 0;
}
