#include <iostream>

int main(void)
{
  int a, b, c, d;
  std::cout << "Enter numerator:" << std::flush;
  std::cin >> a;
  std::cout << "Enter denominator:" << std::flush;
  std::cin >> b;
  
  c = a/b;
  std::cout << 
  a << " divided by " << b << " is " << c 
  << std::endl;

  d = (a/b)*b + a%b;
  if(d == a) 
  {
    std::cout << "Valid" << std::endl;
  }
  return 0;
}
