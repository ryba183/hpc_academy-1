#include <iostream>
#include <cmath>
#include <fstream>

int main(void) {
  
  int a;
  do {
    std::cout << "Enter a positive integer: " << std::flush;
    std::cin >> a;
  } while (a < 0);
  
  std::ofstream output;
  output.open("test.txt");
  
  for (int i = 0; i <= a ; i++) {
    output << i << "\t" << i*i << std::endl;
  }
  output.close();
}
