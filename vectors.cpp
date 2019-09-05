#include <iostream>
#include <cmath>
#include <complex>
#include <array>
#include <vector>
#include <algorithm>

int main(void) {

  int n;
  std::vector<double> container;
  while(true) {
    std::cout << "Please enter a list of positive numbers, ending with a negative one: " << std::flush;
    std::cin >> n;
    if (n < 0) {
      break;
    }
    container.push_back(n);
  }
  
  std::sort(container.begin(), container.end());
  std::vector<double>::iterator vectorIter;
  std::vector<double> every_other;
  int i = -1;
  for (vectorIter = container.begin(); vectorIter != container.end(); ++vectorIter) {
    i++;
    if (i % 2) continue;
    every_other.push_back(*vectorIter);
  }
  
  std::cout << "Sorted vector, every other element: ";
  for (vectorIter = every_other.begin(); vectorIter != every_other.end(); ++vectorIter) {
   std::cout << *vectorIter << std::endl;
  }
}
