#include <iostream>
#include <cmath>
#include <complex>
#include <array>
#include <vector>
int main(void) {
  std::array<std::array<int, 3>, 3> matrix;
  int v;
  // Matrix filled from left to right??
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      std::cout << "Enter m(" << i << "," << j << "): " << std::flush;
      std::cin >> v;
      matrix[i][j] = v;
    }
  }
  int g = matrix[1][0];
  std::cout << g;
}
