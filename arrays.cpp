#include <iostream>
#include <cmath>
#include <complex>
#include <array>

void partial_sum(std::array<int, 20> input,
                 std::array<int, 20> output) {
  int p = 0;
  for(size_t i = 0; i < input.size(); i++) {
    p += input[i];
    output[i] = p;
    std::cout << output[i] << std::endl;
  }
}

int main(void) {
  std::array<int, 20>
    input{1,2,3,4,5,6,7,8,9,10,
      11,12,13,14,15,16,17,18,19,20};
  std::array<int, 20> output;
  partial_sum(input, output);
}
