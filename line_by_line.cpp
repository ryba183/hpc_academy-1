#include <iostream>
#include <string>
#include <sstream>

int main() {
  for (std::string line; std::getline(std::cin, line);) {
      std::string header;
      if (line.rfind(">", 0) == 0) {
        header = line;
        //std::stringstream ss(line);
      } else {
        std::string token;
        std::stringstream ss(header);
        for(char& c : line) {
          while (getline(ss,token, ';')) {
            std::cout << token; 
          }
          std::cout << c;
        }
      }
  }
  return 0;
}
