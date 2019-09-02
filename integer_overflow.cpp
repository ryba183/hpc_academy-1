#include <iostream>

int main(void)
{
    signed int s = 0;
    for(unsigned int i = 0; i < 10; i++)
    {
      for(unsigned int j = 0; j < 1048576; j++)
      {
        s += 1024;
      }
      std::cout << "s = " << s << std::endl;
    }
    return 0;
}
