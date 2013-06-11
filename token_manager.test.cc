#include "token_manager.h"

#include <iostream>
#include <fstream>

int main()
{
  // first without condensing small POS
  {
    std::ifstream fs("/Users/bob/C/text/test_in");
    if (!fs)
      std::cout << "TEST: Could not open file.\n";
    else
    { TokenManager tm (fs);
      tm.print_to_stream(std::cout);
    }
  }

  // with reducing small POS to OTH
  {
    std::ifstream fs("/Users/bob/C/text/test_in");
    if (!fs)
      std::cout << "TEST: Could not open file.\n";
    else
    { TokenManager tm (fs, 0.01);     // 1% threshold for combining
      tm.print_to_stream(std::cout);
    }
  }

}
