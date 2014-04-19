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
    {
      while (fs)
      { string str;
	std::cout << " ------------- HERE -------main------" << std::endl;
	fs >> str;
	std::cout << " Read the string... " << str << std::endl;
      }
      std::cout << " ------------- HERE -------main------" << std::endl;

      return 0;
    
      TokenManager tm (fs);
      std::cout << " ------------- HERE -------------" << std::endl;
      std::vector<Type> typeVector = tm.type_vector();   // check indexing
      std::cout << "TEST: first 10 types in type vector from Token Manager with index checks (should match)..." << std::endl;
      for(int i=0; i<10; ++i)
	std::cout << "i=" << i << " type=" << typeVector[i] << "  index=" << tm.type_index(typeVector[i]) << "   frequency=" << tm.type_freq(typeVector[i]) << std::endl;
      tm.print_to_stream(std::cout);
      std::cout << "TEST: filling bigram map\n";
      std::map<std::pair<int,int>,int> bigramMap;
      int skip = 1;
      const bool transpose = false;
      tm.fill_bigram_map(bigramMap, skip, tm, transpose);
    }
  }

  // with reducing small POS to OTH
  {
    std::ifstream fs("/Users/bob/C/text/test_in");
    if (!fs)
      std::cout << "TEST: Could not open file.\n";
    else
    { TokenManager tm (fs, 0.01);                       // 1% threshold for combining
      std::vector<Type> typeVector = tm.type_vector();  // check indexing
      std::cout << "TEST: first 10 types in type vector from Token Manager with index checks (should match)..." << std::endl;
      for(int i=0; i<10; ++i)
	std::cout << "i=" << i << " type=" << typeVector[i] << " type index of this type=" << tm.type_index(typeVector[i]) << std::endl;
      tm.print_to_stream(std::cout);
    }
  }

}
