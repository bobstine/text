#include "cross_tab.h"


#include <iostream>

using std::cout;
using std::endl;

int main()
{
  int nR = 4;
  int nC = 3;
  
  CrossTab table(nR,nC);

  for(int i=0; i<nR; ++i)
    for(int j=0; j<nC; ++j)
    { int lim = rand() % 10 + 1;
      for(int k=0; k<lim; ++k)
	table.increment(i,j);
    }
  
  std::cout << "MAIN:  Table printing \n -------------------- \n \n";
  table.print_to_stream(std::cout);
  std::cout << "\n \n -------------------- \n \n";
}

