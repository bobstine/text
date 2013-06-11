#include "confusion_matrix.h"

#include <vector>
#include <iostream>


int main()
{
  using std::string;
  
  std::vector<string> uniqueLabels;
  uniqueLabels.push_back("AAA");
  uniqueLabels.push_back("BBB");
  uniqueLabels.push_back("CCC");
  uniqueLabels.push_back("DDD");
  uniqueLabels.push_back("EEE");

  const int      n         =  20;
  const float  matchProb   = 0.75;
  std::vector<string>  trueLabels(n);
  std::vector<string>   estLabels(n);

  srand(2032);
  for(int i=0; i<n; ++i)
  {  int labelIndex = rand() % 5;
     trueLabels[i] = uniqueLabels[labelIndex];
     float m = ((float)(rand() % 1000))/1000;
     if (m < matchProb)
       estLabels[i] = uniqueLabels[labelIndex];
     else
       estLabels[i] = uniqueLabels[rand() % 5];
  }

  /*
    std::cout << "TEST: data labels...\n";
    for(int i=0; i<n; ++i)
    std::cout << trueLabels[i] << "  " << estLabels[i] << std::endl;
  */
  
  ConfusionMatrix mat(trueLabels.cbegin(), trueLabels.cend(), estLabels.cbegin());

  std::cout << "TEST:  Accuracy=" << mat.accuracy() << "  Purity=" << mat.purity() << std::endl;
  mat.print_to_stream(std::cout);  
}
