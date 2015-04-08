#include <iostream>
#include <cmath>
#include "simple_vocabulary.h"
#include "simple_eigenword_dictionary.h"

#include "nan_utils.h"

inline
std::ostream&
operator<<(std::ostream& os, std::vector<float>  const& x)
{
  for(float xi : x)
    os << " " << xi;
  return os;
}

int main()
{
  Text::SimpleVocabulary vocab;
  vocab.insert("word1");
  vocab.insert("word2");
  vocab.insert("word3");
  vocab.insert("word4");
  vocab.insert("word5");
  vocab.insert("word6");

  const int    seed = 1234;
  const size_t  dim = 5;
  Text::SimpleEigenwordDictionary dict = Text::make_random_simple_eigenword_dictionary(seed, dim, vocab);

  std::cout << "Dict[OOV] = "   << dict["OOV"]   << std::endl;

  std::cout << "Dict[word1] = " << dict["word1"] << std::endl;
  std::cout << "Dict[word2] = " << dict["word2"] << std::endl;

  std::vector<float> coordNA = dict["NA"];
  std::cout << "Dict[NA] = "    << coordNA << std::endl;
  if (std::isnan(coordNA[1]))
    std::cout << "Std Test: NA coor = NaN\n";
  else
    std::cout << "Std Test: NA coor " << coordNA[1] << " is *not* equal to Nan\n";
  if (IsNan(coordNA[1]))
    std::cout << "Utils Test: NA coor = NaN\n";
  else
    std::cout << "Utils Test: NA coor " << coordNA[1] << " is *not* equal to Nan\n";

  return 0;
}
  
