#ifndef _SIMPLE_EIGENWORD_DICTIONARY_H_
#define _SIMPLE_EIGENWORD_DICTIONARY_H_

// An eigen dictionary is a map read from a file

#include "simple_vocabulary.h"

#include <string>
#include <map>
#include <vector>

namespace Text {
  typedef  std::map<std::string, std::vector<float>>  SimpleEigenwordDictionary;
  
  SimpleEigenwordDictionary
    make_simple_eigenword_dictionary(std::string filename, size_t dimToUse, SimpleVocabulary const& vocabulary, bool downcase);

  SimpleEigenwordDictionary
    make_random_simple_eigenword_dictionary(int seed, size_t dimToUse, SimpleVocabulary const& vocabulary);

  void compare_dictionary_to_vocabulary(SimpleEigenwordDictionary const& eigenDict, SimpleVocabulary const& vocab);
}

#endif
