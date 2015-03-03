#ifndef _simple_eigen_dictionary_h
#define _simple_eigen_dictionary_h

// An eigen dictionary is a map read from a file

#include "simple_vocabulary.h"

#include <string>
#include <map>
#include <vector>

namespace Text {
  typedef  std::map<std::string, std::vector<float>>  SimpleEigenDictionary;
  
  SimpleEigenDictionary
    make_simple_eigen_dictionary(std::string filename, size_t dimToUse, SimpleVocabulary const& vocabulary, bool downcase);
  
  void compare_dictionary_to_vocabulary(SimpleEigenDictionary const& eigenDict, SimpleVocabulary const& vocab);
}

#endif
