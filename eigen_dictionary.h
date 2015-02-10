#ifndef _eigen_dictionary_h
#define _eigen_dictionary_h

// An eigen dictionary is a map read from a file

#include "simple_vocabulary.h"

#include <string>
#include <map>
#include <vector>

typedef  std::map<std::string, std::vector<float>>  EigenDictionary;

EigenDictionary
make_eigen_dictionary(std::string filename, size_t dimToUse, SimpleVocabulary const& vocabulary);

void compare_dictionary_to_vocabulary(EigenDictionary const& eigenDict, SimpleVocabulary const& vocab);


#endif
