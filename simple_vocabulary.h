#ifndef _simple_vocabulary_h
#define _simple_vocabulary_h

// Simple vocabulary is just a set of unique words constructed from file
#include <set>
#include <string>

typedef std::set<std::string> SimpleVocabulary;


SimpleVocabulary
make_simple_vocabulary(std::string filename);


#endif
