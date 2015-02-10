#include "simple_vocabulary.h"

#include <iostream>
#include <fstream>
#include <assert.h>

const std::string TAG = "SVOC: ";

SimpleVocabulary
make_simple_vocabulary(std::string filename)
{
  SimpleVocabulary vocab;
  std::ifstream vocabStream (filename.c_str(), std::ifstream::in);
  if (vocabStream.fail())
  { std::cerr << TAG << " *** ERROR ***   Vocabulary file `" << filename << "' not found; terminating.\n";
    assert(false);
  }
  while (vocabStream.good())
  { std::string token;
    vocabStream >> token;
    vocab.insert(token);
  }
  std::clog << TAG << "Read vocabulary of " << vocab.size() << " tokens from file " << filename << std::endl;
  return vocab;
}
