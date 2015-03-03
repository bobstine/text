#include "simple_vocabulary.h"

#include <iostream>
#include <fstream>
#include <algorithm>

#include <assert.h>

const std::string TAG = "SVOC: ";

Text::SimpleVocabulary
Text::make_simple_vocabulary(std::string filename, bool downcase)
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
    if (downcase) std::transform(token.begin(), token.end(), token.begin(), ::tolower); // lower case
    vocab.insert(token);
  }
  std::clog << TAG << "Read vocabulary of " << vocab.size();
  if (downcase)
    std::clog << " unique lower-case tokens from file " << filename << std::endl;
  else
    std::clog << " unique mixed-case tokens from file " << filename << std::endl;
  return vocab;
}
