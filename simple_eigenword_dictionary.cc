#include "simple_eigenword_dictionary.h"

#include <random>
#include <algorithm>

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>   // transform
#include <cmath>       // nan
/*
  ASSUME
    - coordinates for OOV come first
    - dictionary is in inverse Zipf order (so overwrite less with more common)
*/

typedef float Scalar;

const bool verbose = true;

Text::SimpleEigenwordDictionary
Text::make_simple_eigenword_dictionary(std::string filename, size_t dimToUse, Text::SimpleVocabulary const& vocabulary, bool downcase)
{
  using std::string;
  bool needToFlushEigenStream = false;
  SimpleEigenwordDictionary dict;
  {
    std::ifstream eigenStream (filename.c_str(), std::ifstream::in);
    if (eigenStream.fail())
    { std::cerr << "ERROR: Eigen dictionary file `" << filename << "' not found; terminating.\n";
      return dict;
    }
    std::clog << "DICT: Build eigen dictionary with initial dim " << dimToUse << " from file " << filename << std::endl;
    { // special handling for first line, assumed to define OOV coordinates; determine if use a subset of coefs
      std::vector<float> oov;
      string token;
      eigenStream >> token;  // toss file label for OOV; use "OOV"
      string line;
      std::getline(eigenStream, line);
      std::istringstream ss(line);
      float x;
      while (ss >> x) oov.push_back(x);
      if (dimToUse == 0)
      { dimToUse = (int) oov.size();
	if (verbose) std::clog << "Found eigen dimension = " << oov.size() << " based on found size for OOV.\n";
	needToFlushEigenStream = false;
      }
      else
      { if (dimToUse <= oov.size())
	{ oov.resize(dimToUse);
	  needToFlushEigenStream = true;
	  if (verbose && (dimToUse < oov.size())) std::clog << "      Will need to flush trailing elements in eigen stream.\n";
	}
	else
	{ std::cerr << "DICT:  *** ERROR *** Found eigen dimension " << oov.size() << " which is smaller than requested " << dimToUse << std::endl;
	  return dict;
	}
      }
      dict["OOV"] = oov;
    }
    while(!eigenStream.eof())
    { string token;
      string junk;
      eigenStream >> token;
      if(token.size() == 0) break;                                                        // file has trailing blank line
      if (downcase) std::transform(token.begin(), token.end(), token.begin(), ::tolower); // lower case
      if(vocabulary.count(token) != 0)                                                    // word type found in vocabulary
      { std::vector<float> coor(dimToUse);
	for(size_t i=0; i<dimToUse; ++i)
	  eigenStream >> coor[i];
	if (needToFlushEigenStream) std::getline(eigenStream, junk);
	dict[token] = coor;        // over-writes if token present
      }
      else // flush rest of line
	std::getline(eigenStream, junk);
    }
  }
  // propagate missing values; assign nan to missing
  {
    std::vector<float> missing(dimToUse);
    for(size_t i=0; i<dimToUse; ++i)
      missing[i] = std::nanf("missing");
    dict["NA"] = missing;
  }
  return dict;
}


Text::SimpleEigenwordDictionary
Text::make_random_simple_eigenword_dictionary(int seed, size_t dimToUse, Text::SimpleVocabulary const& vocabulary)
{
  std::mt19937 generator;   
  std::normal_distribution<Scalar> normal_dist(0, 1);
  SimpleEigenwordDictionary dict;

  generator.seed(seed); 
  std::vector<Scalar> randomVector (dimToUse);
  std::clog << "DICT: Build *random* simple eigen dictionary of dim " << dimToUse << std::endl;
  // too messy
  // std::generate(randomVector.begin(), randomVector.end(),  [&generator,&normal_dist]()->Scalar { return normal_dist(generator); });
  // warning message about the use of the random generator being uninitialized???
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
  for(size_t i=0; i<dimToUse; ++i)
    randomVector[i] = normal_dist(generator);
  dict["OOV"] = randomVector;   // special handling for first line, assumed to define OOV coordinates
  for(std::string type : vocabulary)
  { for(size_t i=0; i<dimToUse; ++i)
      randomVector[i] = normal_dist(generator);
    dict[type] = randomVector;
  }
#pragma GCC diagnostic pop
  // propagate missing values; assign nan to missing
  {
    std::vector<float> missing(dimToUse);
    for(size_t i=0; i<dimToUse; ++i)
      missing[i] = std::nanf("missing");
    dict["NA"] = missing;
  }
  return dict;
}

void
Text::compare_dictionary_to_vocabulary(Text::SimpleEigenwordDictionary const& dict, Text::SimpleVocabulary const& vocab) 
{
  if (dict.size() < vocab.size())
  { std::clog << "EDIC:  Dictionary coordinates not found for " << vocab.size()-dict.size()
	      << " tokens (some shown below); writing all to `missing_words.txt'.\n";
    std::ofstream os{"missing_words.txt"};
    if(os.good())
    { int counter = 0;
      for(auto x : vocab)
      { if (dict.count(x) == 0)
	{ ++counter;
	  os << " " << x;
	  if (0 == (counter%128))
	  { os << std::endl;
	    std::clog << " " << x;
	    if (0 == (counter%4096))
	      std::clog << std::endl;
	  }
	}
      }
      os << std::endl;
      std::clog << std::endl;
    }
    else std::cerr << "EDIC: Could not open file to write missing words.\n";
  }
}
