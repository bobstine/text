#include "eigen_dictionary.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>   // transform
#include <cmath>       // nanf
/*
  ASSUME
    - coordinates for OOV come first
    - dictionary is in inverse Zipf order (so overwrite less with more common)
*/

const bool verbose = true;

EigenDictionary
make_eigen_dictionary(std::string filename, size_t dimToUse, SimpleVocabulary const& vocabulary)
{
  using std::string;
  bool needToFlushEigenStream = false;
  EigenDictionary dict;
  {
    std::ifstream eigenStream (filename.c_str(), std::ifstream::in);
    if (eigenStream.fail())
    { std::cerr << "ERROR: Eigen dictionary file `" << filename << "' not found; terminating.\n";
      return dict;
    }
    std::clog << "DICT: Build eigen dictionary with initial dim " << dimToUse << " from file " << filename << std::endl;
    { // special handling for first line, assumed to define OOV coordinates
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
      { if (dimToUse < oov.size())
	{ oov.resize(dimToUse);
	  needToFlushEigenStream = true;
	  if (verbose) std::clog << "      Will need to flush trailing elements in eigen stream.\n";
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
      if(token.size() == 0) break;       // file has trailing blank line
      std::transform(token.begin(), token.end(), token.begin(), ::tolower);
      if(vocabulary.count(token) != 0)   // word type found in vocabulary
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
    std::vector<float> missing;
    for(size_t i=0; i<dimToUse; ++i)
      missing.push_back(std::nanf("missing"));
    dict["NA"] = missing;
  }
  return dict;
}



void
compare_dictionary_to_vocabulary(EigenDictionary const& dict, SimpleVocabulary const& vocab) 
{
  if (dict.size() < vocab.size())
  { std::clog << "EDIC:  Dictionary coordinates not found for the following " << vocab.size()-dict.size() << " tokens:\n";
    for(auto x : vocab)
    { if (dict.count(x) == 0)
	std::clog << " " << x;
    }
    std::clog << std::endl;
  }
}
