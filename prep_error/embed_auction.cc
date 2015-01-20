/*
  Embeds the output of convert.cc into an eigenword representation
  Builds the eigenword dictionary with prior input from a given
  vocabulary.

  Currently assuming no POS in the text

  Input: Rectangular style (tab delimited, with uniform set of columns)
       Y        BGL      BGR   FF   FH         PV
       into     hurry    the   ,    classroom  [blank]
       around   look     the   IN   [blank]    look

  Output: embeds explanatory features into eigenwords in streaming layout
          wanted by auction

	      response
	      cv_indicator
	      eigen_stream_1
	      eigen_stream_2
	       ...


*/

#include "read_utils.h"
#include <getopt.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include <string>
#include <vector>
#include <map>
#include <set>

#include <algorithm>
#include <iterator>


void
parse_arguments(int argc, char** argv, bool &placeholder,
		std::string& vocabularyFileName,
		std::string& eigenwordFileName, int& eigenwordDimension,
		std::string& outputDirectory);

template <class T>
inline
std::ostream& operator<< (std::ostream& os, std::vector<T> vec)
// write with tabs *between* elements
{
  // would be nice, but appends a trailing separator
  //   std::copy(vec.begin(), vec.end(), std::ostream_iterator<T>(os,"\t"));
  if(vec.size() > 0)
  { os << vec[0];
    for (size_t i=1; i<vec.size(); ++i) os << "\t" << vec[i];
  }
  return os;
}

///

const bool verbose = false;

///

int main(int argc, char** argv)
{
  using std::string;
  using std::vector;
  using std::set;
  using std::map;

  // defaults
  string vocabFileName ("vocabulary.txt");
  string eigenFileName ("eigenwords.test");
  int    nEigenDim     (0);                 // use all that are found
  string outputDir     ("data_dir");
  parse_arguments(argc, argv, vocabFileName, eigenFileName, nEigenDim, outputDir);

  // read vocabulary
  set<string> vocabulary;
  {
    std::ifstream vocabStream (vocabFileName.c_str(), std::ifstream::in);
    if (vocabStream.fail())
    { std::cerr << "ERROR: Vocabulary file `" << vocabFileName << "' not found; terminating.\n";
      return 0;
    }
    while (! vocabStream.eof() )
    { string token;
      vocabStream >> token;
      vocabulary.insert(token);
    }
    std::clog << "MAIN: Read vocabulary of " << vocabulary.size() << " tokens\n";
  }

  // build dictionary for this vocab from eigenwords  (Paramveer's is space delimited, *UNKNOWN* first)
  // ASSUME coordinates for OOV come first
  map<string, vector<float>> dictionary;
  bool needToFlushEigenStream = false;
  {
    std::ifstream eigenStream (eigenFileName.c_str(), std::ifstream::in);
    if (eigenStream.fail())
    { std::cerr << "ERROR: Eigen dictionary file `" << eigenFileName << "' not found; terminating.\n";
      return 0;
    }
    std::clog << "MAIN: Starting to build dictionary with initial dim " << nEigenDim << " from file " << eigenFileName << std::endl;
    { // special handling for first line, assumed to define OOV coordinates
      vector<float> oov;
      string token;
      eigenStream >> token;  // toss file label for OOV; use "OOV"
      string line;
      std::getline(eigenStream, line);
      std::istringstream ss(line);
      float x;
      while (ss >> x) oov.push_back(x);
      if (nEigenDim == 0)
      { nEigenDim = (int) oov.size();
	if (verbose) std::clog << "Found eigen dimension = " << oov.size() << " based on found size for OOV.\n";
	needToFlushEigenStream = false;
      }
      else
      { if (nEigenDim < (int)oov.size())
	{ oov.resize(nEigenDim);
	  needToFlushEigenStream = true;
	  if (verbose) std::clog << "       Will need to flush trailing elements in eigen stream.\n";
	}
	else
	{ std::cerr << "ERROR: Found eigen dimension " << oov.size() << " which is smaller than requested " << nEigenDim << std::endl;
	  return 0;
	}
      }
      dictionary["OOV"] = oov;
    }
    while(!eigenStream.eof())
    { string token;
      string junk;
      eigenStream >> token;
      if(token.size() == 0) break;       // file has trailing blank line
      std::transform(token.begin(), token.end(), token.begin(), ::tolower);
      if(vocabulary.count(token) != 0)   // word type found in vocabulary
      { std::vector<float> coor(nEigenDim);
	for(int i=0; i<nEigenDim; ++i)
	  eigenStream >> coor[i];
	if (needToFlushEigenStream) std::getline(eigenStream, junk);
	dictionary[token] = coor;        // over-writes if token present
	if (verbose) std::clog << "      Added dictionary coor for token " << token << std::endl;
      }
      else // flush rest of line
	std::getline(eigenStream, junk);
    }
  }
  std::clog << "      Completed eigenword dictionary of size " << dictionary.size() << std::endl;
  if (dictionary.size() < vocabulary.size())
  { std::clog << "      Dictionary coordinates not found for the following " << vocabulary.size()-dictionary.size() << " tokens:\n";
    for(auto x : vocabulary)
    { if (dictionary.count(x) == 0)
	std::clog << " " << x;
    }
    std::clog << std::endl;
  }

  // process each *column* to std output; use header line to count fields
  string responseName;
  vector<string> bundleNames;             
  if (!std::cin.eof())
  { std::cin >> responseName;       
    string headerLine;
    std::getline(std::cin, headerLine);
    std::istringstream ss(headerLine);
    string name
    while (ss >> name)
      bundleNames.push_back(name);
  }
  std::clog << "MAIN: Reading response and " << bundleNames.size() << " words to define eigen coordinates.\n";
  std::vector<             string>   response;
  std::vector< std::vector<string> > theWords(bundleNames.size());
  while (!std::cin.eof())
  { string thePrep;
    std::cin >> thePrep;
    if (thePrep.size() == 0) break;
    response.push_back(thePrep);
    for(int i=0; i<(int)bundleNames.size(); ++i)
    { string token;
      std::cin >> token;
      theWords[i].push_back(token);
    }
  }
  std::clog << "MAIN: Read " << response.size() << " cases for response and " << theWords[0].size() << " words for first predictor.\n";
  std::clog << "MAIN: Embedding " << bundleNames.size() << " blocks of eigen coordinates.\n";
  {
    std::ofstream file (outputDir + "response");
    file << responseName << std::endl;
    file << "role y" << std::endl;
    file << response << std::endl;
  }
  {
    int n = (int)response.size();
    for(int bundle=0; bundle<(int)bundleNames.size(); ++i)
    { std::vector<std::ofstream&> files;
      for(i=0; i<nEigenDim; ++i)
      { std::ofstream file (outputDir + bundleNames[bundle] + "_" + std::to_string(i));
	files.push_back(file);
      }
      for (int i=0; i<n; ++i)
      { string token = theWords[bundle][i];
	if (dictionary.count(token) == 0)
	{ if (verbose) std::clog << "WARNING: Token " << token << " was not found. Treating as OOV.\n";
	  token = "OOV";
	}
	std::vector<float> coord = dictionary[token];
	if(i < n-1)
	  for(j=0; j<nEigenDim; ++j) files[j] << coord[j] << "\t";
	else // avoid trailing tab
	  for(j=0; j<nEigenDim; ++j) files[j] << coord[j];
      }
    }
  }
}

void
parse_arguments(int argc, char** argv,
		std::string& vocabFileName, std::string& eigenFileName, int& eigenDim, std::string outputDir)
{
  static struct option long_options[] = {
    {"eigen_dim",  required_argument, 0, 'd'},
    {"eigen_file", required_argument, 0, 'e'},
    {"output_dir", required_argument, 0, 'o'},
    {"keep_POS",   no_argument,       0, 'p'},
    {"vocab",      required_argument, 0, 'v'},
    {0, 0, 0, 0}                             // terminator
  };
  int key;
  int option_index = 0;
  while (-1 !=(key = getopt_long (argc, argv, "d:e:o:pv:", long_options, &option_index))) // colon means has argument
  {
    // std::cout << "PARSE: Key " << char(key) << " for option " << long_options[option_index].name << ", option_index=" << option_index << std::endl;
    switch (key)
    {
    case 'd' :  { eigenDim = read_utils::lexical_cast<int>(optarg); break; }
    case 'e' :  { eigenFileName = optarg; break;      }
    case 'o' :  { outputDir = optarg;     break;      }
    case 'p' :  { keepPOS = true;         break;      }
    case 'v' :  { vocabFileName = optarg; break;      }
    default:
      {
	std::cout << "PARSE: Option not recognized; returning.\n";
      }
    } // switch
  } // while
}

