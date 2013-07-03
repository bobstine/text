#include "read_utils.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <getopt.h>

using std::string;
using std::endl;

class TypeIndex  // alas cannot inherit from int, so must have one
{
  int mInt;
  
 public:
  explicit TypeIndex(int i) : mInt(i)       { }

  operator int() const { return mInt; }
  
  void  print_to_stream(std::ostream &os) const   { os << "[" << mInt << "]" ; }
 };


void
parse_arguments(int argc, char** argv, string &vFileName, int &minFrequency);


  
int main(int argc, char** argv)
{
  // read input options
  string inputFileName ( "");
  int    minFrequency  ( 3 );
  parse_arguments(argc, argv, inputFileName, minFrequency);
  std::clog << "MAIN: regressor --input_file " << inputFileName << " --min_frequency " << minFrequency << endl;
  
  // build the initial vocabulary then trim oov
  std::map<string,int> vocabulary;
  {
    std::map<string,int> fullVocabulary;
    std::ifstream inputFile (inputFileName);
    if(!inputFile)
    { std::cerr << "MAIN: input file " << inputFileName << " not found; exiting." << endl;
      return -1;
    }
    string token;
    int nTokens = 0;
    while (inputFile >> token)
    { ++nTokens;
      ++fullVocabulary[token];
    }
    std::clog << "MAIN: Obtain full vocabulary of " << fullVocabulary.size() << " types from input of " << nTokens << " tokens." << endl;
    string oov("OOV");
    std::clog << "MAIN: OOV threshold at frequency " << minFrequency << " tokens." << endl;
    for (auto x : fullVocabulary)
    { if (x.second < minFrequency)  // mark as oov
	vocabulary[oov] += x.second;
      else
	vocabulary.insert(x);       // transfer
    }
    int checkSum (0);
    for (auto x: vocabulary)
      checkSum += x.second;
    std::clog << "MAIN: Thresholded vocabulary of " << vocabulary.size() << " types with token count " << checkSum << " tokens." << endl;    
  }

    
  return 0;
}


void
parse_arguments(int argc, char** argv,
		string &fileName, int &oovThreshold)
{
  static struct option long_options[] = {
    {"input_file",   required_argument, 0, 'i'},
    {"min_frequency",required_argument, 0, 'f'},
    {0, 0, 0, 0}                             // terminator 
  };
  int key;
  int option_index = 0;
  while (-1 !=(key = getopt_long (argc, argv, "i:f:", long_options, &option_index))) // colon means has argument
  {
    // std::cout << "Option key " << char(key) << " for option " << long_options[option_index].name << ", option_index=" << option_index << std::endl;
    switch (key)
    {
    case 'i' : { fileName       = optarg;                                  break; }
    case 'f' : { oovThreshold   = read_utils::lexical_cast<int>(optarg);   break; }
    default  : { std::cout << "PARSE: Option not recognized; returning.\n";       }
    } // switch
  } 
}
