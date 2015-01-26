/*
  Converts the base data directory into a new 'analysis' directory that has a
  binary response as well as subsets the data to match cases used in the binary
  response.

*/

#include <dirent.h>
#include <getopt.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include <string>
#include <vector>
#include <set>
#include <numeric>   // accumulate

#include "string_trim.h"

/////

const bool verbose = true;

const std::string tag = "RCDD: ";

/////

std::vector<bool>
write_binary_response(std::string word0, std::string word1, std::string inputDir, std::string outputDir);

std::vector<std::string>
files_in_directory (std::string dir);

int
rewrite_predictor_file (std::string inputFile, std::vector<bool> const& selector, std::string outputFile);

template <class T>
inline
std::ostream& operator<< (std::ostream& os, std::vector<T> vec)   // /t separated
{
  if(vec.size() > 0)
  { os << vec[0];
    for (size_t i=1; i<vec.size(); ++i) os << "\t" << vec[i];
  }
  return os;
}

void
parse_arguments(int argc, char** argv,
		std::string& inputDataDir, std::string& word0, std::string& word1,
		std::string& outputDataDir);

/////   

int main(int argc, char** argv)
{
  using std::string;

  // defaults
  string inputDataDir  ("input_data_dir/");
  string word0         ("");           // all but word1
  string word1         ("word");
  string outputDataDir ("output_data_dir/");
  
  parse_arguments(argc, argv, inputDataDir, word0, word1, outputDataDir);
  if ( inputDataDir[ inputDataDir.size()-1] != '/')  inputDataDir += '/';
  if (outputDataDir[outputDataDir.size()-1] != '/') outputDataDir += '/';

  std::clog << "recode_data --input_dir=" << inputDataDir << " --output_dir=" << outputDataDir
	    << " --word0=" << word0 << " --word1=" << word1 << std::endl;

  
  // selector records which cases match word0 or word1
  std::vector<bool> selector = write_binary_response(word0, word1, inputDataDir, outputDataDir);
  int nObsSelected = std::accumulate(selector.begin(), selector.end(), 0, [](int tot, bool x) { if(x) return tot+1; else return tot; });
  if (verbose) std::clog << tag << "Writing " << nObsSelected << " cases for word pair " << word0 << "-" << word1
			 << " from input dir " << inputDataDir << std::endl;

  { // write the file with the number of observations
    std::ofstream countFile (outputDataDir + "n_obs");
    if (!countFile.good())
    { std::cerr << tag << "Could not open file `n_obs' for the count.\n";
      return -20;
    }
    countFile << nObsSelected << std::endl;
  }
    
  { // process the rest of the files
    std::vector<string> allfilenames = files_in_directory(inputDataDir);
    std::set<string> removeNames;
    removeNames.insert(".");   removeNames.insert("n_obs");    removeNames.insert("Y");
    removeNames.insert("..");  removeNames.insert("index.sh"); removeNames.insert(word0+"_"+word1);
    std::vector<string> filenames;
    for (auto filename : allfilenames)
      if (0==removeNames.count(filename))
	filenames.push_back(filename);
    for (auto filename : filenames)
    { if (verbose) std::clog << "RECODE: Recoding data file " << filename << std::endl;
      int nCasesWritten = rewrite_predictor_file(inputDataDir+filename, selector, outputDataDir + filename);
      if (nCasesWritten != nObsSelected)
      { std::cerr << "ERROR: Number cases written for " << filename << " was " << nCasesWritten
		  << " != " << nObsSelected << std::endl;
	return -11;
      }
    }
    return 0;
  }
}

//     write_binary_response     write_binary_response     write_binary_response     write_binary_response     

//       converts intput text into 0/1, selecting only appropriate cases identified in selector
//       writes n_obs on first line followed by 3 line auction format.

std::vector<bool>
write_binary_response(std::string word0, std::string word1, std::string inputDir, std::string outputDir)
{
  using std::string;
  
  std::vector<bool> selector;
  std::ifstream input (inputDir + "Y");
  if (!input.good())
    { std::cerr << "ERROR: Cannot open input file " << inputDir << " Y text to convert to binary.\n";
      return selector;
    }
  std::ofstream output (outputDir + "Y");
  if (!output.good())
    { std::cerr << "ERROR: Cannot open output file " << outputDir << " for binary response.\n";
      return selector;
    }
  // read 3-line file, echoing
  string line;
  std::getline(input, line);             // handle names
  std::istringstream ss(line);
  string name;
  ss >> name;
  std::getline(input, line);             // attributes
  std::vector<int> binaryY;
  string word;
  int nObs = 0;
  if (word0.size() == 0) // code all but word1 as 0
    while(input.good() && (input >> word))
    { ++nObs;
      selector.push_back(true);
      if (word == word1)
	binaryY.push_back(1);
      else
	binaryY.push_back(0);
    }
  else // code only two selected words
    while(input.good() && (input >> word))
    { if ((word == word0) || (word == word1))
      { ++nObs;
	selector.push_back(true);
	if (word == word1)
	  binaryY.push_back(1);
	else
	  binaryY.push_back(0);
      }
      else selector.push_back(false);
    }
  // write to output
  output << nObs << std::endl;
  output << "Y" << std::endl;
  output << line << " word0 " << word0 << " word1 " << word1 << " name " << name << std::endl;
  output << binaryY << std::endl;
  return selector;
}


//     rewrite_predictor_file     rewrite_predictor_file     rewrite_predictor_file     rewrite_predictor_file
int
rewrite_predictor_file (std::string inputFile, std::vector<bool> const& selector, std::string outputFile)
{
  std::ifstream input (inputFile);
  if (!input.good())
  { std::cerr << "ERROR: Attempting to rewrite predictor; could not open file " << inputFile << std::endl;
    return 0;
  }
  std::ofstream output (outputFile);
  if (!output.good())
  { std::cerr << "ERROR: In rewrite, could not open output file " << outputFile << std::endl;
    return 0;
  }
  std::string line;
  std::getline(input, line);    // copy first two lines
  output << line << std::endl;
  std::getline(input, line);
  output << line << std::endl;
  int count = 0;
  for(int i=0; i<(int)selector.size(); ++i)
  { float x;
    input >> x;
    if(selector[i])
    { output << x << "\t";
      ++count;
    }
  }
  output << std::endl;   
  
  return count;
}
  


//     files_in_directory     files_in_directory     files_in_directory     files_in_directory     files_in_directory
std::vector<std::string>
files_in_directory (std::string dir)
{
  DIR *dp;
  struct dirent *dirp;
  
  std::vector<std::string> files;
  if((dp  = opendir(dir.c_str())) == NULL)
  { std::cerr << "Error(" << errno << ") opening " << dir << std::endl;
    return files;
  }
  while ((dirp = readdir(dp)) != NULL)
    files.push_back(std::string(dirp->d_name));
  closedir(dp);
  if (verbose)
  { std::clog << "     Found the following files: " ;
    for(auto f : files) std::clog << f << ",";
    std::clog << std::endl;
  }
  return files;
}


//     parse_arguments     parse_arguments     parse_arguments     parse_arguments     parse_arguments     parse_arguments       
void
parse_arguments(int argc, char** argv,
		std::string& inputDir, std::string& word0, std::string& word1, std::string& outputDir)
{
  static struct option long_options[] = {
    {"input_dir",  required_argument, 0, 'i'},
    {"output_dir", required_argument, 0, 'o'},
    {"word1",      required_argument, 0, 'w'},
    {"word0",      required_argument, 0, 'z'},
    {0, 0, 0, 0}                             // terminator
  };
  int key;
  int option_index = 0;
  while (-1 !=(key = getopt_long (argc, argv, "i:o:w:z:", long_options, &option_index))) // colon means has argument
  {
    switch (key)
    {
    case 'i' :  { inputDir  = optarg; break;      }
    case 'o' :  { outputDir = optarg; break;      }
    case 'w' :  { word1     = optarg; break;      }
    case 'z' :  { word0     = optarg; break;      }
      //    case 'n' :      { 	nRounds = read_utils::lexical_cast<int>(optarg);	break;      }
    default:
      {
	std::cout << "PARSE: Option not recognized; returning.\n";
      }
    } // switch
  } // while
}

