#include "eigenword_dictionary.h"
#include "file_utils.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <getopt.h>

typedef EigenwordDictionary::Matrix  Matrix;
typedef EigenwordDictionary::Vector  Vector;

void
parse_arguments(int argc, char** argv, string &name);

int main(int argc, char** argv)
{
  // build eigenword dictionary
  const std::string dictionaryFileName ("text_src/eigenwords/google_tri.txt");
  const int nLines = FileUtils::count_lines(dictionaryFileName);
  const int nDictCols = FileUtils::count_fields(dictionaryFileName);
  std::cout << "MAIN: Found " << nDictCols << " columns in file "
	    << dictionaryFileName << " with " << nLines << " lines in eigenword dictionary." << std::endl;
  EigenwordDictionary dict(dictionaryFileName, nLines-1, nDictCols-1);
  std::cout << "MAIN: Dictionary is " << dict << std::endl;

  // read document strings
  std::string baseName = "";
  parse_arguments(argc, argv, baseName);
  const std::string documentFileName ("text_src/anes/" + baseName + ".txt");
  const int nDocs = FileUtils::count_lines(documentFileName);
  std::cout << "MAIN: Found " << nDocs << " lines of text in file " << documentFileName << " for encoding centroids." << std::endl;
  std::vector<std::string> docs;
  Matrix docCentroids(nDocs, dict.eigen_dim());
  std::ifstream is(documentFileName);

  std::string line;
  int lineNumber = 0;
  while (getline(is,line))
  { docs.push_back(line);
    std::vector<std::string> words;
    std::string w;
    for(int i=0; line[i]!= '\0'; ++i)
    { if(line[i]==' ')
      {	if(w.size())
	{ words.push_back(w);
	  w.clear();
	}
      } else
	w.push_back(std::tolower(line[i]));
    }
    if (w.size())
      words.push_back(w);
    
    std::clog << "MAIN: words for line " << lineNumber << " are ";
    for(int i=0; i<(int)words.size(); ++i) std::clog << " `" << words[i] << "'  ";
    std::clog << std::endl;

    Vector centroid = Vector::Zero(docCentroids.cols());
    int nFound = 0;
    for(auto w : words)
    { Vector v = dict[Type(w)];
      if (v.size())
      { ++nFound;
	centroid = centroid + v;
      }
    }
    if (1 < nFound)
      docCentroids.row(lineNumber) = centroid.array() / nFound;
    else
      docCentroids.row(lineNumber) = centroid;    
    ++lineNumber;
  }
  
  // write to output file
  std::clog << "MAIN: Leading rows of centroid matrix are " << std::endl << docCentroids.topRows(5) << std::endl;
  const std::string outFileName ("text_src/anes/" + baseName + "_centroids.txt");
  std::cout << "MAIN: Writing output centroids to file " << outFileName << std::endl;
  std::ofstream of (outFileName);
  of << "Source\t";
  for(int i=0; i<docCentroids.cols(); ++i)
    of << "C" << i << "\t";
  of << std::endl;
  Eigen::IOFormat fmt(Eigen::StreamPrecision,Eigen::DontAlignCols,"\t","\n","","","","");
  for(int i=0; i<docCentroids.rows(); ++i)
    of << docs[i] << "\t" << docCentroids.row(i).format(fmt) << std::endl;
}



void
parse_arguments(int argc, char** argv, string &name)
{
  static struct option long_options[] = {
    {"name",          required_argument, 0, 'n'},
    {0, 0, 0, 0}                             // terminator 
  };
  int key;
  int option_index = 0;
  while (-1 !=(key = getopt_long (argc, argv, "n:", long_options, &option_index))) // colon means has argument
  {
    switch (key)
    {
    case 'n' : { name      = optarg;                                        break; }
    default  : { std::cout << "PARSE: Option not recognized; returning.\n";       }
    } // switch
  } 
}


      //    case 'n' : { n              = read_utils::lexical_cast<int>(optarg);      break; }
