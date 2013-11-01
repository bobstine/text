/*
  This application reads a data file, then fits a sequence of
  regression models, with Y as the first variable read then each
  of the following variables.
*/

#include "helpers.h"
#include "vocabulary.h"
#include "file_utils.h"

#include <iostream>
#include <fstream>
#include <vector>

#include <Eigen/Core>


#define MIN(A,B) (((A)<(B)) ? (A) : (B))

using std::string;
using std::endl;

typedef Eigen::VectorXd Vector;
typedef Eigen::MatrixXd Matrix;


void
parse_arguments(int argc, char** argv,	string &inputFileName, string &outputFileName);
  
int main(int argc, char** argv)
{
  // read input options
  string inputFileName   (  ""   ); 
  string outputFileName  (  ""   ); 
  
  parse_arguments(argc, argv, inputFileName, outputFileName);
  std::clog << "MAIN: regressor --input_file=" << inputFileName << " --output_file_name=" << endl;

  const int nCols  = count_fields(inputFileName);
  const int nLines = count_lines(inputFileName);

  string              yName;
  std::vector<string> xNames;
  std::ifstream input(inputFileName);
  std::ofstream output(file);
  if ((file.size()>0) && (!output)) std::clog << "MAIN: Could not open file " << file << " for writing r2 sequence.\n";
  os << "Name  r2  RSS AICc\n";

  const int n = nLines-1;
  const int k = nCols -1;
  {
    get_line(input);                       // first line has column names
    std::stringstream ss (input);
    ss >> yName;
    for(int j=0; j<k; ++j)
    { string name;
      ss >> name;
      xNames.push_back(name);
    }
  }
  Vector Y(n);
  Matrix X(n, k);
  for(int i = 0; i < n; ++i)
  { string line;
    get_line(input, line);
    std::stringstream ss(line);
    ss >> Y(i);
    for (int j=0; j<k; ++j)
      ss >> X.coeff(i,j);
  }

  std::clog << "MAIN: Calculate fits loop; fitting " << k << " regressors labeled with Y = " << yName << endl;
  LinearRegression regr(yName, Y, 0);
  for (int j=0; j<k; ++j)
  { FStatistic f=regr.f_test_predictor(xNames[j], X.col(j));
    if(f.f_stat() > 0.0001) regr.add_predictors();
    else std::clog << "MAIN: F = " << f.f_stat() << " for " << xNames[j] << " is (near) singular in sequence r2 and skipped.\n";
    if (os)
      os << xName[j] << " " << regr.r_squared() << " " << regr.residual_ss() << " " << regr.aic_c() << endl;
  }
  std::clog << "MAIN: Regression completed with results written to " << outputFileName << endl;
}


void
parse_arguments(int argc, char** argv, string &inputFile, string &outputFile)
{
  static struct option long_options[] = {
    {"input_file",    required_argument, 0, 'i'},
    {"output_file",   required_argument, 0, 'o'},
    {0, 0, 0, 0}                             // terminator 
  };
  int key;
  int option_index = 0;
  while (-1 !=(key = getopt_long (argc, argv, "i:o:", long_options, &option_index))) // colon means has argument
  {
    switch (key)
    {
    case 'i' : { inputFile   = optarg;                                     break; }
    case 'o' : { outputFile  = optarg;                                     break; }
    default  : { std::cout << "PARSE: Option not recognized; returning.\n";       }
    } // switch
  } 
}
  
