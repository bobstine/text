/*
  This application reads a data file, then fits a sequence of
  regression models, with Y as the first variable read then each
  of the following variables.
*/

#include "file_utils.h"
#include "regression.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <getopt.h>

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
  std::clog << "MAIN: seq_regressor --input_file=" << inputFileName << " --output_file_name=" << outputFileName << endl;

  const int nCols  = FileUtils::count_fields(inputFileName);
  const int nLines = FileUtils::count_lines(inputFileName);

  string              yName;
  std::vector<string> xNames;
  std::ifstream input(inputFileName);
  std::ofstream output(outputFileName);
  if ((input) && (!output))
  { std::clog << "MAIN: Could not open file " << outputFileName << " for writing r2 sequence.\n";
    return 0;
  }
  output << "Name  r2  RSS AICc\n";

  const int n = nLines-1;
  const int k = nCols -1;
  // read variable names
  input >> yName;
  for(int j=0; j<k; ++j)
  { string name;
    input >> name;
    xNames.push_back(name);
  }
  // read numerical values
  Vector Y(n);
  Matrix X(n, k);
  for(int i = 0; i < n; ++i)
  { input >> Y(i);
    for (int j=0; j<k; ++j)
    { double x;
      input >> x;
      X(i,j) = x;
    }
  }
  // fit models
  std::clog << "MAIN: Calculate fits loop; fitting " << k << " regressors labeled with Y = " << yName << endl;
  LinearRegression regr(yName, Y, 0);
  for (int j=0; j<k; ++j)
  { FStatistic f=regr.f_test_predictor(xNames[j], X.col(j));
    if(f.f_stat() > 0.0001) regr.add_predictors();
    else std::clog << "MAIN: F = " << f.f_stat() << " for " << xNames[j] << " is (near) singular in sequence r2 and skipped.\n";
    if (output)
      output << xNames[j] << " " << regr.r_squared() << " " << regr.residual_ss() << " " << regr.aic_c() << endl;
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
  
