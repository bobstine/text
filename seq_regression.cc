/*
  This application...

  -- reads a tab delimited data file that begins with the
  names of the variables, then
              nRows of data and nCols of data
  The number of rows does not count the leading line of
  variable names.

  -- fits a sequence of regression models, with Y as the first
  variable read, followed by each of the predictors, then

  -- writes the sequence r2 obtained
  
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
parse_arguments(int argc, char** argv,	int &nRow, int &nCol, string &inputFileName, int &initModelSize, string &outputFileName);

void
fit_models(std::istream& input, int initSize, std::ostream& output);


int main(int argc, char** argv)
{
  // read input options
  int    nRows           ( 0  );
  int    nCols           ( 0  );
  string inputFileName   ( "" ); 
  string outputFileName  ( "" ); 
  int    initModelSize   (  0 );
  parse_arguments(argc, argv, nRows, nCols, inputFileName, initModelSize, outputFileName);
  std::clog << "MAIN: seq_regressor --nrows=" << nRows << " --ncols=" << nCols << " --input_file="
	    << inputFileName << " --init_size=" << initModelSize << " --output_file_name=" << outputFileName << endl;
  std::ofstream output(outputFileName);
  if (!output)
  { std::clog << "MAIN: Could not open file " << outputFileName << " for writing r2 sequence.\n";
    return 0;
  }
  if (inputFileName.size() > 0)
  { std::ifstream input(inputFileName);
    if (!input)
    { std::clot << "MAIN: Could not open file " << inputFileName << " for reading input.\n";
      return 0;
    }
    fit_models(input, nRows, nCols, initModelSize, output);
  }
  else fit_models(std::cin, nRows, nCols, initModelSize, output);
}

void
fit_models(std::istream &input, int n, int nCols, int initModelSize, std::ostream &output)
{
  output << "Name  r2  RSS AICc\n";
  string              yName;
  std::vector<string> iNames;
  std::vector<string> xNames;
  // read variable names
  input >> yName;
  for(int j=0; j<initModelSize; ++j)
  { string name;
    input >> name;
    iNames.push_back(name);
  }
  for(int j=0; j<k; ++j)
  { string name;
    input >> name;
    xNames.push_back(name);
  }
  // read numerical values
  Vector Y(n);
  Matrix Z(n,         initModelSize    );
  Matrix X(n, nCols - initModelSize - 1);
  for(int i = 0; i < n; ++i)
  { input >> Y(i);
    for (int j=0; j<Z.cols(); ++j)
      input >> Z(i,j);
    for (int j=0; j<k; ++j)
      input >> X(i,j);
  }
  // fit initial model
  const int blockSize = 0;
  LinearRegression regr(yName, Y, blockSize);
  if (initModelSize > 0)
  { std::clog << "MAIN: Calculate initial model, fitting " << initModelSize << " regressors to response " << yName << endl;
    FStatistic f=regr.f_test_predictors(iNames, Z);
    if(f.f_stat() > 0.0001) regr.add_predictors();
    else std::clog << "MAIN: F = " << f.f_stat() << " initial regressors (near) singular and skipped.\n";
  }
  std::clog << "MAIN: Calculation loop; fitting " << k << " regressors to response " << yName << endl;
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
parse_arguments(int argc, char** argv, int &rows, int&cols, string &inputFile, int &initModelSize, string &outputFile)
{
  static struct option long_options[] = {
    {"nrows",         required_argument, 0, 'r'},
    {"ncols",         required_argument, 0, 'c'},
    {"input_file",    required_argument, 0, 'i'},
    {"init_size",     required_argument, 0, 'k'},
    {"output_file",   required_argument, 0, 'o'},
    {0, 0, 0, 0}                             // terminator 
  };
  int key;
  int option_index = 0;
  while (-1 !=(key = getopt_long (argc, argv, "r:c:i:k:o:", long_options, &option_index))) // colon means has argument
  {
    switch (key)
    {
    case 'r' : { rows           = read_utils::lexical_cast<int>(optarg);      break; }
    case 'c' : { cols           = read_utils::lexical_cast<int>(optarg);      break; }
    case 'i' : { inputFile      = optarg;                                     break; }
    case 'k' : { initModelSize  = read_utils::lexical_cast<int>(optarg);      break; }
    case 'o' : { outputFile     = optarg;                                     break; }
    default  : { std::cout << "PARSE: Option not recognized; returning.\n";       }
    } // switch
  } 
}
  
