/*
  This application...

  -- reads 3 tab delimited data files that begin with the names of the
  variables, then have the indicated number of data rows and columns.
  The number of rows does not count the leading line of variable
  names.

  -- fits a sequence of regression models, with Y as the first
  variable read, followed the initializers, then by each of the
  predictors in X.

  -- writes the sequence r2 obtained by the Xs
  
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
parse_arguments(int argc, char** argv,	int &n,
		int &ny, string &yFileName,
		int &ni, string &iFileName,
		int &nx, string &xFileName, string &outputFileName);

void
fit_models(int n, int ny, std::istream& yStream,
	   int ni, std::istream &iStream,
	   int nx, std::istream &xStream, std::ostream& output);


int main(int argc, char** argv)
{
  // read input options
  int    n         ( 0  );
  int    ny        ( 2  );    // Y, m
  int    ni        ( 0  );
  int    nx        ( 0  );
  string yFileName ( "" ); 
  string iFileName ( "" );
  string xFileName ( "" );
  string outputFileName("");

  parse_arguments(argc, argv, n, ny, yFileName, ni, iFileName, nx, xFileName, outputFileName);
  std::clog << "MAIN: seq_regressor --n=" << n
	    << " --ny=" << ny << " --y_file=" << yFileName 
	    << " --ni=" << ni << " --i_file=" << iFileName
            << " --nx=" << nx << " --x_file=" << xFileName
	    << " --output_file_name=" << outputFileName << endl;
  // assume all but the iFile must be named
  std::ifstream yStream (yFileName);
  std::ifstream xStream (xFileName);
  std::ofstream output(outputFileName);
  if (0 == (n * ny * ni * nx))
  { std::cerr << "MAIN: Illegal length, n=" << n << "  ny=" << ny << "  ni=" << ni << "  nx=" << nx << endl;
    return 0;
  }
  if ((!yStream) || (!xStream) || (!output))
  { std::clog << "MAIN: Could not open required file for writing r2 sequence.\n";
    return 0;
  }
  if (iFileName.size() > 0)
  { std::ifstream iStream(iFileName);
    if (!iStream)
    { std::clot << "MAIN: Could not open file " << iFileName << " for reading initialization data.\n";
      return 0;
    }
    fit_models(n, ny, yStream, ni, iStream , nx, xStream, output);
  }
  else fit_models(n, ny, yStream, ni, std::cin, nx, xStream, output);
}

void
fit_models(int n, int ny, std::istream& yStream,
	   int ni, std::istream &iStream,
	   int nx, std::istream &xStream, std::ostream& output);
{
  output << "Name  r2  RSS AICc\n";
  // read Y data
  std::vector<string> yNames;
  Matrix Y(n,ny);
  for(int j=0; j<ny; ++j)
  { string name;
    yStream >> name;
    yNames.push_back(name);
  }
  for(int i = 0; i<n; ++i)
    for(int j=0; j<ny; ++j)
      yStream >> Y(i,j);
  // read initialization data
  std::vector<string> iNames;
  Matrix Xi(n,ni);
  for(int j=0; j<ni; ++j)
  { string name;
    iStream >> name;
    iNames.push_back(name);
  }
  for(int i=0; i<n; ++i)
    for(int j=0; j<ni; ++j)
      iStream >> Xi(i,j);
  // read sequence of predictors
  std::vector<string> xNames;
  Matrix X(n,ni);
  for(int j=0; j<nx; ++j)
  { string name;
    xStream >> name;
    xNames.push_back(name);
  }
  for(int i=0; i<n; ++i)
    for(int j=0; j<nx; ++j)
      xStream >> X(i,j);
  // fit initial model
  const int blockSize = 0;
  LinearRegression regr(yName, Y.col(0), blockSize);
  if (initModelSize > 0)
  { std::clog << "MAIN: Calculate initial model, fitting " << 1+ni << " regressors to response " << yName << endl;
    FStatistic f = regr.f_test_predictor("nTokens",Y.col(1));
    if(f.f_stat() > 0.0001) regr.add_predictors();
    f = regr.f_test_predictors(iNames, Xi);
    if(f.f_stat() > 0.0001) regr.add_predictors();
    else std::clog << "MAIN: F = " << f.f_stat() << " initial regressors (near) singular and skipped.\n";
  }
  std::clog << "MAIN: Calculation loop; fitting " << nx << " regressors to response " << yName << endl;
  for (int j=0; j<nx; ++j)
  { FStatistic f=regr.f_test_predictor(xNames[j], X.col(j));
    if(f.f_stat() > 0.0001) regr.add_predictors();
    else std::clog << "MAIN: F = " << f.f_stat() << " for " << xNames[j] << " is (near) singular in sequence r2 and skipped.\n";
    if (output)
      output << xNames[j] << " " << regr.r_squared() << " " << regr.residual_ss() << " " << regr.aic_c() << endl;
  }
  std::clog << "MAIN: Regression completed with results written to " << outputFileName << endl;
}


void
parse_arguments(int argc, char** argv, int &n,
		int &ny, string &yFileName,
		int &ni, string &iFileName,
		int &nx, string &xFileName, string &outputFileName)
{
  static struct option long_options[] = {
    {"n",             required_argument, 0, 'n'},
    {"n_y",           required_argument, 0, 'y'},
    {"y_file",        required_argument, 0, 'Y'},
    {"n_i",           required_argument, 0, 'i'},
    {"i_file",        required_argument, 0, 'I'},
    {"n_x",           required_argument, 0, 'x'},
    {"x_file",        required_argument, 0, 'X'},
    {"output_file",   required_argument, 0, 'o'},
    {0, 0, 0, 0}                             // terminator 
  };
  int key;
  int option_index = 0;
  while (-1 !=(key = getopt_long (argc, argv, "n:y:Y:i:I:x:X:o:", long_options, &option_index))) // colon means has argument
  {
    switch (key)
    {
    case 'n' : { n           = read_utils::lexical_cast<int>(optarg);      break; }
    case 'y' : { ny          = read_utils::lexical_cast<int>(optarg);      break; }
    case 'Y' : { yFileName   = optarg;                                     break; }
    case 'i' : { ni          = read_utils::lexical_cast<int>(optarg);      break; }
    case 'I' : { iFileName   = optarg;                                     break; }
    case 'x' : { nx          = read_utils::lexical_cast<int>(optarg);      break; }
    case 'X' : { xFileName   = optarg;                                     break; }
    case 'o' : { outputFile     = optarg;                                     break; }
    default  : { std::cout << "PARSE: Option not recognized; returning.\n";       }
    } // switch
  } 
}
  
