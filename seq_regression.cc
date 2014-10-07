/*
  This application...

  -- reads 3 tab delimited data files that begin with the names of the
     variables on the first line. Following lines have the indicated
     (on the command line) number of data rows and columns.  The
     number of rows does not count the leading line of variable names.

  -- fits a sequence of regression models, with Y as the first
     variable read, followed the initializers, then by each of the
     predictors in X.

  -- writes the sequence of r2, AICc, and optionally obtained by the
     Xs to the output file

  -- random seed controls the splits used in CV.  Setting the seed to
     0 uses a deterministic split (see code in eigen/regression.cc)
     along with output details of fitted models.
  
*/

#include "read_utils.h"
#include "file_utils.h"
#include "regression.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <getopt.h>

#include <Eigen/Core>
#include <Eigen/QR>


#define MIN(A,B) (((A)<(B)) ? (A) : (B))

using std::string;
using std::endl;

typedef Eigen::VectorXd Vector;
typedef Eigen::MatrixXd Matrix;


void
parse_arguments(int argc, char** argv,	int &n, string &yFileName,
		int &ni, string &iFileName, int &nx, string &xFileName,
		int &nFoldsCV, int&randomSeed, bool &orthoRotate, string &outputFileName);

void
fit_models(int  n, std::istream &yStream,
	   int ni, std::istream &iStream,
	   int nx, std::istream &xStream,
	   int randomSeed, int nFoldsCV, bool rotate, std::ostream& output);


int main(int argc, char** argv)
{
  // prepare eigen for multiple threads
  Eigen::initParallel();

  // read input options
  int    n         ( 0  );
  int    ni        ( 0  );       // number X's used to initialize the model (0 means none)
  int    nx        ( 0  );       // number X's used to extend the model (get AIC, CVSS for these)
  int    cvFolds   ( 0  );       // number folds for cross-validation (0 means no cross validation)
  int    randomSeed(2837);       // random seed controls validation slices
  string yFileName ( "" );       // two columns, y is first
  string iFileName ( "" );       // only used if ni>0, preconditioning variables
  string xFileName ( "" );       // sequence to explore in CV
  bool   rotate    (false);
  string outputFileName("");

  parse_arguments(argc, argv, n, yFileName, ni, iFileName, nx, xFileName, cvFolds, randomSeed, rotate, outputFileName);
  std::clog << "MAIN: seq_regressor --n=" << n << " --folds=" << cvFolds;
  if (0 < cvFolds) std::clog << " --seed=" << randomSeed;
  if (rotate)      std::clog << " --rotate";
  std::clog << endl
	    << "  Files are        --y_file=" << yFileName << endl
            << "                   --x_file=" << xFileName << " --nx=" << nx << endl
	    << "                   --i_file=" << iFileName << " --ni=" << ni << endl
	    << "                   --output=" << outputFileName << endl;
  // assume all but the iFile must be named
  std::ifstream yStream (yFileName);
  std::ifstream xStream (xFileName);
  std::ofstream output(outputFileName);
  if (0 == (n * nx))
  { std::cerr << "MAIN: Illegal length, n=" << n << "  nx=" << nx << endl;
    return 0;
  }
  if ((!yStream) || (!xStream) || (!output))
  { std::clog << "MAIN: Could not open required files for r2 sequence ("
	      << yFileName << "," << xFileName << "," << outputFileName << ")." << endl;
    return 0;
  }
  if (iFileName.size() > 0)
  { std::ifstream iStream(iFileName);
    if (!iStream)
    { std::cerr << "MAIN: Could not open file " << iFileName << " for reading initialization data.\n";
      return 0;
    }
    fit_models   (n, yStream, ni, iStream , nx, xStream, cvFolds, randomSeed, rotate, output);
  }
  else fit_models(n, yStream, ni, std::cin, nx, xStream, cvFolds, randomSeed, rotate, output);
}

void
fit_models(int  n, std::istream &yStream,                  // ystream has 1 col, with name first
	   int ni, std::istream &iStream,                  // additional optional conditioning variables
	   int nx, std::istream &xStream,                  // stepwise sequence
	   int nFolds, int seed, bool rotate, std::ostream& output)
{
  std::clog << "MAIN: Reading data for y, Xi, and X from source listed source files.\n";
  // read Y data and total word count (put counts into first col of Xi)
  string yName, countName;
  std::vector<string> iNames;
  Vector Y{n};
  yStream >> yName;
  for(int i = 0; i<n; ++i)
    yStream >> Y(i,0); 
  // read preconditioning X data
  Matrix Xi{n,ni};
  for(int j=0; j<ni; ++j)
  { string name;
    iStream >> name;
    iNames.push_back(name);
  }
  for(int i=0; i<n; ++i)
    for(int j=0; j<ni; ++j)
      iStream >> Xi(i,j);
  // read predictors from x file
  std::vector<string> xNames;
  Matrix X(n,nx);
  for(int j=0; j<nx; ++j)
  { string name;
    xStream >> name;
    xNames.push_back(name);
  }
  for(int i=0; i<n; ++i)
    for(int j=0; j<nx; ++j)
      xStream >> X(i,j);
  // following applies an orthogonal rotation to y, Xi, and X prior to regression
  if (rotate)
  { std::clog << "MAIN: Initialize random matrix with dimensions " << n << "x" << n << std::endl;
    Matrix R = Eigen::MatrixXd::Random(n,n);
    std::clog << "MAIN: Computing Householder orthogonal matrix Q from random start." << std::endl;
    Matrix Q = Eigen::HouseholderQR<Matrix>(R).householderQ();
    std::clog << "MAIN: Computed orthgonal matrix with dimensions " << Q.rows() << "x" << Q.cols() << std::endl;
    // center prior to rotations
    double initialSum = Y.sum();
    Y = Y.array() - Y.sum()/n;
    Y = Q * Y;
    Vector mXi = Xi.colwise().sum()/n;
    for(int i=0; i<Xi.rows(); ++i)
      Xi.row(i) -= mXi;
    Xi= Q * Xi;
    Vector mX = X.colwise().sum()/n;
    for(int i=0; i<X.rows(); ++i)
      X.row(i) -= mX;
    X = Q * X;
    std::clog << "Check of column sum of Y dropped  from " << initialSum << " to " << Y.sum() << ".  For Xi, " << Xi.colwise().sum() << std::endl;
    std::clog << "MAIN: Orthogonal rotations completed\n";
  }
  // call code to validate using threads
  Eigen::MatrixXd results(1+X.cols(),4);            // R2, RSS, AICc, CVSS, starting with preconditions
  validate_regression(Y, Xi, X, nFolds, results, seed);
  Eigen::IOFormat fmt(Eigen::StreamPrecision,Eigen::DontAlignCols,"\t","\n","","","","");
  output << "R2\tRSS\tAICc\tCVSS\n" << results.format(fmt);
  std::clog << "MAIN: Regression completed with results written to output stream."<< endl;
}
 
 
void
parse_arguments(int argc, char** argv, int &n, string &yFileName,
		int &ni, string &iFileName, int &nx, string &xFileName, int &nFolds, int& seed, bool& rotate, string &outputFileName)
{
  static struct option long_options[] = {
    {"n",             required_argument, 0, 'n'},
    {"y_file",        required_argument, 0, 'Y'},
    {"n_i",           required_argument, 0, 'i'},
    {"i_file",        required_argument, 0, 'I'},
    {"n_x",           required_argument, 0, 'x'},
    {"x_file",        required_argument, 0, 'X'},
    {"folds",         required_argument, 0, 'v'},
    {"seed",          required_argument, 0, 's'},
    {"rotate",              no_argument, 0, 'r'},
    {"output_file",   required_argument, 0, 'o'},
    {0, 0, 0, 0}                             // terminator 
  };
  int key;
  int option_index = 0;
  while (-1 !=(key = getopt_long (argc, argv, "n:Y:i:I:x:X:v:s:ro:", long_options, &option_index))) // colon means has argument
  {
    switch (key)
    {
    case 'n' : { n              = read_utils::lexical_cast<int>(optarg);      break; }
    case 'Y' : { yFileName      = optarg;                                     break; }
    case 'i' : { ni             = read_utils::lexical_cast<int>(optarg);      break; }
    case 'I' : { iFileName      = optarg;                                     break; } 
    case 'x' : { nx             = read_utils::lexical_cast<int>(optarg);      break; }
    case 'X' : { xFileName      = optarg;                                     break; }
    case 'v' : { nFolds         = read_utils::lexical_cast<int>(optarg);      break; }
    case 's' : { seed           = read_utils::lexical_cast<int>(optarg);      break; }
    case 'r' : { rotate         = true;                                       break; }
    case 'o' : { outputFileName = optarg;                                     break; }
    default  : { std::cout << "PARSE: Option not recognized; returning.\n";       }
    } // switch
  } 
}
  
