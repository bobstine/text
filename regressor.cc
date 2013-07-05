#include "vocabulary.h"
#include "read_utils.h"
#include "file_utils.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <getopt.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SVD>

using std::string;
using std::endl;

typedef Eigen::VectorXf Vector;
typedef Eigen::MatrixXf Matrix;


void
parse_arguments(int argc, char** argv,
		string &vFileName, int &minFrequency, bool &bidirectional, int &nProjections, string &rFileName);

  
int main(int argc, char** argv)
{
  // read input options
  string vocabFileName (  ""  );  // text used to build bigram, eigenwords
  string regrFileName  (  ""  );  // modeling text, with leading y followed by text
  int    minFrequency  (   3  );
  bool   bidirectional ( true );
  int    nProjections  (  50  );
  parse_arguments(argc, argv, vocabFileName, minFrequency, bidirectional, nProjections, regrFileName);
  std::clog << "MAIN: regressor --vocab_file=" << vocabFileName << " --regr_file=" << regrFileName
	    << " --min_frequency=" << minFrequency << " --n_projections=" << nProjections;
  if (bidirectional) std::clog << "--bidirectional ";
  std::clog << endl;
  
  // build vocabulary from source
  Vocabulary vocabulary(vocabFileName, minFrequency);
  std::clog << "MAIN: " << vocabulary << endl;

  // compute bigram matrix
  Vocabulary::SparseMatrix B(vocabulary.n_types(), vocabulary.n_types());
  vocabulary.fill_sparse_bigram_matrix(B,0);
  std::clog << "MAIN: Init sparse bigram B[" << B.rows() << "x" << B.cols() << "] from map; sum +/,B= " << B.sum() << endl;
  std::clog << "MAIN: Bigram row sums are (0)"  // these should check with those printed by the vocabulary (but for 1 pairs, not singles)
	    << B.row(0).sum() << "  (1)" << B.row(1).sum() << "  (2)"
	    << B.row(2).sum() << "  (3)" << B.row(3).sum() << "  (4)" << B.row(4).sum() << std::endl;

  // form random projections
  Matrix RP, rightR, leftR;
  if (!bidirectional)
  { rightR = Matrix::Random(B.cols(), nProjections);
    RP = B * rightR;
  }
  else
  { assert (nProjections % 2 == 0);
    RP.resize(B.rows(),nProjections);
    rightR = Matrix::Random(B.cols(), nProjections/2);
    RP.rightCols(nProjections/2) = B * rightR;
    leftR = Matrix::Random(B.cols(), nProjections/2);
    RP.leftCols(nProjections/2) = B.transpose() * leftR;
  }
  std::clog << "MAIN: Completed random projection RP[" << RP.rows() << "x" << RP.cols() << "]" << endl;

  // convert data text lines into vectors
  int nLines (FileUtils::count_lines(regrFileName));
  std::clog << "MAIN: Building regressor matrix from " << nLines << " lines of input in file " << regrFileName << ".\n";
  std::ifstream is(regrFileName);
  Vocabulary::SparseMatrix S(nLines,vocabulary.n_types());
  Vector Y (nLines);
  vocabulary.fill_sparse_regr_design(Y, S, is);
  std::clog << "MAIN : sum of row 0 of S is " << S.row(0).sum() << endl;
  std::clog << "MAIN : sum of row 1 of S is " << S.row(1).sum() << endl;

  // compute dense projection coefficients, with room for Y and leading count column
  Matrix X (S.rows(), 2+RP.cols());
  X.col(0) = Y;
  X.col(1) = S * Vector::Ones(S.cols());  // total tokens in a line
  Vector irNorm (RP.colwise().norm().array().inverse());
  RP = RP * irNorm.asDiagonal();         // normalize projection vectors
  Vector isNorm (S.rows());
  for (int i=0; i<S.rows(); ++i)
    isNorm(i) = 1/S.row(i).norm();
  S = isNorm.asDiagonal() * S;
  X.rightCols(RP.cols()) = S * RP;
  std::clog << "MAIN: First 3 rows of X are\n" << X.topRows(3);

  // write to file
  std::ofstream os("/Users/bob/Desktop/regr.txt");
  os << " Y n";
  for(int i=0; i<RP.cols(); ++i) os << " X" << i;
  os << endl << X << endl;
  
  return 0;
}



void
parse_arguments(int argc, char** argv,
		string &fileName, int &oovThreshold, bool &bidirectional, int &nProjections, string &regrFileName)
{
  static struct option long_options[] = {
    {"vocab_file",    required_argument, 0, 'v'},
    {"regr_file",     required_argument, 0, 'i'},
    {"min_frequency", required_argument, 0, 'f'},
    {"bidirectional",       no_argument, 0, 'b'},
    {"n_projections", required_argument, 0, 'r'},
    {0, 0, 0, 0}                             // terminator 
  };
  int key;
  int option_index = 0;
  while (-1 !=(key = getopt_long (argc, argv, "v:i:f:br:", long_options, &option_index))) // colon means has argument
  {
    // std::cout << "Option key " << char(key) << " for option " << long_options[option_index].name << ", option_index=" << option_index << std::endl;
    switch (key)
    {
    case 'v' : { fileName       = optarg;                                  break; }
    case 'i' : { regrFileName   = optarg;                                  break; }
    case 'f' : { oovThreshold   = read_utils::lexical_cast<int>(optarg);   break; }
    case 'b' : { bidirectional  = true ;                                   break; }
    case 'r' : { nProjections   = read_utils::lexical_cast<int>(optarg);   break; }
    default  : { std::cout << "PARSE: Option not recognized; returning.\n";       }
    } // switch
  } 
}
