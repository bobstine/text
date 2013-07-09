#include "vocabulary.h"
#include "read_utils.h"
#include "file_utils.h"
#include "regex.h"

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
		string &vFileName, int &minFrequency, bool &bidirectional, int &nProjections, int &powerIterations, string &rFileName);

  
int main(int argc, char** argv)
{
  // read input options
  string vocabFileName   (  ""   );  // text used to build bigram, eigenwords
  string regrFileName    (  ""   );  // modeling text, with leading y followed by text
  int    powerIterations (   0   );
  int    minFrequency    (   3   );
  bool   bidirectional   ( false );
  int    nProjections    (  50   );
  int    bigramSkip      (   0   );
  parse_arguments(argc, argv, vocabFileName, minFrequency, bidirectional, nProjections, powerIterations, regrFileName);
  std::clog << "MAIN: regressor --vocab_file=" << vocabFileName << " --regr_file=" << regrFileName
	    << " --min_frequency=" << minFrequency << " --n_projections=" << nProjections << " --power_iter " << powerIterations;
  if (bidirectional) std::clog << " --bidirectional ";
  std::clog << endl;
  
  // build vocabulary from source
  Vocabulary vocabulary(vocabFileName, minFrequency);
  std::clog << "MAIN: " << vocabulary << endl;

  // compute bigram matrix
  Vocabulary::SparseMatrix B (vocabulary.n_types(), vocabulary.n_types());
  vocabulary.fill_sparse_bigram_matrix(B, bigramSkip);
  std::clog << "MAIN: Init sparse bigram B[" << B.rows() << "x" << B.cols() << "] from map; sum +/,B= " << B.sum() << endl;
  std::clog << "MAIN: Bigram row sums are (0)"  // these should check with those printed by the vocabulary (but for 1 pairs, not singles)
	    << B.row(0).sum() << "  (1)" << B.row(1).sum() << "  (2)"
	    << B.row(2).sum() << "  (3)" << B.row(3).sum() << "  (4)" << B.row(4).sum() << std::endl;

  // form random projections
  Matrix RP;
  {
    if (!bidirectional)
    { Matrix O = Matrix::Random(B.cols(), nProjections);
      RP = B * O;
      if (powerIterations)
      { Vocabulary::SparseMatrix BBt = B * B.transpose();
	while (powerIterations--)
	{ Matrix Q = Eigen::HouseholderQR<Matrix>(RP).householderQ() * Matrix::Identity(RP.rows(),RP.cols());  // Eigen trick for thin Q
	  RP = BBt * Q;
	}
      }
    }
    else
    { Matrix rightR, leftR;
      assert (nProjections % 2 == 0);
      int half = nProjections/2;
      RP.resize(B.rows(),nProjections);
      rightR = Matrix::Random(B.cols(), half);
      RP.rightCols(half) = B * rightR;
      leftR = Matrix::Random(B.cols(), half);
      RP.leftCols(half) = B.transpose() * leftR;
      if (powerIterations)
      { Vocabulary::SparseMatrix BBt = B * B.transpose();
	Vocabulary::SparseMatrix BtB = B.transpose() * B;
	while (powerIterations--)
	{ Matrix Q = Eigen::HouseholderQR<Matrix>(RP.rightCols(half)).householderQ()  * Matrix::Identity(RP.rows(),half);  // Eigen trick for thin Q
	  RP.rightCols(half) = BBt * Q;
	  Q = Eigen::HouseholderQR<Matrix>(RP.leftCols(half)).householderQ() * Matrix::Identity(RP.rows(), half);
	  RP.leftCols(half) = BtB * Q;
	 }
       }
     }
    std::clog << "MAIN: Completed random projection RP[" << RP.rows() << "x" << RP.cols() << "]";
    if (powerIterations) std::clog << " with power iterations.";
    std::clog << endl;
  }
  
  { // write eigenwords to file
    Vocabulary::TypeVector names = vocabulary.types();
    std::ofstream os("/Users/bob/Desktop/dictionary.txt");
    os << "Type";
    for (int i=0; i<RP.cols(); ++i) os << " RP" << i;
    os << std::endl;
    for (int i=0; i<RP.rows(); ++i)
      os << names[i] << " " << RP.row(i) << std::endl;
  }

  // convert data text lines into vectors
  int nLines (FileUtils::count_lines(regrFileName));
  std::clog << "MAIN: Building regressor matrix from " << nLines << " lines of input in file " << regrFileName << ".\n";
  std::ifstream is(regrFileName);
  Vocabulary::SparseMatrix S(nLines,vocabulary.n_types());
  Vector Y (nLines);
  vocabulary.fill_sparse_regr_design(Y, S, is);
  std::clog << "MAIN : sum of row 0 of S is " << S.row(0).sum() << endl;
  std::clog << "MAIN : sum of row 1 of S is " << S.row(1).sum() << endl;
  is.close();
  
  // parse for domain-specific attributes
  Vector sqft (nLines);  Vector sqftObserved = Vector::Zero(nLines);
  Vector bdrm (nLines);  Vector bdrmObserved = Vector::Zero(nLines);
  Vector bath (nLines);  Vector bathObserved = Vector::Zero(nLines);
  {
    double sqftSum (0.0); int nSqft (0);
    float  bdrmSum (0.0); int nBdrm (0);
    float  bathSum (0.0); int nBath (0);
    is.open(regrFileName);
    std::string line;
    for (int i=0; i<nLines; ++i)
    { std::getline(is, line);
      sqft(i) = square_footage(line);
      if (sqft(i)>0)
      { sqftObserved(i) = 1; sqftSum += sqft(i); ++nSqft; }      
      bdrm(i) = number_bedrooms(line);
      if (bdrm(i)>0)
      { bdrmObserved(i) = 1; bdrmSum += bdrm(i); ++nBdrm; }
      bath(i) = number_bathrooms(line);
      if (bath(i)>0)
      { bathObserved(i) = 1; bathSum += bath(i); ++nBath; }
    }
    float sqftMean = sqftSum/nSqft;
    float bdrmMean = bdrmSum/nBdrm;
    float bathMean = bathSum/nBath;
    for (int i=0; i<nLines; ++i)
    { if (0 == sqftObserved(i)) sqft(i) = sqftMean;
      if (0 == bdrmObserved(i)) bdrm(i) = bdrmMean;
      if (0 == bathObserved(i)) bath(i) = bathMean;
    }
  }

  // compute dense projection coefficients for common words
  Matrix X (S.rows(), 2+6+RP.cols());
  X.col(0) = Y;                                  // stuff Y into first column for output
  X.col(1) = S * Vector::Ones(S.cols());         // put total count n of type into second col
  X.col(2) = sqft;
  X.col(3) = sqftObserved;
  X.col(4) = bdrm;
  X.col(5) = bdrmObserved;
  X.col(6) = bath;
  X.col(7) = bathObserved;
  Vector irNorm (RP.colwise().norm().array().inverse());
  RP = RP * irNorm.asDiagonal();                 // normalize projection vector, sparse coefs so X has corr
  Vector isNorm (S.rows());
  for (int i=0; i<S.rows(); ++i)
    isNorm(i) = 1/S.row(i).norm();
  S = isNorm.asDiagonal() * S;
  X.rightCols(RP.cols()) = S * RP;
  std::clog << "MAIN: First 10 rows of X are\n" << X.topRows(10) << endl;
  
  // handle oov words
  { std::clog << "MAIN: 200 OOV words are...\n    ";
    int i=0;
    for(auto x : vocabulary.oov_map())
    { ++i;
      std::clog << " (" << x.first << "," << x.second << ")";
      if (i > 200) break;
    }
    std::clog << endl;
  }
  // write to file
  std::ofstream os("/Users/bob/Desktop/regr.txt");
  os << " Y n SqFt SqFt_Obs Bedrooms Bedroom_Obs Bathrooms Bathroom_Obs";
  for(int i=0; i<RP.cols(); ++i) os << " X" << i;
  os << endl << X << endl;
  
  return 0;
}



void
parse_arguments(int argc, char** argv,
		string &fileName, int &oovThreshold, bool &bidirectional, int &nProjections, int &powerIterations, string &regrFileName)
{
  static struct option long_options[] = {
    {"vocab_file",    required_argument, 0, 'v'},
    {"regr_file",     required_argument, 0, 'i'},
    {"min_frequency", required_argument, 0, 'f'},
    {"bidirectional",       no_argument, 0, 'b'},
    {"power_iter",    required_argument, 0, 'p'},
    {"n_projections", required_argument, 0, 'r'},
    {0, 0, 0, 0}                             // terminator 
  };
  int key;
  int option_index = 0;
  while (-1 !=(key = getopt_long (argc, argv, "v:i:f:bp:r:", long_options, &option_index))) // colon means has argument
  {
    // std::cout << "Option key " << char(key) << " for option " << long_options[option_index].name << ", option_index=" << option_index << std::endl;
    switch (key)
    {
    case 'v' : { fileName       = optarg;                                  break; }
    case 'i' : { regrFileName   = optarg;                                  break; }
    case 'f' : { oovThreshold   = read_utils::lexical_cast<int>(optarg);   break; }
    case 'b' : { bidirectional  = true ;                                   break; }
    case 'p' : { powerIterations= read_utils::lexical_cast<int>(optarg);   break; }
    case 'r' : { nProjections   = read_utils::lexical_cast<int>(optarg);   break; }
    default  : { std::cout << "PARSE: Option not recognized; returning.\n";       }
    } // switch
  } 
}
