/*
  This version of the regressor code 'unifies' the two approaches by computing a single SVD
  from the type x (doc type) matrix W'[D W_1].  That unifies the LSA analysis with the
  bigram to generate a single set of regressors rather than two (or more).

  Original regressor code includes many examples that explore other calculation options,
  such used to compare the exact SVD to the random projection.
*/

#include "helpers.h"
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


#define MIN(A,B) (((A)<(B)) ? (A) : (B))

using std::string;
using std::endl;

typedef Eigen::VectorXf Vector;
typedef Eigen::MatrixXf Matrix;


void
parse_arguments(int argc, char** argv,
		string &vocabFileName, string &regrFileName,
		int &minFrequency, int &nSkip, int &nProjections, int &powerIterations, int &seed,
		string &outputPath);
  
int main(int argc, char** argv)
{
  // read input options
  string vocabFileName   (  ""   );  // text used to build bigram, eigenwords with following options
  int    nSkipInitTokens (   1   );  // regression response at start of line
  bool   markEndOfLine   ( true  );  // adds end-of-line token
  string regrFileName    (  ""   );  // used to build regression variables (may be same as vocab file)
  string outputPath      (  ""   );  // text files for parsed_#, lsa_#, bigram_#
  int    powerIterations (   0   );  // in Tropp algo for SVD
  int    minFrequency    (   3   );  // lower freq are treated as OOV
  int    nProjections    (  50   );
  int    bigramSkip      (   0   );  // 0 is standard bigram
  int    randomSeed      ( 77777 );  // used to replicate random projection
  
  parse_arguments(argc, argv,
		  vocabFileName, regrFileName,
		  minFrequency, bigramSkip, nProjections, powerIterations, randomSeed, outputPath);
  std::clog << "MAIN: regressor --vocab_file=" << vocabFileName << " --regr_file=" << regrFileName << " --output_path=" << outputPath
	    << " --min_frequency=" << minFrequency << " --bigram_skip=" << bigramSkip << " --n_projections=" << nProjections
	    << " --power_iter " << powerIterations << " --random_seed=" << randomSeed << endl;
  
  // global random seed set here (controls random projections)
  srand(randomSeed);

  // build vocabulary
  Vocabulary vocabulary(vocabFileName, nSkipInitTokens, markEndOfLine, minFrequency);
  std::clog << "MAIN: " << vocabulary << endl;
  {
    std::ofstream os (outputPath + "type_freq.txt");                        // write frequencies to file
    vocabulary.write_type_freq(os);
    int const maxNumberToWrite (200);                                       // write oov to screen
    vocabulary.write_oov_to_stream(std::clog, maxNumberToWrite);
  }  

  if (false)
    Helper::scan_google_vocabulary_for_oov (vocabulary);                    // check to see if oov words are in google

  // compute the type x document matrix, with averages used; also reads off house prices (response)
  int nDocs (FileUtils::count_lines(regrFileName));
  std::clog << "MAIN: Building word count matrix from " << nDocs << " lines of input in file " << regrFileName << ".\n";
  Vocabulary::SparseMatrix W(nDocs,vocabulary.n_types());
  Vector Y(nDocs), nTokens(nDocs);
  {
    bool sumToOne = false;
    std::ifstream is(regrFileName);
    vocabulary.fill_sparse_regr_design_from_stream(Y, W, nTokens, is, sumToOne);
  }
  std::clog << "MAIN: Number tokens in first docs are "        << nTokens.head(5).transpose() << endl;
  {
    Vector checkSum = W * Vector::Ones(W.cols());                           // total for each document
    std::clog << "MAIN: Sums of leading rows of LSA matrix are " << checkSum.head(5).transpose() << endl;
  }
  
  // compute bigram matrix from vocabulary, weighted to have column sum 1
  Vocabulary::SparseMatrix B (vocabulary.n_types(), vocabulary.n_types());
  vocabulary.fill_sparse_bigram_matrix(B, bigramSkip);
  Vector bCounts = B * Vector::Ones(B.cols());
  Vector scaleWeights = vocabulary.type_frequency_vector().array().inverse();
  B = B * scaleWeights.asDiagonal();

  // transpose the LSA matrix
  Vocabulary::SparseMatrix Wt = W.transpose();
  assert(Wt.rows() /* types */ == B.rows());
  Vector wCounts = Wt * Vector::Ones(Wt.cols());
  scaleWeights = nTokens.array().inverse();
  Wt = Wt * scaleWeights.asDiagonal();
  
  // join two sparse matrices, [B Wt]
  Vocabulary::SparseMatrix C (B.rows(), B.cols()+W.rows());
  C.reserve(bCounts + wCounts);
  for (int r=0; r<B.outerSize(); ++r)
    for (Vocabulary::SparseMatrix::InnerIterator it(B,r); it; ++it)
      C.coeffRef(it.row(), it.col()) = it.value();
  for (int r=0; r<Wt.outerSize(); ++r)
    for (Vocabulary::SparseMatrix::InnerIterator it(Wt,r); it; ++it)
      C.coeffRef(it.row(), it.col()+B.cols()) = it.value();
  {
    Vector sumB     = B  * Vector::Ones(B.cols());
    Vector sumWt    = Wt * Vector::Ones(Wt.cols());
    Vector checkSum = C  * Vector::Ones(C.cols());
    std::clog << "DEBUG: sum(B)=" << sumB.sum() << "  sum(Wt)=" << sumWt.sum() << "  sum(C)=" << checkSum.sum() << endl;
    std::clog << "DEBUG: row sums of B    " << sumB.head(10).transpose() << endl;
    std::clog << "                   Wt   " << sumWt.head(10).transpose() << endl;
    std::clog << "                   C    " << checkSum.head(10).transpose() << endl;
  }
  
  // form random projections of unified matrix
  std::clog << "MAIN: Building random projection for unified matrix.\n";
  Matrix P(C.rows(), nProjections);
  Vector noWeights = Vector::Zero(0);
  Helper::fill_random_projection(P, C, noWeights, noWeights, powerIterations);
  std::clog << "MAIN: Completed random projection of bigram P[" << P.rows() << "x" << P.cols() << "]";
  if (powerIterations) std::clog << " with power iteration.";
  std::clog << endl;

  if (false)
  { std::clog << "MAIN: Writing eigenword matrix P to file eigenwords.txt." << std::endl;
    Helper::write_eigenwords_to_file (outputPath + "eigenwords.txt", P, vocabulary);
  }
  else std::clog << "MAIN: Skipping output of eigenword matrix to file.\n";
  
  // build centroid predictors eigenwords as average word position
  std::clog << "MAIN: Building centroid matrix A [" << nDocs << "," << P.cols() << "]\n";
  Matrix A (nDocs, P.cols());
  for (int i=0; i<W.outerSize(); ++i)            // W is in row major order
  { Vector centroid = Vector::Zero(P.cols());
    for (Vocabulary::SparseMatrix::InnerIterator it(W,i); it; ++it)
      centroid += P.row(it.col()) * it.value();
    A.row(i) = centroid.array() / nTokens(i);
  }

  // compute dense projection coefficients for common words
  Matrix YX (W.rows(), 2);                         // y , m_i
  YX.col(0) = Y.array().log();                     // stuff log Y into first column for output
  YX.col(1) = nTokens;

  // write to tab delimited output files if path assigned
  if (outputPath.size() > 0)
  { // prec, align, col sep, row sep, row pre, row suf, file pre, file suff
    Eigen::IOFormat fmt(Eigen::StreamPrecision,Eigen::DontAlignCols,"\t","\n","","","","");
    {
      std::ofstream os (outputPath + "unified_y_m.txt");
      os << "Y\tm" << endl;
      os <<  YX.format(fmt) << endl;
    }
    {
      string dim (std::to_string(nProjections));
      std::ofstream os (outputPath + "unified_svd_" + dim + ".txt");
      os << "Un0";
      for(int i=1; i<A.cols(); ++i) os << "\tUn" << i;
      os << endl << A.format(fmt) << endl;
    }
  }

  return 0;
}



void
parse_arguments(int argc, char** argv,
		string &fileName, string &regrFileName,
		int &oovThreshold, int &bigramSkip, int &nProjections, int &powerIterations, int &seed, string &outputPath)
{
  static struct option long_options[] = {
    {"vocab_file",    required_argument, 0, 'v'},
    {"regr_file",     required_argument, 0, 'i'},
    {"min_frequency", required_argument, 0, 'f'},
    {"bigram_skip",   required_argument, 0, 'k'},
    {"power_iter",    required_argument, 0, 'p'},
    {"n_projections", required_argument, 0, 'r'},
    {"random_seed",   required_argument, 0, 's'},
    {"output_path",   required_argument, 0, 'o'},
    {0, 0, 0, 0}                             // terminator 
  };
  int key;
  int option_index = 0;
  while (-1 !=(key = getopt_long (argc, argv, "v:i:f:k:p:r:s:o:", long_options, &option_index))) // colon means has argument
  {
    // std::cout << "Option key " << char(key) << " for option " << long_options[option_index].name << ", option_index=" << option_index << std::endl;
    switch (key)
    {
    case 'v' : { fileName       = optarg;                                  break; }
    case 'i' : { regrFileName   = optarg;                                  break; }
    case 'f' : { oovThreshold   = read_utils::lexical_cast<int>(optarg);   break; }
    case 'k' : { bigramSkip     = read_utils::lexical_cast<int>(optarg);   break; }
    case 'p' : { powerIterations= read_utils::lexical_cast<int>(optarg);   break; }
    case 'r' : { nProjections   = read_utils::lexical_cast<int>(optarg);   break; }
    case 's' : { seed           = read_utils::lexical_cast<int>(optarg);   break; }
    case 'o' : { outputPath     = optarg;                                  break; }
    default  : { std::cout << "PARSE: Option not recognized; returning.\n";       }
    } // switch
  } 
}
