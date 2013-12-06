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
		int &minFrequency, bool &bidirectional, int &nProjections, int &powerIterations, int &seed,
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
  bool   bidirectional   ( false );  // left and right eigenvectors of B
  int    nProjections    (  50   );
  int    bigramSkip      (   0   );  // 0 is standard bigram
  int    randomSeed      ( 77777 );  // used to replicate random projection
  
  parse_arguments(argc, argv, vocabFileName, regrFileName, minFrequency, bidirectional, nProjections, powerIterations, randomSeed, outputPath);
  std::clog << "MAIN: regressor --vocab_file=" << vocabFileName << " --regr_file=" << regrFileName << " --output_file=" << outputPath
	    << " --min_frequency=" << minFrequency << " --n_projections=" << nProjections << " --power_iter " << powerIterations
	    << " --random_seed=" << randomSeed;
  srand(randomSeed);
  if (bidirectional) std::clog << " --bidirectional ";
  std::clog << endl;
  
  // build vocabulary
  Vocabulary vocabulary(vocabFileName, nSkipInitTokens, markEndOfLine, minFrequency);
  std::clog << "MAIN: " << vocabulary << endl;
  {
    std::ofstream os ("text_src/temp/type_freq.txt");                       // write frequencies to file
    vocabulary.write_type_freq(os);
  }  
  {
    std::clog << "MAIN: 200 OOV words are...\n    ";                        // display OOV tokens
    int i=0;
    for(auto x : vocabulary.oov_map())
    { ++i;
      std::clog << " (" << x.first << "," << x.second << ")";
      if (i > 200) break;
    }
    std::clog << endl;
  }
  if (false)
    Helper::scan_google_vocabulary_for_oov (vocabulary);                    // check to see if oov words are in google
  
  // compute bigram matrix from vocabulary
  bool const useCorrScaling = true;                                         // normalize by diagonal counts
  Vocabulary::SparseMatrix B (vocabulary.n_types(), vocabulary.n_types());
  vocabulary.fill_sparse_bigram_matrix(B, bigramSkip, useCorrScaling);
  std::clog << "MAIN: Bigram row sums are (0)"                              // match those from vocabulary if not scaled
	    << B.row(0).sum() << "  (1)" << B.row(1).sum() << "  (2)"
	    << B.row(2).sum() << "  (3)" << B.row(3).sum() << "  (4)" << B.row(4).sum() << std::endl;
  
  if (false) // exact decomposition via SVD
  {
    std::clog << "MAIN: Computing exact SVD of bigram matrix begins.\n";
    std::ofstream os1 ("/home/bob/C/text/text_src/temp/ChicagoOld3/UDV.d.txt");
    if(!os1)
      std::clog << "MAIN: Could not open files for reporting exact SVD of bigram; skipping.\n";
    else
    { Eigen::JacobiSVD<Matrix> svd(B, Eigen::ComputeThinU|Eigen::ComputeThinV);
      Matrix U = svd.matrixU() * Matrix::Identity(B.rows(), nProjections);
      Matrix V = svd.matrixV() * Matrix::Identity(B.rows(), nProjections);
      Vector s = svd.singularValues();
      std::clog << "MAIN: Leading singular values are " << s.transpose().head(20) << endl;
      Eigen::IOFormat fmt(Eigen::StreamPrecision,Eigen::DontAlignCols,"\t","\n","","","","");
      os1 << s.transpose() << std::endl;
      os1.close();
      std::ofstream os2 ("/home/bob/C/text/text_src/temp/ChicagoOld3/UDV.u.txt");
      os2 << U.format(fmt) << std::endl;
      os2.close();
      std::ofstream os3 ("/home/bob/C/text/text_src/temp/ChicagoOld3/UDV.v.txt");
      os3 << V.format(fmt) << std::endl;
      os3.close();
    }
  }

  // form random projections of bigram matrix
  Matrix P(B.rows(), 2*nProjections);
  Vector noWeights = Vector::Zero(0);
  if (!bidirectional)
    Helper::fill_random_projection(P, B, noWeights, powerIterations);
  else
  { Matrix Pl(B.rows(), nProjections);
    Helper::fill_random_projection(Pl, B, noWeights, powerIterations);
    Matrix Pr(B.rows(), nProjections);
    Vocabulary::SparseMatrix Bt = B.transpose();
    Helper::fill_random_projection(Pr, Bt, noWeights, powerIterations);
    P.leftCols (nProjections) = Pl;
    P.rightCols(nProjections) = Pr;
  }
  std::clog << "MAIN: Completed random projection of bigram P[" << P.rows() << "x" << P.cols() << "]";
  if (powerIterations) std::clog << " with power iterations.";
  std::clog << endl;
  if (true)
    Helper::write_eigenwords_to_file ("text_src/temp/eigenwords.txt", P, vocabulary);
  else std::clog << "MAIN: Skipping output of eigenword matrix to file.\n";
  
  // convert data text lines into vectors for regression
  int nLines (FileUtils::count_lines(regrFileName));
  std::clog << "MAIN: Building regressor matrix from " << nLines << " lines of input in file " << regrFileName << ".\n";
  Vocabulary::SparseMatrix W(nLines,vocabulary.n_types());
  Vector Y (nLines);
  {
    std::ifstream is(regrFileName);
    vocabulary.fill_sparse_regr_design(Y, W, is);
  }
  Vector m = W * Vector::Ones(W.cols());                   // token count for each document
  std::clog << "MAIN: Sum of row 0 of W is " << m(0) << "  sum of row 1 of W is " << m(1) << endl;

  // optionally track R2 sequence of regression models for words
  if (false)
  { const int nColsRegr = 3000;
    Eigen::VectorXd YY(nLines), mm(nLines);                // put into double and take log for regression code
    for(int i=0; i<nLines; ++i)
    { YY(i) = log (Y(i));
      mm(i) = m(i);
    }
    bool reverse (false);                                  // reverse tests low frequency words first
    std::string fileName ("text_src/temp/word_regr_fit_");
    if (m.size() > 0) fileName += "with_m";
    else              fileName += "no_m";
    if (reverse) fileName += "_rev.txt";
    else         fileName += "_for.txt";
    Helper::calculate_sequence_r2 (YY, mm, reverse, W, vocabulary, nColsRegr, fileName);
  }
  
  // optionally write W to file
  if (false)
  { Helper::write_word_counts_to_file ("text_src/temp/w3000.txt", W, MIN(W.cols(), 3000), vocabulary);
    std::clog << "MAIN: Wrote W matrix to file w.txt." << std::endl;
  }
  else std::clog << "MAIN: Skipping output of W matrix to file.\n";

  // build regressor matrix A from eigenwords as average word position (see below for correlation)
  Matrix A (nLines, P.cols());
  for (int i=0; i<W.outerSize(); ++i)            // W is in row major order
  { Vector dest (Vector::Zero(P.cols()));
    int mi (0);
    for (Vocabulary::SparseMatrix::InnerIterator it(W,i); it; ++it)
    { mi += it.value();
      dest += P.row(it.col()) * it.value();
    }
    if(mi != m(i)) std::clog << "MAIN: ********  Word count mismatched !!!\n";
    A.row(i) = dest.array() / mi;
  }

  // form random projections of LSA variables
  Matrix L(nLines, nProjections);
  if (false) // exact, no weighting
  { std::clog << "MAIN: Computing left singular vectors of L by exact SVD.\n";
    Helper::fill_left_SVD(P, W);
  }
  else
  { std::clog << "MAIN: Computing left singular vectors of L by random projection";
    if (powerIterations) std::clog << " with power iterations.\n" ; else std::clog << "." << endl;
    Helper::fill_random_projection(L, W, noWeights /* m.array().inverse() */, powerIterations);
  }
  std::clog << "MAIN: Completed LSA projection.  L[" << L.rows() << "x" << L.cols() << "]\n";
  if (false)                                                // compute sequence of regressions for LSA variables
  { std::clog << "MAIN: Fitting regressions on singular vectors.\n";
    Eigen::VectorXd YY(nLines), mm(nLines);                // copy into double and take log for regression code
    for(int i=0; i<nLines; ++i)
    { YY(i) = log (Y(i));
      mm(i) = m(i);
    }
    std::string fileName ("text_src/temp/lsa_regr_fit_");
    if (m.size() > 0) fileName += "with_m";
    else              fileName += "no_m";
    Helper::calculate_sequence_r2 (YY, mm, "LSA_", L, fileName+".txt");
    fileName = "text_src/temp/bigram_regr_fit_";
    if (m.size() > 0) fileName += "with_m";
    else              fileName += "no_m";
    Helper::calculate_sequence_r2 (YY, mm, "Bigram_", A.leftCols(nProjections), fileName+".txt");
  }
  
  // parse domain-specific attributes
  const std::vector<Parser> parsers { Parser("SqFt", square_footage), Parser("Bedrooms", number_bedrooms), Parser("Bathrooms", number_bathrooms) };
  Matrix parsed = Matrix::Zero(nLines, parsers.size());
  {
    std::vector<int> counts (parsers.size(),0);
    std::ifstream is(regrFileName);
    std::string line;
    for (int i=0; i<nLines; ++i)
    { std::getline(is, line);
      int j=0;
      for (auto f : parsers)
      {	parsed(i,j) = f(line);
	if(parsed(i,j) != 0) ++counts[j];
	++j;
      }
    }
    std::clog << "MAIN: Observed listings generate (  ";
    for(size_t j=0; j<counts.size(); ++j)
      std::clog << counts[j] << "  ";
    std::clog << ") parsed counts." << std::endl;
  }

  // compute dense projection coefficients for common words
  Matrix YX (W.rows(), 2+parsers.size());         // y , m_i, (parsed)
  YX.col(0) = Y.array().log();                    // stuff log Y into first column for output
  YX.col(1) = m;                                  // put total count n of type into second col (rowwise not avail for sparse)
  YX.block(0,2, nLines, parsers.size()) = parsed; // custom variables
  std::clog << "MAIN: First 5 rows and columns of YX are\n"      << YX.block(0,0,5,5) << endl;
  std::clog << "MAIN: First 5 rows and columns of LSA are\n"     <<  L.block(0,0,5,5) << endl;
  std::clog << "MAIN: First 5 rows and columns of bigram are\n"  <<  A.block(0,0,5,5) << endl;

  // build variables from bigram correlations
  /*  Vector irNorm (P.cols());
      irNorm = P.colwise().norm().array().inverse();
      P = P * irNorm.asDiagonal();                  // normalize projection vector, sparse coefs so X has corr
      Vector isNorm (W.rows());
      for (int i=0; i<W.rows(); ++i)
      isNorm(i) = 1/W.row(i).norm();
      W = isNorm.asDiagonal() * W;
      X.block(0,offset,W.rows(),P.cols()) = W * P;
  */
  
  // write to tab delimited output files if path assigned
  if (outputPath.size() > 0)
  { // prec, align, col sep, row sep, row pre, row suf, file pre, file suff
    Eigen::IOFormat fmt(Eigen::StreamPrecision,Eigen::DontAlignCols,"\t","\n","","","","");
    string dim (std::to_string(nProjections));
    {
      string fileName (outputPath + "parsed.txt");
      std::ofstream os(fileName);
      if (!os)
      { std::clog << "MAIN: Could not open output file " << fileName << std::endl;
	return 0;
      }
      std::clog << "MAIN: Writing data file to " << fileName << std::endl;
      os << "Y\tm";
      for (auto f : parsers) os << "\t" << f.name() ;
      os << endl << YX.format(fmt) << endl;
    }
    {
      std::ofstream os (outputPath + "LSA_" + dim + ".txt");
      os << "D0";
      for(int i=1; i<L.cols(); ++i) os <<    "\tD" << i;
      os << endl << L.format(fmt) << endl;
    }
    {
      std::ofstream os (outputPath + "bigram_" + dim + ".txt");
      os << "BL0";
      for(int i=1; i<A.cols()/2; ++i) os << "\tBL" << i;
      for(int i=0; i<A.cols()/2; ++i) os << "\tBR" << i;
      os << endl << A.format(fmt) << endl;
    }
  }
  return 0;
}



void
parse_arguments(int argc, char** argv,
		string &fileName, string &regrFileName,
		int &oovThreshold, bool &bidirectional, int &nProjections, int &powerIterations, int &seed, string &outputPath)
{
  static struct option long_options[] = {
    {"vocab_file",    required_argument, 0, 'v'},
    {"regr_file",     required_argument, 0, 'i'},
    {"min_frequency", required_argument, 0, 'f'},
    {"bidirectional",       no_argument, 0, 'b'},
    {"power_iter",    required_argument, 0, 'p'},
    {"n_projections", required_argument, 0, 'r'},
    {"random_seed",   required_argument, 0, 's'},
    {"output_path",   required_argument, 0, 'o'},
    {0, 0, 0, 0}                             // terminator 
  };
  int key;
  int option_index = 0;
  while (-1 !=(key = getopt_long (argc, argv, "v:i:f:bp:r:s:o:", long_options, &option_index))) // colon means has argument
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
    case 's' : { seed           = read_utils::lexical_cast<int>(optarg);   break; }
    case 'o' : { outputPath     = optarg;                                  break; }
    default  : { std::cout << "PARSE: Option not recognized; returning.\n";       }
    } // switch
  } 
}

  // exact decomposition via SVD
  /*
  std::clog << "MAIN: Computing SVD of bigram matrix begins.\n";
  Eigen::JacobiSVD<Matrix> svd(B, Eigen::ComputeThinU|Eigen::ComputeThinV);
  Matrix U = svd.matrixU() * Matrix::Identity(B.rows(), nProjections);
  Matrix V = svd.matrixV() * Matrix::Identity(B.rows(), nProjections);
  Vector s = svd.singularValues();
  std::clog << "MAIN: Leading singular values are " << s.transpose().head(20) << endl;
  */
  
