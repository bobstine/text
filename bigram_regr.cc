/*
  This version of the regressor code does the bigram regression analysis
*/

#include "helpers.Template.h"
#include "vocabulary.h"
#include "read_utils.h"
#include "file_utils.h"
#include "timing.h"

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
		string &fileName,
		int &minFrequency, int &nProjections, char &adjust, bool &quadratic, int &powerIterations, int &seed,
		string &outputPath);
  
int main(int argc, char** argv)
{
  std::clog << "MAIN: Eigen will use " << Eigen::nbThreads() << " threads.\n";
  
  // read input options
  string fileName        (  ""   );  // used to build regression variables from document text (leading y_i)
  int    nSkipInitTokens (   1   );  // regression response at start of line (to isolate y_i from vocab)
  bool   markEndOfLine   ( false );  // avoid end-of-document token
  char   adjust          (  'w'  );  // use raw counts  (ignored if quadratic space)
  bool   quadratic       ( false );  // use second order token variables
  string outputPath      (  ""   );  // text files for parsed_#, lsa_#, bigram_#
  int    powerIterations (   0   );  // in Tropp algo for SVD
  int    minFrequency    (   3   );  // lower freq are treated as OOV
  int    nProjections    (  50   );
  int    randomSeed      ( 77777 );  // used to replicate random projection
  
  parse_arguments(argc, argv, fileName, minFrequency, nProjections, adjust, quadratic, powerIterations, randomSeed, outputPath);
  {
    std::string qStr = (quadratic) ? " --quadratic" : "";
    std::clog << "MAIN: regressor --file=" << fileName << " --output_path=" << outputPath << qStr
	      << " --min_frequency=" << minFrequency << " --n_projections=" << nProjections << " ---adjustment=" << adjust
	      << " --power_iter " << powerIterations << " --random_seed=" << randomSeed << endl;
  }
  string powerTag = "_p" + std::to_string(powerIterations);
  string wTag;
  switch (adjust)
    {
    case 'w' : { wTag = "raw"  ; break; }
    case 'r' : { wTag = "row"  ; break; }  // was sqrt
    case 'c' : { wTag = "col"  ; break; }
    case 'n' : { wTag = "recip"; break; }
    case 'b' : { wTag = "cca"  ; break; }
    default:   { wTag = "error"; std::clog << "MAIN: Unrecognized adjustment " << adjust << " given.\n"; return 0;}
    }
  
  // global random seed set here (controls random projections)
  srand(randomSeed);

  // build vocabulary
  Vocabulary vocabulary(fileName, nSkipInitTokens, markEndOfLine, minFrequency);
  std::clog << "MAIN: " << vocabulary << endl;
  {
    std::ofstream os (outputPath + "type_freq.txt");                        // write frequencies to file
    vocabulary.write_type_freq(os);
    int const maxNumberToWrite (200);                                       // write oov to screen
    vocabulary.write_oov_to_stream(std::clog, maxNumberToWrite);
  }

  // compute the document x type matrix W; also reads house prices (response) at head of input line
  int nDocs (FileUtils::count_lines(fileName));
  std::clog << "MAIN: Building word count matrix from " << nDocs << " lines of input in file " << fileName << ".\n";
  Vocabulary::SparseMatrix W(nDocs,vocabulary.n_types());
  Vector Y(nDocs), nTokens(nDocs);
  {
    bool sumToOne = false;            // avg of types (row sums to one; used to find centroids)
    if (adjust == 'n') sumToOne = true;
    std::ifstream is(fileName);
    vocabulary.fill_sparse_regr_design_from_stream(Y, W, nTokens, is, sumToOne);
    std::clog << "MAIN: Number tokens in first docs are "        << nTokens.head(5).transpose() << endl;
    std::clog << "MAIN: Leading block of the raw doc/word matrix W[" << W.rows() << "," << W.cols() << "] \n" << W.block(0,0,5,10) << endl;
  }

  // compute bigram matrix from vocabulary
  Vocabulary::SparseMatrix B (vocabulary.n_types(), vocabulary.n_types());
  const int bigramSkip = 0;
  vocabulary.fill_sparse_bigram_matrix(B, bigramSkip);
  std::clog << "MAIN: Bigram row sums are (0)"                              // match those from vocabulary if not scaled
	    << B.row(0).sum() << "  (1)" << B.row(1).sum() << "  (2)"
	    << B.row(2).sum() << "  (3)" << B.row(3).sum() << "  (4)" << B.row(4).sum() << std::endl;
  
  // adjustments to the bigram matrix
  if (adjust == 'r')
  { Vector sr = vocabulary.type_frequency_vector().array().sqrt().inverse();
    B = sr.asDiagonal() * B;
    std::clog << "MAIN: Leading block of the bigram matrix after ROW sqrt adjustment: \n" << B.block(0,0,5,10) << endl;
  }
  else if (adjust == 'c')
  { Vector st = vocabulary.type_frequency_vector().array().sqrt().inverse();
    B = B * st.asDiagonal();
    std::clog << "MAIN: Leading block of the bigram matrix after COL sqrt adjustment: \n" << B.block(0,0,5,10) << endl;
  }     
  else if (adjust == 'b')
  { Vector sc = vocabulary.type_frequency_vector().array().sqrt().inverse();
    B = sc.asDiagonal() * B * sc.asDiagonal();
    std::clog << "MAIN: Leading block of the bigram matrix after CCA adjustment: \n" << B.block(0,0,5,10) << endl;
  }

  // compute exact SVD decomposition
  if (false)
  { std::clog << "MAIN: Computing exact SVD of document-term matrix B begins.\n";
    Helper::write_exact_svd_to_path(B, nProjections, outputPath, wTag);
  }

  // P holds random projections version of SVD of bigram
  Matrix P(nDocs, nProjections);
  //  Helper::fill_random_projection(P, B, powerIterations);
  std::clog << "MAIN: Completed bigram projection.  P[" << P.rows() << "x" << P.cols() << "]\n";

  if (false)
  { std::clog << "MAIN: Writing eigenword matrix P to file eigenwords.txt." << std::endl;
    Helper::write_eigenwords_to_file (outputPath + "eigenwords.txt", P, vocabulary);
  }
  else std::clog << "MAIN: Skipping output of eigenword matrix to file.\n";

  // build centroid predictors A from eigenwords as average word position (see below for correlation)
  std::clog << "MAIN: Building centroid matrix A\n";
  Matrix A (nDocs, P.cols());
  for (int i=0; i<W.outerSize(); ++i)            // W is in row major order
  { Vector dest (Vector::Zero(P.cols()));
    int mi (0);
    for (Vocabulary::SparseMatrix::InnerIterator it(W,i); it; ++it)
    { mi += it.value();
      dest += P.row(it.col()) * it.value();
    }
    if(mi != nTokens(i)) std::clog << "MAIN: ********  Word count mismatched !!!\n";
    A.row(i) = dest.array() / mi;
  }

  // optionally compute R2 sequence of regression models
  if (false)
  { std::clog << "MAIN: Fitting regressions on bigram singular vectors.\n";
    Eigen::VectorXd YY(nDocs), mm(nDocs);                  // convert into double and take log for regression code
    YY = Y.cast<double>().array().log();
    mm = nTokens.cast<double>().array().log();
    mm = mm.array() - mm.sum()/mm.size();                  // center to reduce collinearity
    bool reverse (false);                                  // reverse tests low frequency words first
    std::string fileName (outputPath + "bigram_regr_fit_");
    if (nTokens.size() > 0) fileName += "with_m";
    else              fileName += "no_m";
    if (reverse) fileName += "_rev.txt";
    else         fileName += "_for.txt";
    const int degree = 5;
    Helper::calculate_sequence_r2 (YY, mm, degree, reverse, A, vocabulary, P.cols(), fileName);  // need to have strings, not words
  }


  // write to tab delimited output files if path assigned
  if (outputPath.size() > 0)
  { // prec, align, col sep, row sep, row pre, row suf, file pre, file suff
    Eigen::IOFormat fmt(Eigen::StreamPrecision,Eigen::DontAlignCols,"\t","\n","","","","");
    {
      Matrix YX (nDocs, 2);                             // y , m_i
      YX.col(0) = Y.array().log();                      // stuff log Y into first column for output
      YX.col(1) = nTokens;
      std::ofstream os (outputPath + "bigram_ym.txt");
      os << "Y\tm" << endl;
      os <<  YX.format(fmt) << endl;
    }
    {
      string label ="bigram_"  + wTag;
      char varSymbol = 'B';
      string dim  (std::to_string(nProjections));
      std::ofstream os (outputPath + label + "_" + dim + powerTag + ".txt");
      os << varSymbol << "0";
      for(int i=1; i<A.cols(); ++i) os << "\t" << varSymbol << i;
      os << endl << A.format(fmt) << endl;
    }
  }
  return 0;
}



void
parse_arguments(int argc, char** argv,
		string &fileName,
		int &oovThreshold, int &nProjections, char &adjust, bool &quadratic, int &powerIterations, int &seed, string &outputPath)
{
  static struct option long_options[] = {
    {"file",          required_argument, 0, 'i'},
    {"min_frequency", required_argument, 0, 'f'},
    {"adjustment",    required_argument, 0, 'a'},
    {"quadratic",     no_argument,       0, 'q'},
    {"power_iter",    required_argument, 0, 'p'},
    {"n_projections", required_argument, 0, 'r'},
    {"random_seed",   required_argument, 0, 's'},
    {"output_path",   required_argument, 0, 'o'},
    {0, 0, 0, 0}                             // terminator 
  };
  int key;
  int option_index = 0;
  while (-1 !=(key = getopt_long (argc, argv, "i:f:a:qp:r:s:o:", long_options, &option_index))) // colon means has argument
  {
    // std::cout << "Option key " << char(key) << " for option " << long_options[option_index].name << ", option_index=" << option_index << std::endl;
    switch (key)
    {
    case 'i' : { fileName       = optarg;                                  break; }
    case 'f' : { oovThreshold   = read_utils::lexical_cast<int>(optarg);   break; }
    case 'a' : { adjust         = *optarg;                                 break; }
    case 'q' : { quadratic      = true;                                    break; }
    case 'p' : { powerIterations= read_utils::lexical_cast<int>(optarg);   break; }
    case 'r' : { nProjections   = read_utils::lexical_cast<int>(optarg);   break; }
    case 's' : { seed           = read_utils::lexical_cast<int>(optarg);   break; }
    case 'o' : { outputPath     = optarg;                                  break; }
    default  : { std::cout << "PARSE: Option not recognized; returning.\n";       }
    } // switch
  } 
}
