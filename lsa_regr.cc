/*
  This version of the regressor code does the LSA regression analysis,
  optionally with quadratic expansion of the type space.
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

typedef Helper::ScalarType ScalarType;
typedef Eigen::VectorXf Vector;
typedef Eigen::MatrixXf Matrix;


void
parse_arguments(int argc, char** argv,
		string &fileName,  bool& hasY, bool &useLog, 
		int &minFrequency, int &nProjections, char &adjust, bool &quadratic, int &powerIterations, int &seed,
		string &outputPath);
  
int main(int argc, char** argv)
{
  std::clog << "MAIN: Eigen will use " << Eigen::nbThreads() << " threads.\n";
  
  // read input options
  string fileName        (  ""   );  // used to build regression variables from document text (leading y_i)
  bool   markEndOfLine   ( false );  // avoid end-of-document token
  char   adjust          (  'w'  );  // use raw counts  (ignored if quadratic space)
  bool   quadratic       ( false );  // use second order token variables
  bool   hasY            ( false );  // data includes a response
  bool   useLog          ( false );  // log y
  string outputPath      (  ""   );  // text files for parsed_#, lsa_#, bigram_#
  int    powerIterations (   0   );  // in Tropp algo for SVD
  int    minFrequency    (   3   );  // lower freq are treated as OOV
  int    nProjections    (  50   );
  int    randomSeed      ( 77777 );  // used to replicate random projection

  int    nSkipInitToken  = (hasY) ? 1 : 0;  // regression response at start of line (to isolate y_i from vocab)
  
  parse_arguments(argc, argv, fileName, hasY, useLog, minFrequency, nProjections, adjust, quadratic, powerIterations, randomSeed, outputPath);
  {
    std::string qStr = (quadratic) ? " --quadratic" : "";
    std::string yStr = (hasY)      ? " --has_y" : "";
    std::string lStr = (useLog)    ? " --use_log" : "";
    std::clog << "MAIN: lsa_regr " << yStr << lStr << " --file=" << fileName << " --output_path=" << outputPath << qStr 
	      << " --min_frequency=" << minFrequency << " --n_projections=" << nProjections << " ---adjustment=" << adjust
	      << " --power_iter " << powerIterations << " --random_seed=" << randomSeed << endl;
  }

  string powerTag = "_p" + std::to_string(powerIterations);
  string wTag;
  switch (adjust)
    {
    case 'w' : { wTag = "raw"  ; break; }
    case 't' : { wTag = "tfidf"; break; }
    case 'r' : { wTag = "row"  ; break; }  // was sqrt
    case 'c' : { wTag = "col"  ; break; }
    case 'n' : { wTag = "recip"; break; }
    case 'b' : { wTag = "cca"  ; break; }
    default:   { wTag = "error"; std::clog << "MAIN: Unrecognized adjustment " << adjust << " given.\n"; }
    }
  
  // global random seed set here (controls random projections)
  srand(randomSeed);

  // build vocabulary
  Vocabulary vocabulary(fileName, nSkipInitToken, markEndOfLine, minFrequency);
  std::clog << "MAIN: " << vocabulary << endl;
  {
    std::ofstream os (outputPath + "type_freq.txt");                        // write frequencies to file
    vocabulary.write_type_freq(os);
    int const maxNumberToWrite (200);                                       // write oov to screen
    vocabulary.write_oov_to_stream(std::clog, maxNumberToWrite);
  }

  // compute the document x type matrix W; optionally reads response at head of input line
  int nDocs (FileUtils::count_lines(fileName));
  std::clog << "MAIN: Building word count matrix from " << nDocs << " lines of input in file " << fileName << ".\n";
  Vocabulary::SparseMatrix W(nDocs,vocabulary.n_types());
  Vector Y(nDocs), nTokens(nDocs);
  {
    bool sumToOne = false;            // avg of types (row sums to one; used to find centroids)
    if (adjust == 'n') sumToOne = true;
    std::ifstream is(fileName);
    if (hasY)
      vocabulary.fill_sparse_regr_design_from_stream(Y, W, nTokens, is, sumToOne);
    else
    { Y.setZero();
      vocabulary.fill_sparse_regr_design_from_stream(W, nTokens, is, sumToOne);
    }
    std::clog << "MAIN: Number tokens in first docs are "        << nTokens.head(5).transpose() << endl;
    std::clog << "MAIN: Leading block of the raw LSA matrix W[" << W.rows() << "," << W.cols() << "] \n" << W.block(0,0,5,10) << endl;
  }

  // optionally write W to file
  if (false)
    { const int numWords (MIN(W.cols(), 6000));
      std::string name ("w" + std::to_string(numWords) + ".txt");
      Helper::write_word_counts_to_file (outputPath + name, W, numWords, vocabulary);
      std::clog << "MAIN: Wrote W matrix to file " << name << std::endl;
  }
  else std::clog << "MAIN: Skipping output of W matrix to file.\n";

  // adjustments to the elements of the term-document matrix
  Helper::scale_doc_term_matrix(adjust, &W, vocabulary, nTokens);

  if (false)   // compute negative of tf-idf, sort columns of W on this variable (so biggest tf-idf are first)
  { std::clog << "MAIN: Sorting top 5000 columns of W using tf-idf." << endl;
    Helper::sort_columns_using_tfidf (&W, vocabulary);
  }
  
  if (false) // compute exact SVD decomposition
  { std::clog << "MAIN: Computing exact SVD of document-term matrix W begins.\n";
    Helper::write_exact_svd_to_path(W, nProjections, outputPath, wTag);
  }

  // P holds random projection SVD of LSA variables
  std::clog << "MAIN: Preparing for random projection of W[" << W.rows() << "," << W.cols() << "]\n";
  Matrix P(nDocs, nProjections);
  Vector sv(nProjections);
  Matrix V (nDocs,25);
  if (! quadratic)                                                                          
    Helper::fill_random_projection_svd(&P,&sv,&V, W,powerIterations);
  else  // quadratic
  { std::clog << "MAIN: Computing random projection of " << (W.cols()*(W.cols()+1))/2 << " quadratics (excludes linear).";
    if (powerIterations) std::clog << " with power iterations.\n" ; else std::clog << ".\n";
    print_with_time_stamp("Starting base projection", std::clog);
    Matrix R = Matrix::Zero(nDocs, nProjections);
    Eigen::SparseMatrix<ScalarType, Eigen::ColMajor> X = W; 
    for (int j=0; j<W.cols(); ++j)
      for (int k=j; k<W.cols(); ++k)         // X'X = sum x_i x_i'
      { Eigen::SparseVector<ScalarType>  cp  = X.col(j).cwiseProduct(X.col(k));
	Vector                          rand = Vector::Random(nProjections);
	for (int i=0; i<nProjections; ++i)   // same as P += cp * rand.transpose(), but faster this way
	  R.col(i) += cp * rand(i);
      }
    P = Eigen::HouseholderQR<Matrix>(R).householderQ() * Matrix::Identity(P.rows(),P.cols());
    print_with_time_stamp("Complete base projection", std::clog);
    if (powerIterations)
    { print_with_time_stamp("Starting power iterations", std::clog);
      while (powerIterations--)
      { R = X * X.transpose() * P;
	P = Eigen::HouseholderQR<Matrix>(R).householderQ() * Matrix::Identity(P.rows(),P.cols());
      }
      print_with_time_stamp("Complete power iterations of quadratic projection", std::clog);
    }
  }
  std::clog << "MAIN: Completed LSA projection.  P[" << P.rows() << "x" << P.cols() << "]\n";

  
  // optionally compute R2 sequence of regression models
  if (true && hasY)
  { bool useM (false);                                     // use log token count in fitted model
    bool reverse (false);                                  // reverse tests low frequency words first
    std::clog << "MAIN: Fitting regressions on singular vectors.\n";
    Eigen::VectorXd YY(nDocs), mm(nDocs);                  // convert into double and take log for regression code
    YY = Y.cast<double>();
    if(useLog)
    { std::clog << "MAIN: Transforming response to log scale.\n";
      YY = YY.array().log();
    }
    if(useM)
    { std::clog << "MAIN: Adding log of tokens to fitted regression models.\n";
      mm = nTokens.cast<double>().array().log();
      mm = mm.array() - mm.sum()/mm.size();                // center to reduce collinearity
    }
    else
      mm.resize(0);
    std::string fileName (outputPath + "lsa_regr_fit_");
    if (mm.size() > 0) fileName += "with_m";
    else              fileName += "no_m";
    if (reverse) fileName += "_rev.txt";
    else         fileName += "_for.txt";
    const int degree = 5;
    Helper::calculate_sequence_r2 (YY, mm, degree, reverse, P, vocabulary, P.cols(), fileName);  // need to have strings, not words
  }


  // write to tab delimited output files if path assigned
  if (outputPath.size() > 0)
  { // prec, align, col sep, row sep, row pre, row suf, file pre, file suff
    Eigen::IOFormat fmt(Eigen::StreamPrecision,Eigen::DontAlignCols,"\t","\n","","","","");
    {
      Matrix YX (W.rows(), 2);                         // y , m_i
      if(useLog)                                       // stuff log Y into first column for output
      { std::clog << "MAIN: Writing file with response on log scale.\n";
	YX.col(0) = Y.array().log();
      }
      else
      { std::clog << "MAIN: Writing response on original scale.\n";
	YX.col(0) = Y;
      }
      YX.col(1) = nTokens;
      std::ofstream os (outputPath + "lsa_ym.txt");
      os << "Y\tm" << endl;
      os <<  YX.format(fmt) << endl;
    }
    {
      string label;
      char varSymbol;
      if (quadratic)
      { label = "lsaq_" + wTag;  varSymbol = 'Q'; }
      else
      { label = "lsa_"  + wTag; varSymbol = 'L'; }
      string dim  (std::to_string(nProjections));
      string prefix = outputPath + label + "_" + dim + powerTag;
      {
	std::ofstream os (prefix + ".txt");      // write the singular vectors
	os << varSymbol << "0";                                                     // write col labels for R
	for(int i=1; i<P.cols(); ++i) os << "\t" << varSymbol << i;
	os << endl << P.format(fmt) << endl;
      }
      {
	std::ofstream os2 (prefix + "_d.txt");  // write singular values
	os2 << sv.transpose().format(fmt) << endl;
      }
      {
	std::ofstream os3 (prefix + "_v.txt");
	os3 << "V0";                            // write col labels for V
	for(int i=1; i<V.cols(); ++i) os3 << "\t" << "V" << i;
	os3 << endl << V.format(fmt) << endl;
      }
    }
  }
  return 0;
}



void
parse_arguments(int argc, char** argv,
		string &fileName, bool& hasY, bool& useLog, 
		int &oovThreshold, int &nProjections, char &adjust, bool &quadratic, int &powerIterations, int &seed, string &outputPath)
{
  static struct option long_options[] = {
    {"file",          required_argument, 0, 'i'},
    {"min_frequency", required_argument, 0, 'f'},
    {"adjustment",    required_argument, 0, 'a'},
    {"quadratic",     no_argument,       0, 'q'},
    {"has_y",         no_argument,       0, 'y'},
    {"log_transform", no_argument,       0, 'l'},
    {"power_iter",    required_argument, 0, 'p'},
    {"n_projections", required_argument, 0, 'r'},
    {"random_seed",   required_argument, 0, 's'},
    {"output_path",   required_argument, 0, 'o'},
    {0, 0, 0, 0}                             // terminator 
  };
  int key;
  int option_index = 0;
  while (-1 !=(key = getopt_long (argc, argv, "i:f:a:qlyp:r:s:o:", long_options, &option_index))) // colon means has argument
  {
    // std::cout << "Option key " << char(key) << " for option " << long_options[option_index].name << ", option_index=" << option_index << std::endl;
    switch (key)
    {
    case 'i' : { fileName       = optarg;                                  break; }
    case 'f' : { oovThreshold   = read_utils::lexical_cast<int>(optarg);   break; }
    case 'a' : { adjust         = *optarg;                                 break; }
    case 'q' : { quadratic      = true;                                    break; }
    case 'y' : { hasY           = true;                                    break; }
    case 'l' : { useLog         = true;                                    break; }
    case 'p' : { powerIterations= read_utils::lexical_cast<int>(optarg);   break; }
    case 'r' : { nProjections   = read_utils::lexical_cast<int>(optarg);   break; }
    case 's' : { seed           = read_utils::lexical_cast<int>(optarg);   break; }
    case 'o' : { outputPath     = optarg;                                  break; }
    default  : { std::cout << "PARSE: Option not recognized; returning.\n";       }
    } // switch
  } 
}
