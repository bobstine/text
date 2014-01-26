/*
  This version of the regressor code does the LSA regression analysis,
  optionally with quadratic expansion of the type space.
*/

#include "helpers.h"
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


/*
  Adjustments to raw frequencies:
    r = raw counts
    t = tf-idf
    s = divide by sqrt n_i
    n = divide by n_i
    c = CCA adjustment
*/

void
parse_arguments(int argc, char** argv,
		string &fileName,
		int &minFrequency, int &nProjections, char &adjust, bool &quadratic, int &powerIterations, int &seed,
		string &outputPath);
  
int main(int argc, char** argv)
{
  // read input options
  string fileName        (  ""   );  // used to build regression variables from document text (leading y_i)
  int    nSkipInitTokens (   1   );  // regression response at start of line (to isolate y_i from vocab)
  bool   markEndOfLine   ( false );  // avoid end-of-document token
  char   adjust          (  'r'  );  // use raw counts  (ignored if quadratic space)
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
    case 'r' : { wTag = "raw"  ; break; }
    case 't' : { wTag = "tfidf"; break; }
    case 's' : { wTag = "sqrt" ; break; }  
    case 'n' : { wTag = "recip"; break; }
    case 'c' : { wTag = "cca"  ; break; }
    default:   { wTag = "error"; std::clog << "MAIN: Unrecognized adjustment " << adjust << " given.\n"; }
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
    std::clog << "MAIN: Leading block of the LSA matrix after tf-idf adjustment: \n" << W.block(0,0,5,10) << endl;
  }

  // adjustments to the elements of the term-document matrix
  if (adjust == 's')
  { Vector sr = (W * Vector::Ones(W.cols())).array().sqrt().inverse();
    W = sr.asDiagonal() * W;
    std::clog << "MAIN: Leading block of the LSA matrix after sqrt adjustment: \n" << W.block(0,0,5,10) << endl;
  }
  else if (adjust == 'c')
  { Vector sr = nTokens.array().sqrt().inverse();
    Vector st = vocabulary.type_frequency_vector().array().sqrt().inverse();
    W = sr.asDiagonal() * W * st.asDiagonal();
    std::clog << "MAIN: Leading block of the LSA matrix after CCA adjustment: \n" << W.block(0,0,5,10) << endl;
  }     
  else if (adjust == 't')
  { Vector termCts = Vector::Zero(W.cols());   // count docs in which terms appear
    for (int doc=0; doc<W.outerSize(); ++doc)
      for (Vocabulary::SparseMatrix::InnerIterator it(W,doc); it; ++it)
	++termCts(it.col());
    for (int doc=0; doc<W.outerSize(); ++doc)
      for (Vocabulary::SparseMatrix::InnerIterator it(W,doc); it; ++it)
	W.coeffRef(doc, it.col()) = it.value() * log(W.rows()/termCts(it.col()));
    std::clog << "MAIN: Leading block of the LSA matrix after tf-idf adjustment: \n" << W.block(0,0,5,10) << endl;
  }
  
  if (true) // compute exact SVD decomposition
  { std::clog << "MAIN: Computing exact SVD of bigram matrix begins.\n";
    Helper::write_exact_svd_to_path(W, nProjections, outputPath, wTag);
  }

  // P holds random projections of LSA variables
  Matrix P(nDocs, nProjections);
  if (! quadratic)                                                                             // adapted from Helper::fill_random_projection
  { std::clog << "MAIN: Computing left singular vectors of L by random projection";
    if (powerIterations) std::clog << " with power iterations.\n" ; else std::clog << ".\n";
    print_with_time_stamp("Starting base linear random projection", std::clog);
    Matrix R = W * Matrix::Random(W.cols(), P.cols());    
    P = Eigen::HouseholderQR<Matrix>(R).householderQ() * Matrix::Identity(P.rows(),P.cols());  // block does not work; use to get left P.cols()
    print_with_time_stamp("Completed base linear random projection", std::clog);
    if (powerIterations)
      while (powerIterations--)
      { R = W * W.transpose() * P;   // IS THIS RIGHT IF W'W != I???
	print_with_time_stamp("Preparing Householder step of iterated random projection", std::clog);
	P = Eigen::HouseholderQR<Matrix>(R).householderQ() * Matrix::Identity(P.rows(),P.cols());
      }
    print_with_time_stamp("Completed power iterations of linear projection", std::clog);
    std::clog << "MAIN: Checking norms of leading terms in random projection; 0'0="
	      << P.col(0).dot(P.col(0)) << "   0'1=" << P.col(0).dot(P.col(1)) << "   1'1=" << P.col(1).dot(P.col(1)) << std::endl;
  }
  else  // quadratic
  { std::clog << "MAIN: Computing random projection of " << (W.cols()*(W.cols()+1))/2 << " quadratics (excludes linear).";
    if (powerIterations) std::clog << " with power iterations.\n" ; else std::clog << ".\n";
    print_with_time_stamp("Starting base projection", std::clog);
    Matrix R = Matrix::Zero(nDocs, nProjections);
    Eigen::SparseMatrix<float, Eigen::ColMajor> X = W; 
    for (int j=0; j<W.cols(); ++j)
      for (int k=j; k<W.cols(); ++k)         // X'X = sum x_i x_i'
      { Eigen::SparseVector<float>  cp  = X.col(j).cwiseProduct(X.col(k));
	Vector                     rand = Vector::Random(nProjections);
	for (int i=0; i<nProjections; ++i)   // same as P += cp * rand.transpose(), but faster this way
	  R.col(i) += cp * rand(i);
      }
    P = Eigen::HouseholderQR<Matrix>(R).householderQ() * Matrix::Identity(P.rows(),P.cols());
    print_with_time_stamp("Complete base projection", std::clog);
    if (powerIterations)
    { print_with_time_stamp("Starting power iterations", std::clog);
      /* Much slower if form XXt and multiply unless have very very few rows in X
	 Eigen::SparseMatrix<float,Eigen::ColMajor> XXt(X.rows(),X.rows());
	 XXt.setZero();
	 for (int j=0; j<W.cols(); ++j)
	 {	for (int k=j; k<W.cols(); ++k)
	 { Eigen::SparseVector<float> cp = X.col(j).cwiseProduct(X.col(k));
	 for (Eigen::SparseVector<float>::InnerIterator it(cp); it; ++it)
	 XXt.col(it.index()) += cp * it.value();
	 }	 }  */
      while (powerIterations--)
      { R = X * X.transpose() * P;
	P = Eigen::HouseholderQR<Matrix>(R).householderQ() * Matrix::Identity(P.rows(),P.cols());
      }
      print_with_time_stamp("Complete power iterations of quadratic projection", std::clog);
    }
  }
  if (adjust == 'c')
  { Vector sr = nTokens.array().sqrt().inverse();
    P = sr.asDiagonal() * P;
  }

  std::clog << "MAIN: Completed LSA projection.  P[" << P.rows() << "x" << P.cols() << "]\n";
  if (false)                                                // compute sequence of regressions for LSA variables
  { std::clog << "MAIN: Fitting regressions on singular vectors.\n";
    Eigen::VectorXd YY(nDocs), mm(nDocs);                // copy into double and take log for regression code
    for(int i=0; i<nDocs; ++i)
    { YY(i) = log (Y(i));
      mm(i) = nTokens(i);
    }
    std::string fileName (outputPath + "lsa_regr_fit_");
    if (nTokens.size() > 0) fileName += "with_m";
    else              fileName += "no_m";
    Helper::calculate_sequence_r2 (YY, mm, "LSA_", P, fileName+".txt");
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
      std::ofstream os (outputPath + "lsa_ym.txt");
      os << "Y\tm" << endl;
      os <<  YX.format(fmt) << endl;
    }
    {
      string label;
      char varSymbol;
      if (quadratic)
      { label = "lsaq_" + wTag + powerTag;  varSymbol = 'Q'; }
      else
      { label = "lsa_"  + wTag + powerTag; varSymbol = 'L'; }
      string dim  (std::to_string(nProjections));
      std::ofstream os (outputPath + label + "_" + dim + ".txt");
      os << varSymbol << "0";
      for(int i=1; i<P.cols(); ++i) os << "\t" << varSymbol << i;
      os << endl << P.format(fmt) << endl;
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
