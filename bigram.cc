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
		string &textFileName, int &minFrequency, int &nProjections, int &powerIterations, int &seed,
		string &outputPath);
  
int main(int argc, char** argv)
{
  // read input options
  string textFileName   (  ""   );   // text used to build bigram, eigenwords 
  int    nSkipInitTokens (   1   );  // regression response at start of line
  bool   markEndOfLine   ( true  );  // adds end-of-line tokens to separate text on distinct lines (as if paragraphs)
  string regrFileName    (  ""   );  // used to build regression variables (may be same as vocab file)
  string outputPath      (  ""   );  // text files for parsed_#, lsa_#, bigram_#
  int    powerIterations (   0   );  // in Tropp algo for SVD
  int    minFrequency    (   3   );  // lower freq are treated as OOV
  int    nProjections    (  50   );
  int    bigramSkip      (   0   );  // 0 is standard bigram
  int    randomSeed      ( 77777 );  // used to replicate random projection
  
  parse_arguments(argc, argv, textFileName, minFrequency, nProjections, powerIterations, randomSeed, outputPath);
  std::clog << "MAIN: bigram --text_file=" << textFileName << " --regr_file=" << regrFileName << " --output_file=" << outputPath
	    << " --min_frequency=" << minFrequency << " --n_projections=" << nProjections << " --power_iter " << powerIterations
	    << " --random_seed=" << randomSeed;
  srand(randomSeed);
  
  // build vocabulary
  Vocabulary vocabulary(textFileName, nSkipInitTokens, markEndOfLine, minFrequency);
  std::clog << "MAIN: " << vocabulary << endl;
  {
    std::ofstream os ("text_src/temp/type_freq.txt");                       // write frequencies to file
    vocabulary.write_type_freq(os);
  }
  int const maxNumberToWrite (200);
  vocabulary.write_oov_to_stream(std::clog, maxNumberToWrite);

  if (false)
    Helper::scan_google_vocabulary_for_oov (vocabulary);                    // check to see if oov words are in google
  
  // compute bigram matrix from vocabulary
  Vocabulary::SparseMatrix B (vocabulary.n_types(), vocabulary.n_types());
  vocabulary.fill_sparse_bigram_matrix(B, bigramSkip);
  std::clog << "MAIN: Bigram row sums are (0)"                              // match those from vocabulary if not scaled
	    << B.row(0).sum() << "  (1)" << B.row(1).sum() << "  (2)"
	    << B.row(2).sum() << "  (3)" << B.row(3).sum() << "  (4)" << B.row(4).sum() << std::endl;
  
  if (false) // exact decomposition via SVD
  { std::clog << "MAIN: Computing exact SVD of bigram matrix begins.\n";
    Helper::write_exact_svd_to_path(B, nProjections, outputPath);
  }
  
  Vocabulary::SparseMatrix Bt = B.transpose();
  Matrix Pl(B.rows(), nProjections);
  Matrix Pr(B.rows(), nProjections);
  Vector leftWeights  = Vector::Zero(0);
  Vector rightWeights = Vector::Zero(0);

  // first without weights
  std::clog << "MAIN: Writing raw components.\n";
  Helper::fill_random_projection(Pl, B, leftWeights, rightWeights, powerIterations);
  Helper::write_matrix_to_file(Pl, outputPath + "u_raw.txt", "L");
  Helper::fill_random_projection(Pr, Bt, rightWeights, leftWeights, powerIterations);
  Helper::write_matrix_to_file(Pr, outputPath + "v_raw.txt", "R");

  //  left weights
  std::clog << "MAIN: Writing left weighted components.\n";
  leftWeights = vocabulary.type_frequency_vector().array().inverse();
  rightWeights = Vector::Zero(0);
  Helper::fill_random_projection(Pl, B, leftWeights, rightWeights, powerIterations);
  Helper::write_matrix_to_file(Pl, outputPath + "u_left.txt", "L");
  Helper::fill_random_projection(Pr, Bt, rightWeights, leftWeights, powerIterations);
  Helper::write_matrix_to_file(Pr, outputPath + "v_left.txt", "R");

  //  right weights
  std::clog << "MAIN: Writing right weighted components.\n";
  leftWeights = Vector::Zero(0);
  rightWeights = vocabulary.type_frequency_vector().array().inverse();
  Helper::fill_random_projection(Pl, B, leftWeights, rightWeights, powerIterations);
  Helper::write_matrix_to_file(Pl, outputPath + "u_right.txt", "L");
  Helper::fill_random_projection(Pr, Bt, rightWeights, leftWeights, powerIterations);
  Helper::write_matrix_to_file(Pr, outputPath + "v_right.txt", "R");

  //  symmetric weights
  std::clog << "MAIN: Writing symmetrically weighted components.\n";
  rightWeights = vocabulary.type_frequency_vector().array().sqrt().inverse();
  Helper::fill_random_projection(Pl, B, rightWeights, rightWeights, powerIterations);
  Helper::write_matrix_to_file(Pl, outputPath + "u_sym.txt", "L");
  Helper::fill_random_projection(Pr, Bt, rightWeights, rightWeights, powerIterations);
  Helper::write_matrix_to_file(Pr, outputPath + "v_sym.txt", "R");

  return 0;
}



void
parse_arguments(int argc, char** argv, string &fileName, 
		int &oovThreshold, int &nProjections, int &powerIterations, int &seed, string &outputPath)
{
  static struct option long_options[] = {
    {"text_file",     required_argument, 0, 't'},
    {"min_frequency", required_argument, 0, 'f'},
    {"power_iter",    required_argument, 0, 'p'},
    {"n_projections", required_argument, 0, 'r'},
    {"random_seed",   required_argument, 0, 's'},
    {"output_path",   required_argument, 0, 'o'},
    {0, 0, 0, 0}                             // terminator 
  };
  int key;
  int option_index = 0;
  while (-1 !=(key = getopt_long (argc, argv, "t:f:bp:r:s:o:", long_options, &option_index))) // colon means has argument
  {
    // std::cout << "Option key " << char(key) << " for option " << long_options[option_index].name << ", option_index=" << option_index << std::endl;
    switch (key)
    {
    case 't' : { fileName       = optarg;                                  break; }
    case 'f' : { oovThreshold   = read_utils::lexical_cast<int>(optarg);   break; }
    case 'p' : { powerIterations= read_utils::lexical_cast<int>(optarg);   break; }
    case 'r' : { nProjections   = read_utils::lexical_cast<int>(optarg);   break; }
    case 's' : { seed           = read_utils::lexical_cast<int>(optarg);   break; }
    case 'o' : { outputPath     = optarg;                                  break; }
    default  : { std::cout << "PARSE: Option not recognized; returning.\n";       }
    } // switch
  } 
}

  
