#include "vocabulary.h"
#include "eigenword_dictionary.h"
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

void fill_random_bigram_projection(Matrix &P, Vocabulary::SparseMatrix const&B, int power)
{
  assert (B.rows() == P.rows());
  Matrix R = B * Matrix::Random(B.cols(), P.cols());                                         //  return orthog system
  P = Eigen::HouseholderQR<Matrix>(R).householderQ() * Matrix::Identity(P.rows(),P.cols());  //  Eigen trick for thin Q
  if (power > 0)
  { Vocabulary::SparseMatrix BBt = B * B.transpose();
    while (power--)
    { R = BBt * P;
      P = Eigen::HouseholderQR<Matrix>(R).householderQ() * Matrix::Identity(P.rows(),P.cols());
    }
  }
}

void
parse_arguments(int argc, char** argv,
		string &vFileName, string &rFileName,
		int &minFrequency, int &nProjections, int &powerIterations, string &outputFileName);

  
int main(int argc, char** argv)
{
  // should set random seed from input options...
  srand(2763);

  // read input options
  string vocabFileName   (  ""   );  // text used to build bigram, eigenwords
  string regrFileName    (  ""   );  // modeling text, with leading y followed by text
  string outputFileName  (  ""   );
  int    powerIterations (   0   );
  int    minFrequency    (   3   );
  int    nProjections    (  50   );
  parse_arguments(argc, argv, vocabFileName, regrFileName, minFrequency, nProjections, powerIterations, outputFileName);
  std::clog << "MAIN: lsa --vocab_file=" << vocabFileName << " --regr_file=" << regrFileName << " --output_file=" << outputFileName
	    << " --min_frequency=" << minFrequency << " --n_projections=" << nProjections << " --power_iter " << powerIterations;
  std::clog << endl;
  
  // build vocabulary
  Vocabulary vocabulary(vocabFileName, minFrequency);
  std::clog << "MAIN: " << vocabulary << endl;

  // compute context matrix from vocabulary and source lines (bags of word for each)
  Eigen::SparseMatrix<int,Eigen::RowMajor> C;
  {
    std::list< Eigen::Triplet<int> > triplets;
    std::ifstream is(vocabFileName);
    string line;
    int lineIndex = 0;
    while (getline(is, line))
    { std::istringstream ss (line);
      string token;
      std::map<int,int> lineMap;
      while (ss >> token)
	++lineMap[vocabulary.type_index(Type(token))];
      for(auto x : lineMap)
	triplets.push_back(Eigen::Triplet<int>(x.first, lineIndex, x.second));  // word index, context index, count
      ++lineIndex;
    }
    C.resize(vocabulary.n_types(), lineIndex);
    C.setFromTriplets(triplets.begin(), triplets.end());
  }
  std::clog << "MAIN: Dim of context matrix C is " << C.rows() << "x" << C.cols() << " with row sums "
	    << C.row(0).sum() << " " << C.row(1).sum() << " " << C.row(2).sum() << std::endl;
  
  // form random projections
  int half = nProjections/2;     assert(0 == nProjections%2);
  Matrix L(C.rows(), half);
  Matrix R(C.cols(), half);
  {
    fill_random_bigram_projection(L, C, powerIterations);
    Vocabulary::SparseMatrix Ct = C.transpose();
    fill_random_bigram_projection(R, Ct, powerIterations);
  }
  std::clog << "MAIN: Completed random projection.  L[" << L.rows() << "x" << L.cols() << "]    R[" << R.rows() << "x" << R.cols() << "]";
  if (powerIterations) std::clog << " with power iterations.";
  std::clog << endl;

  if (outputFileName.size()>0)
  { 
    if (false)  // this writes svd for the types
    { Vocabulary::TypeVector names = vocabulary.types();
      std::ofstream os("/Users/bob/Desktop/left.txt");
      os << "Type";
      for (int i=0; i<L.cols(); ++i) os << " L" << i;
      os << std::endl;
      for (int i=0; i<L.rows(); ++i)
	os << names[i] << " " << L.row(i) << std::endl;
    }
    std::ofstream os (outputFileName);
    if (os)
    { std::clog << "MAIN: Writing output file to " << outputFileName << std::endl;
      os << "Context";
      for (int i=0; i<R.cols(); ++i) os << " R" << i;
      os << std::endl;
      for (int i=0; i<R.rows(); ++i)
	os << i << " " << R.row(i) << std::endl;
    }
  }
  return 0;
}



void
parse_arguments(int argc, char** argv,
		string &fileName, string &regrFileName,
		int &oovThreshold,  int &nProjections, int &powerIterations, string &outputFileName)
{
  static struct option long_options[] = {
    {"vocab_file",    required_argument, 0, 'v'},
    {"regr_file",     required_argument, 0, 'i'},
    {"min_frequency", required_argument, 0, 'f'},
    {"power_iter",    required_argument, 0, 'p'},
    {"n_projections", required_argument, 0, 'r'},
    {"output_file",   required_argument, 0, 'o'},
    {0, 0, 0, 0}                             // terminator 
  };
  int key;
  int option_index = 0;
  while (-1 !=(key = getopt_long (argc, argv, "v:i:f:p:r:o:", long_options, &option_index))) // colon means has argument
  {
    // std::cout << "Option key " << char(key) << " for option " << long_options[option_index].name << ", option_index=" << option_index << std::endl;
    switch (key)
    {
    case 'v' : { fileName       = optarg;                                  break; }
    case 'i' : { regrFileName   = optarg;                                  break; }
    case 'f' : { oovThreshold   = read_utils::lexical_cast<int>(optarg);   break; }
    case 'p' : { powerIterations= read_utils::lexical_cast<int>(optarg);   break; }
    case 'r' : { nProjections   = read_utils::lexical_cast<int>(optarg);   break; }
    case 'o' : { outputFileName = optarg;                                  break; }
    default  : { std::cout << "PARSE: Option not recognized; returning.\n";       }
    } // switch
  } 
}

  
