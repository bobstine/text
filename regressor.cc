#include "vocabulary.h"
#include "eigenword_dictionary.h"
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
		string &outputFileName);

void
scan_google_vocabulary_for_oov (Vocabulary const& vocab);


void
write_eigenwords_to_file (string fileName, Matrix const& M, Vocabulary const& vocab)
{
  std::ofstream os(fileName);
  if (!os)
  { std::clog << "MAIN: Could not open file " << fileName << " for writing eigenwords.\n";
    return;
  }
  Vocabulary::TypeVector names = vocab.types();
  os << "Type";      // write column headers
  for (int i=0; i<M.cols(); ++i) os << " M" << i;
  os << std::endl;   // now data
  for (int i=0; i<M.rows(); ++i)
    os << names[i] << " " << M.row(i) << std::endl;
}

void
write_word_counts_to_file(std::string fileName, Vocabulary::SparseMatrix const& W, int nCols, Vocabulary const& vocab)
{
  std::ofstream os(fileName);
  if (!os)
  { std::clog << "MAIN: Could not open file " << fileName << " for writing eigenwords.\n";
    return;
  }
  Vocabulary::TypeVector tv = vocab.types();
  for (int i=0; i<nCols; ++i)   // write column headers
    os << " " << tv[i];
  os << std::endl;
  os << W.leftCols(nCols);
}


void
fill_random_projection(Matrix &P, Vocabulary::SparseMatrix const& M, Vector const& wts, int power)
{
  assert (M.rows() == P.rows());
  Matrix R;
  if(wts.size() > 0)
  { std::clog << "MAIN: Weighting in random projection.\n";
    R = wts.asDiagonal() * (M * Matrix::Random(M.cols(), P.cols()));
  }
  else
    R = M * Matrix::Random(M.cols(), P.cols());    
  P = Eigen::HouseholderQR<Matrix>(R).householderQ() * Matrix::Identity(P.rows(),P.cols());  //  Eigen trick for thin Q
  if (power > 0)
  { Vocabulary::SparseMatrix MMt = M * M.transpose();
    while (power--)
    { R = MMt * P;
      P = Eigen::HouseholderQR<Matrix>(R).householderQ() * Matrix::Identity(P.rows(),P.cols());
    }
  }
}

  
int main(int argc, char** argv)
{
  // read input options
  string vocabFileName   (  ""   );  // text used to build bigram, eigenwords with following options
  int    nSkipInitTokens (   1   );  // regression response at start of line
  bool   markEndOfLine   ( true  );  // adds end-of-line token
  string regrFileName    (  ""   );  // used to build regression variables (may be same as vocab file)
  string outputFileName  (  ""   );  // text file suitable for jmp or R
  int    powerIterations (   0   );  // in Tropp algo for SVD
  int    minFrequency    (   3   );  // lower freq are treated as OOV
  bool   bidirectional   ( false );  // left and right eigenvectors of B
  int    nProjections    (  50   );
  int    bigramSkip      (   0   );  // 0 is standard bigram
  int    randomSeed      ( 77777 );  // used to replicate random projection
  
  parse_arguments(argc, argv, vocabFileName, regrFileName, minFrequency, bidirectional, nProjections, powerIterations, randomSeed, outputFileName);
  std::clog << "MAIN: regressor --vocab_file=" << vocabFileName << " --regr_file=" << regrFileName << " --output_file=" << outputFileName
	    << " --min_frequency=" << minFrequency << " --n_projections=" << nProjections << " --power_iter " << powerIterations
	    << " --random_seed=" << randomSeed;
  srand(randomSeed);
  if (bidirectional) std::clog << " --bidirectional ";
  std::clog << endl;
  
  // build vocabulary
  Vocabulary vocabulary(vocabFileName, nSkipInitTokens, markEndOfLine, minFrequency);
  std::clog << "MAIN: " << vocabulary << endl;
  {
    std::ofstream os ("text_src/temp/type_freq.txt");                  // write frequencies to file
    vocabulary.write_type_freq(os);
  }  
  {
    std::clog << "MAIN: 200 OOV words are...\n    ";                     // display OOV tokens
    int i=0;
    for(auto x : vocabulary.oov_map())
    { ++i;
      std::clog << " (" << x.first << "," << x.second << ")";
      if (i > 200) break;
    }
    std::clog << endl;
  }
  if (false)
    scan_google_vocabulary_for_oov (vocabulary);                      // check to see if oov words are in google

  
  // compute bigram matrix from vocabulary
  Vocabulary::SparseMatrix B (vocabulary.n_types(), vocabulary.n_types());
  vocabulary.fill_sparse_bigram_matrix(B, bigramSkip);
  std::clog << "MAIN: Init sparse bigram B[" << B.rows() << "x" << B.cols() << "] from map; sum +/,B= " << B.sum() << endl;
  std::clog << "MAIN: Bigram row sums are (0)"  // these should check with those printed by the vocabulary (but for 1 pairs, not singles)
	    << B.row(0).sum() << "  (1)" << B.row(1).sum() << "  (2)"
	    << B.row(2).sum() << "  (3)" << B.row(3).sum() << "  (4)" << B.row(4).sum() << std::endl;
  
  // form random projections of bigram matrix
  Matrix P(B.rows(), nProjections);
  Vector emptyVec = Vector::Zero(0);
  if (!bidirectional)
    fill_random_projection(P, B, emptyVec, powerIterations);
  else
  { assert(0 == nProjections%2);
    int half = nProjections/2;
    Matrix Pl(B.rows(), half);
    fill_random_projection(Pl, B, emptyVec, powerIterations);
    Matrix Pr(B.rows(), half);
    Vocabulary::SparseMatrix Bt = B.transpose();
    fill_random_projection(Pr, Bt, emptyVec, powerIterations);
    P.leftCols(half) = Pl;
    P.rightCols(half)= Pr;
  }
  std::clog << "MAIN: Completed random projection P[" << P.rows() << "x" << P.cols() << "]";
  if (powerIterations) std::clog << " with power iterations.";
  std::clog << endl;
  if (false)
    write_eigenwords_to_file ("/Users/bob/Desktop/eigenwords.txt", P, vocabulary);
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
  Vector m = W * Vector::Ones(W.cols());                   // total count
  std::clog << "MAIN: Sum of row 0 of W is " << m(0) << "  sum of row 1 of W is " << m(1) << endl;

  // optionally write W to file
  if (false)
  { write_word_counts_to_file ("text_src/temp/w.txt", W, MIN(W.cols(), 2000), vocabulary);
    std::clog << "MAIN: Wrote W matrix to file w.txt." << std::endl;
  }
  else std::clog << "MAIN: Skipping output of W matrix to file.\n";
	
  // form random projections of LSA variables
  int half = nProjections/2;     assert(0 == nProjections%2);
  Matrix L(nLines, half);
  fill_random_projection(L, W, /* emptyVec */ m.array().inverse(), powerIterations);
  std::clog << "MAIN: Completed LSA random projection.  L[" << L.rows() << "x" << L.cols() << "]";
  if (powerIterations) std::clog << " with power iterations.";
  std::clog << endl;
  
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
  int offset = 2 + parsers.size();                // y , m_i, (parsed)
  Matrix X (W.rows(), offset+P.cols()+L.cols());
  X.col(0) = Y;                                  // stuff Y into first column for output
  X.col(1) = m;                                  // put total count n of type into second col (rowwise not avail for sparse)
  X.block(0,2, nLines, parsers.size()) = parsed; // custom variables

  // build variables from bigram correlations
  Vector irNorm (P.cols());
  irNorm = P.colwise().norm().array().inverse();
  P = P * irNorm.asDiagonal();                  // normalize projection vector, sparse coefs so X has corr
  Vector isNorm (W.rows());
  for (int i=0; i<W.rows(); ++i)
    isNorm(i) = 1/W.row(i).norm();
  W = isNorm.asDiagonal() * W;
  X.block(0,offset,W.rows(),P.cols()) = W * P;
  
  // lsa variables
  X.rightCols(L.cols()) = L;
  std::clog << "MAIN: First 10 rows and columns of X are\n" << X.block(0,0,5,10) << endl;
  
  // write to output file if assigned
  if (outputFileName.size() > 0)
  { std::ofstream os(outputFileName);
    if (!os)
    { std::clog << "MAIN: Could not open output file " << outputFileName << std::endl;
      return 0;
    }
    std::clog << "MAIN: Writing data file to " << outputFileName << std::endl;
    os << " Y m ";
    for (auto f : parsers) os << f.name() << " ";
    for(int i=0; i<P.cols()/2; ++i) os << " BL" << i;
    for(int i=0; i<P.cols()/2; ++i) os << " BR" << i;
    for(int i=0; i<L.cols(); ++i) os <<    " D" << i;
    os << endl << X << endl;
  }
  return 0;
}



void
parse_arguments(int argc, char** argv,
		string &fileName, string &regrFileName,
		int &oovThreshold, bool &bidirectional, int &nProjections, int &powerIterations, int &seed, string &outputFileName)
{
  static struct option long_options[] = {
    {"vocab_file",    required_argument, 0, 'v'},
    {"regr_file",     required_argument, 0, 'i'},
    {"min_frequency", required_argument, 0, 'f'},
    {"bidirectional",       no_argument, 0, 'b'},
    {"power_iter",    required_argument, 0, 'p'},
    {"n_projections", required_argument, 0, 'r'},
    {"random_seed",   required_argument, 0, 's'},
    {"output_file",   required_argument, 0, 'o'},
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
    case 'o' : { outputFileName = optarg;                                  break; }
    default  : { std::cout << "PARSE: Option not recognized; returning.\n";       }
    } // switch
  } 
}

void
scan_google_vocabulary_for_oov (Vocabulary const& vocab)
{
  // inspect google eigenwords for rare terms... did not find 940 (for chicago)
  if (false)
  { EigenwordDictionary googleDict ("text_src/eigenwords/google.txt", 100000, 43);
    int count = 0;
    std::clog << "MAIN: Google dictionary " << googleDict << " does not have types:\n   ";
    for (auto t : vocab.types())
      if (googleDict.type_index(t) < 0) { std::clog << " " << t; ++count; }
    std::clog << endl << "MAIN: Did not find " << count << " types in google bigram e-words.\n" << endl << endl;
  }
  // compare vocab to google... did not find 349 types (mostly parsing issues like - $ # @ (as in e-mail))
  if (false)
  { std::ifstream is ("text_src/google/vocab");
    std::map<Type,int> googleVocab;
    int counter = 0;
    while(true)
    { std::string line, token, numStr;
      int         count;
      if (!getline(is, line)) break;
      if (0 == line.size()) break;
      auto it = line.begin();
      while(it != line.end())
      { if(*it == '\t') break;
	token.push_back(*it++);
      }
      while(it != line.end())
	numStr.push_back(*it++);
      count = std::atoi(numStr.c_str());
      googleVocab[Type(token)] = count;
      if (counter++ < 60)
	std::clog << "    Google " << token << " has count " << count << endl;
    }
    std::clog << "MAIN: Read " << googleVocab.size() << " types from google vocabulary.\n";
    std::vector<Type> notFound;
    for (auto t : vocab.types())
    { if(googleVocab.find(t) == googleVocab.end())
	notFound.push_back(t);
    }
    std::clog << "MAIN: Did not find following " << notFound.size() << " types in google vocab:\n    ";
    for (auto x : notFound)
      std::clog << " " << x;
    std::clog << endl << endl;
  }
}

  // exact decomposition via SVD
  /*
  std::clog << "MAIN: Computing SVD of bigram matrix begins.\n";
  Eigen::JacobiSVD<Matrix> svd(B, Eigen::ComputeThinU|Eigen::ComputeThinV);
  Matrix U = svd.matrixU() * Matrix::Identity(B.rows(), nProjections/2);
  Matrix V = svd.matrixV() * Matrix::Identity(B.rows(), nProjections/2);
  Vector s = svd.singularValues();
  std::clog << "MAIN: Leading singular values are " << s.transpose().head(20) << endl;
  */
  
