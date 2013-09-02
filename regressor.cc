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
		int &minFrequency, bool &bidirectional, int &nProjections, int &powerIterations, int &seed,
		string &outputFileName);

  
int main(int argc, char** argv)
{
  // read input options
  string vocabFileName   (  ""   );  // text used to build bigram, eigenwords
  string regrFileName    (  ""   );  // modeling text, with leading y followed by text
  string outputFileName  (  ""   );  // text file suitable for jmp or R
  int    powerIterations (   0   );
  int    minFrequency    (   3   );
  bool   bidirectional   ( false );
  int    nProjections    (  50   );
  int    bigramSkip      (   0   );
  int    randomSeed      ( 77777 );
  parse_arguments(argc, argv, vocabFileName, regrFileName, minFrequency, bidirectional, nProjections, powerIterations, randomSeed, outputFileName);
  std::clog << "MAIN: regressor --vocab_file=" << vocabFileName << " --regr_file=" << regrFileName << " --output_file=" << outputFileName
	    << " --min_frequency=" << minFrequency << " --n_projections=" << nProjections << " --power_iter " << powerIterations
	    << " --random_seed=" << randomSeed;
  srand(randomSeed);
  if (bidirectional) std::clog << " --bidirectional ";
  std::clog << endl;
  
  // build vocabulary
  Vocabulary vocabulary(vocabFileName, minFrequency);
  std::clog << "MAIN: " << vocabulary << endl;
  if(true) // write frequencies to file
  { std::ofstream os ("/Users/bob/Desktop/type_freq.txt");
    vocabulary.write_type_freq(os);
  }

  // google eigenwords for rare terms... did not find 940 (for chicago)
  if (false)
  { EigenwordDictionary googleDict ("text_src/eigenwords/google.txt", 100000, 43);
    int count = 0;
    std::clog << "MAIN: Google dictionary " << googleDict << " does not have types:\n   ";
    for (auto t : vocabulary.types())
      if (googleDict.type_index(t) < 0) { std::clog << " " << t; ++count; }
    std::clog << endl << "MAIN: Did not find " << count << " types in google bigram e-words.\n" << endl << endl;
  }

  // compare real estate vocab to google... did not find 349 types (mostly parsing issues like - $ # @ (as in e-mail))
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
    for (auto t : vocabulary.types())
    { if(googleVocab.find(t) == googleVocab.end())
	notFound.push_back(t);
    }
    std::clog << "MAIN: Did not find following " << notFound.size() << " types in google vocab:\n    ";
    for (auto x : notFound)
      std::clog << " " << x;
    std::clog << endl << endl;
  }
      
  // compute bigram matrix from vocabulary
  Vocabulary::SparseMatrix B (vocabulary.n_types(), vocabulary.n_types());
  vocabulary.fill_sparse_bigram_matrix(B, bigramSkip);
  std::clog << "MAIN: Init sparse bigram B[" << B.rows() << "x" << B.cols() << "] from map; sum +/,B= " << B.sum() << endl;
  std::clog << "MAIN: Bigram row sums are (0)"  // these should check with those printed by the vocabulary (but for 1 pairs, not singles)
	    << B.row(0).sum() << "  (1)" << B.row(1).sum() << "  (2)"
	    << B.row(2).sum() << "  (3)" << B.row(3).sum() << "  (4)" << B.row(4).sum() << std::endl;

  // form random projections
  Matrix P(B.rows(), nProjections);
  if (!bidirectional)
    fill_random_bigram_projection(P, B, powerIterations);
  else
  { assert(0 == nProjections%2);
    int half = nProjections/2;
    Matrix Pl(B.rows(), half);
    fill_random_bigram_projection(Pl, B, powerIterations);
    Matrix Pr(B.rows(), half);
    Vocabulary::SparseMatrix Bt = B.transpose();
    fill_random_bigram_projection(Pr, Bt, powerIterations);
    P.leftCols(half) = Pl;
    P.rightCols(half)= Pr;
  }
  std::clog << "MAIN: Completed random projection P[" << P.rows() << "x" << P.cols() << "]";
  if (powerIterations) std::clog << " with power iterations.";
  std::clog << endl;

  if (false)
  { // write eigenwords to file
    Vocabulary::TypeVector names = vocabulary.types();
    std::ofstream os("/Users/bob/Desktop/eigenwords.txt");
    os << "Type";
    for (int i=0; i<P.cols(); ++i) os << " P" << i;
    os << std::endl;
    for (int i=0; i<P.rows(); ++i)
      os << names[i] << " " << P.row(i) << std::endl;
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

  // form random projections
  int half = nProjections/2;     assert(0 == nProjections%2);
  Matrix L(S.rows(), half);
  fill_random_bigram_projection(L, S, powerIterations);
  std::clog << "MAIN: Completed LSA random projection.  L[" << L.rows() << "x" << L.cols() << "]";
  if (powerIterations) std::clog << " with power iterations.";
  std::clog << endl;
  
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
      sqft(i) = square_footage(line);                       // probably need some filtering limits here
      if ((0<sqft(i)) && (sqft(i)<50000))                   // not too big!
      { sqftObserved(i) = 1; sqftSum += sqft(i); ++nSqft; }      
      bdrm(i) = number_bedrooms(line);
      if (bdrm(i)>0)
      { bdrmObserved(i) = 1; bdrmSum += bdrm(i); ++nBdrm; }
      bath(i) = number_bathrooms(line);
      if (bath(i)>0)
      { bathObserved(i) = 1; bathSum += bath(i); ++nBath; }
    }
    std::clog << "MAIN: Observed counts for " << nSqft << " sq feet, " << nBdrm << " bedroom, " << nBath << " bathrooms.\n";
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
  int offset = 2 + 6;
  Matrix X (S.rows(), offset+P.cols()+L.cols());
  X.col(0) = Y;                                  // stuff Y into first column for output
  // custom variables
  X.col(1) = S * Vector::Ones(S.cols());         // put total count n of type into second col
  X.col(2) = sqft;
  X.col(3) = sqftObserved;
  X.col(4) = bdrm;
  X.col(5) = bdrmObserved;
  X.col(6) = bath;
  X.col(7) = bathObserved;
  // bigram variables
  Vector irNorm (P.colwise().norm().array().inverse());
  P = P * irNorm.asDiagonal();                 // normalize projection vector, sparse coefs so X has corr
  Vector isNorm (S.rows());
  for (int i=0; i<S.rows(); ++i)
    isNorm(i) = 1/S.row(i).norm();
  S = isNorm.asDiagonal() * S;
  X.block(0,offset,S.rows(),P.cols()) = S * P;
  // lsa variables
  X.rightCols(L.cols()) = L;
  std::clog << "MAIN: First 5 rows of X are\n" << X.topRows(5) << endl;
  // handle oov words
  {
    std::clog << "MAIN: 200 OOV words are...\n    ";
    int i=0;
    for(auto x : vocabulary.oov_map())
    { ++i;
      std::clog << " (" << x.first << "," << x.second << ")";
      if (i > 200) break;
    }
    std::clog << endl;
  }
  
  // write to output file if assigned
  if (outputFileName.size() > 0)
  { std::ofstream os(outputFileName);
    if(os)
    { std::clog << "MAIN: Writing data file to " << outputFileName << std::endl;
      os << " Y n SqFt SqFt_Obs Bedrooms Bedroom_Obs Bathrooms Bathroom_Obs";
      for(int i=0; i<P.cols(); ++i) os << " H" << i;
      for(int i=0; i<L.cols(); ++i) os << " L" << i;
      os << endl << X << endl;
    }
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


  // exact decomposition via SVD
  /*
  std::clog << "MAIN: Computing SVD of bigram matrix begins.\n";
  Eigen::JacobiSVD<Matrix> svd(B, Eigen::ComputeThinU|Eigen::ComputeThinV);
  Matrix U = svd.matrixU() * Matrix::Identity(B.rows(), nProjections/2);
  Matrix V = svd.matrixV() * Matrix::Identity(B.rows(), nProjections/2);
  Vector s = svd.singularValues();
  std::clog << "MAIN: Leading singular values are " << s.transpose().head(20) << endl;
  */
  
