#include "helpers.h"

#include "vocabulary.h"
#include "eigenword_dictionary.h"
#include "read_utils.h"
#include "file_utils.h"
#include "regex.h"
#include "regression.h"
#include "timing.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SVD>


using std::string;
using std::endl;


float
Helper::entropy(Vector const& prob)                     // assumes sum p(i) = 1 and p(i) >= 0
{
  float entropy = 0.0;
  for(int i = 0; i<prob.size(); ++i)
    if(0 < prob(i))
      entropy += prob(i) * log(prob(i));
  if (entropy < 0)
    return -entropy/log(prob.size());
  else
    return 0.0;
}

float
Helper::entropy(Vector const& x, float xSum)           // assumes x(i) >= 0
{
  float entropy = 0.0;
  float logSum = log(xSum);
  for(int i = 0; i<x.size(); ++i)
    if(0 < x(i))
      entropy += x(i) * (log(x(i)) - logSum);
  if (entropy != 0)
  { entropy /= xSum;
    return -entropy/log(x.size());
  }
  else
    return 0.0;
}

  
void
Helper::scan_google_vocabulary_for_oov (Vocabulary const& vocab)
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
    { string line, token, numStr;
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


void
Helper::write_eigenwords_to_file (string fileName, Matrix const& M, Vocabulary const& vocab)
{
  std::ofstream os(fileName);
  if (!os)
  { std::clog << "MAIN: Could not open file " << fileName << " for writing eigenwords.\n";
    return;
  }
  Vocabulary::TypeVector names = vocab.types();
  os << "Type";      // write column headers
  for (int i=0; i<M.cols(); ++i) os << " M" << i;
  os << endl;   // now data
  for (int i=0; i<M.rows(); ++i)
    os << names[i] << " " << M.row(i) << endl;
}

void
Helper::write_word_counts_to_file(string fileName, Vocabulary::SparseMatrix const& W, int nCols, Vocabulary const& vocab)
{
  std::ofstream os(fileName);
  if (!os)
  { std::clog << "MAIN: Could not open file " << fileName << " for writing eigenwords.\n";
    return;
  }
  Vocabulary::TypeVector tv = vocab.types();
  for (int i=0; i<nCols; ++i)   // write column headers
    os << " " << tv[i];
  os << endl;
  os << W.leftCols(nCols);
}



void
Helper::write_matrix_to_file(Matrix const& m, std::string fileName, std::string prefix)
{
  std::ofstream os(fileName);
  if (!os)
  { std::clog << "MAIN: Could not open matrix output file " << fileName << std::endl;
    return;
  }
  os << prefix + "0";
  prefix = "\t" + prefix;
  for(int i=1; i<m.cols(); ++i)
    os << prefix << i;
  // prec, align, col sep, row sep, row pre, row suf, file pre, file suff
  Eigen::IOFormat fmt(Eigen::StreamPrecision,Eigen::DontAlignCols,"\t","\n","","","","");
  os << std::endl << m.format(fmt) << endl;
}


void
Helper::write_exact_svd_to_path(Vocabulary::SparseMatrix const& B, int nProjections, std::string path, std::string tag)
{
  std::ofstream os1 (path + "svd_exact_d_" + tag + ".txt");
  if(!os1)
  { std::cerr << "MAIN: Invalid path. Could not open files for reporting exact SVD of bigram; skipping calculation.\n";
    return;
  }
  Eigen::JacobiSVD<Matrix> svd(B, Eigen::ComputeThinU|Eigen::ComputeThinV);
  Matrix U = svd.matrixU() * Matrix::Identity(B.cols(), nProjections);
  Matrix V = svd.matrixV() * Matrix::Identity(B.cols(), nProjections);
  Vector s = svd.singularValues();
  std::clog << "MAIN: Leading singular values are " << s.transpose().head(20) << endl;
  Eigen::IOFormat fmt(Eigen::StreamPrecision,Eigen::DontAlignCols,"\t","\n","","","","");
  os1 << s.transpose() << std::endl;
  os1.close();
  std::ofstream os2 (path + "svd_exact_u_" + tag + ".txt");
  os2 << U.format(fmt) << std::endl;
  os2.close();
  std::ofstream os3 (path + "svd_exact_v_" + tag + ".txt");
  os3 << V.format(fmt) << std::endl;
  os3.close();
}



void
Helper::fill_random_projection(Matrix &P, Vocabulary::SparseMatrix const& B, int powerIterations)
{
  std::clog << "HLPR: Computing left singular vectors of matrix by random projection";
  if (powerIterations) std::clog << " with power iterations.\n" ; else std::clog << ".\n";
  print_with_time_stamp("Starting base linear random projection", std::clog);
  Matrix localP = B * Matrix::Random(B.cols(), P.cols());         // local version to avoid later potential aliasing problem
  while (powerIterations--)
  { print_with_time_stamp("Performing B B' multiplication for power iteration", std::clog);
    Matrix R = B * (B.transpose() * localP);
    print_with_time_stamp("Performing Householder step of iterated random projection", std::clog);
    localP = Eigen::HouseholderQR<Matrix>(R).householderQ() * Matrix::Identity(P.rows(),P.cols());  // block does not work; use to get left P.cols()
  }
  std::clog << "MAIN: Check norms after Householder orthgonalization in random projection; 0'0="
	      << localP.col(0).dot(localP.col(0)) << "   0'1=" << localP.col(0).dot(localP.col(1)) << "   1'1=" << localP.col(1).dot(localP.col(1)) << std::endl;
  Matrix rB = localP.transpose() * B;
  print_with_time_stamp("Computing SVD of reduced matrix", std::clog);
  Eigen::JacobiSVD<Matrix> svd(rB, Eigen::ComputeThinU|Eigen::ComputeThinV);
  Matrix U = svd.matrixU()  ;   // nProjections x nProjections
  P = localP * U;               // beware alias here
  print_with_time_stamp("Completed random projection", std::clog);
}


void
Helper::fill_random_projection(Matrix &P, Vocabulary::SparseMatrix const& M, Vector const& leftWeights, Vector const& rightWeights, int power)
{
  std::cerr << "\n\nHLPR: *********  Using deprecated version of fill_random_projection!!!!  **********\n\n";
  assert (M.rows() == P.rows());
  Matrix R;
  Vocabulary::SparseMatrix wMw = M;
  if (0 < leftWeights.size())
  { std::clog << "MAIN: Weighting on left side in random projection.\n";
    wMw = leftWeights.asDiagonal() * wMw;
  }
  if (0 < rightWeights.size())
  { std::clog << "MAIN: Weighting on right side in random projection.\n";
    wMw = wMw * rightWeights.asDiagonal();  
  }
  R = wMw * Matrix::Random(wMw.cols(), P.cols());    
  P = Eigen::HouseholderQR<Matrix>(R).householderQ() * Matrix::Identity(P.rows(),P.cols());  // block does not work; use to get left P.cols()
  if (power > 0)
  { Vocabulary::SparseMatrix MMt = wMw * wMw.transpose();
    while (power--)
    { R = MMt * P;
      P = Eigen::HouseholderQR<Matrix>(R).householderQ() * Matrix::Identity(P.rows(),P.cols());
    }
  }
  std::clog << "MAIN: Checking norms of leading terms in random projection; 0'0="
	    << P.col(0).dot(P.col(0)) << "   0'1=" << P.col(0).dot(P.col(1)) << "   1'1=" << P.col(1).dot(P.col(1)) << std::endl;
}


void
Helper::fill_left_SVD(Matrix &U, Vocabulary::SparseMatrix const& M)
{
  Eigen::JacobiSVD<Matrix> svd(M, Eigen::ComputeThinU);
  U = svd.matrixU().leftCols(U.cols());
  // Matrix V = svd.matrixV() * Matrix::Identity(M.rows(), nProjections);
  Vector s = svd.singularValues();
  std::clog << "MAIN: Found SVD. Leading singular values are " << s.transpose().head(20) << endl;
}
