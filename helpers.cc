#include "helpers.h"

#include "vocabulary.h"
#include "eigenword_dictionary.h"
#include "read_utils.h"
#include "file_utils.h"
#include "regex.h"
#include "regression.h"

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
Helper::calculate_sequence_r2 (Eigen::VectorXd const& Y, Eigen::VectorXd tokenCount, string xLabel, Eigen::MatrixXf X, string file)
{
  const int k = X.cols();
  std::clog << "MAIN: Top of calculate_sequence_r2 loop; fitting " << k << " regressors labeled " << xLabel << endl;
  std::ofstream os(file);
  if ((file.size()>0) && (!os)) std::clog << "MAIN: Could not open file " << file << " for writing r2 sequence.\n";
  if (os)  os << "Type  r2  RSS AICc\n";
  LinearRegression regr("log price", Y, 0);
  if (tokenCount.size() > 0)
  { std::clog << "MAIN: Calculate_sequence_r2 adjusted for total word count.\n";
    FStatistic f = regr.f_test_predictor("m",tokenCount);
    std::clog << "MAIN: F stat for adding token lengths is " << f.f_stat() << endl;
    if(f.f_stat() > 0.0001) regr.add_predictors();
    if (os)
      os << "nTokens " << regr.r_squared() << " " << regr.residual_ss() << " " << regr.aic_c() << endl;
  }    
  Eigen::VectorXf ind = Eigen::VectorXf::Zero(k);
  for (int j=0; j<k; ++j)
  { std::stringstream ss;
    ss << j;
    string name = xLabel + ss.str();
    Eigen::VectorXd xj(X.rows(),1);                 // transfer to double
    for(int i=0; i<xj.size(); ++i) xj(i) = X.coeff(i,j);
    FStatistic f=regr.f_test_predictor(name, xj);
    if(f.f_stat() > 0.0001) regr.add_predictors();
    else std::clog << "MAIN: F = " << f.f_stat() << " for " << name << " is (near) singular in sequence r2 and skipped.\n";
    if (os)
      os << name << " " << regr.r_squared() << " " << regr.residual_ss() << " " << regr.aic_c() << endl;
  }
  std::clog << "MAIN: Regression on W completed with results written to " << file << ".\n";
}


void
Helper::calculate_sequence_r2 (Eigen::VectorXd const& Y, Eigen::VectorXd tokenCount, bool reverse, Vocabulary::SparseMatrix const& W,
		       Vocabulary const& vocab, int nToFit, string file)
{
  std::clog << "MAIN: Top of calculate_sequence_r2 loop; fitting " << nToFit << " word regressors.\n";
  std::ofstream os(file);
  if ((file.size()>0) && (!os)) std::clog << "MAIN: Could not open file " << file << " for writing r2 sequence.\n";
  if (os)  os << "Type  r2  RSS AICc\n";
  LinearRegression regr("log price", Y, 0);
  if (tokenCount.size() > 0)
  { std::clog << "MAIN: Calculate_sequence_r2 adjusted for total word count.\n";
    FStatistic f = regr.f_test_predictor("m",tokenCount);
    std::clog << "MAIN: F stat for adding token lengths is " << f.f_stat() << endl;
    if(f.f_stat() > 0.0001) regr.add_predictors();
    if (os)
      os << "nTokens " << regr.r_squared() << " " << regr.residual_ss() << " " << regr.aic_c() << endl;
  }    
  const int k = W.cols();
  Vocabulary::TypeVector tv = vocab.types();
  Eigen::VectorXf ind = Eigen::VectorXf::Zero(k);
  for (int j=0; j<nToFit; ++j)
  { int index;
    if(reverse) index = (nToFit-1) - j;
    else        index = j;
    Eigen::VectorXf x(W.rows());
    ind(index) = 1;
    x = W * ind;                      // elaborate way to get column out of sparse matrix
    ind(index) = 0;
    Eigen::VectorXd X(x.size(),1);    // transfer to double
    for(int i=0; i<x.size(); ++i) X(i) = x(i);
    FStatistic f=regr.f_test_predictor("xx", X);
    if(f.f_stat() > 0.0001) regr.add_predictors();
    else std::clog << "MAIN: F = " << f.f_stat() << " for word j=" << index << ", type=" << tv[index] << " is (near) singular in sequence r2 and skipped.\n";
    if (os)
      os << tv[index] << " " << regr.r_squared() << " " << regr.residual_ss() << " " << regr.aic_c() << endl;
  }
  std::clog << "MAIN: Regression on W completed with results written to " << file << ".\n";
}

void
Helper::fill_random_projection(Matrix &P, Vocabulary::SparseMatrix const& M, Vector const& wts, int power)
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


void
Helper::fill_left_SVD(Matrix &U, Vocabulary::SparseMatrix const& M)
{
  Eigen::JacobiSVD<Matrix> svd(M, Eigen::ComputeThinU);
  U = svd.matrixU().leftCols(U.cols());
  // Matrix V = svd.matrixV() * Matrix::Identity(M.rows(), nProjections);
  Vector s = svd.singularValues();
  std::clog << "MAIN: Found SVD. Leading singular values are " << s.transpose().head(20) << endl;
}