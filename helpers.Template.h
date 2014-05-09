#include "helpers.h"

#include "regression.h"

#include <iostream>
#include <fstream>


template <class M>
void
Helper::calculate_sequence_r2 (Eigen::VectorXd const& Y, Eigen::VectorXd tokenCount, int degree, bool reverse, M const& W,
			       Vocabulary const& vocab, int nToFit, string file)
{
  std::clog << "HLPR: Top of calculate_sequence_r2 loop; fitting " << nToFit << " word regressors.\n";
  std::ofstream os(file);
  if ((file.size()>0) && (!os)) std::clog << "MAIN: Could not open file " << file << " for writing r2 sequence.\n";
  if (os)  os << "Type  r2 adjR2  RSS q AICc\n";
  LinearRegression regr("log price", Y, 0);
  if (tokenCount.size() > 0)
  { std::clog << "HLPR: Calculate_sequence_r2 adjusted for total word count.\n";
    Eigen::VectorXd x = tokenCount;
    for(int d=1; d<=degree; ++d)
    { std::string name = "length^" + std::to_string(d);
      FStatistic f = regr.f_test_predictor(name,x);
      std::clog << "MAIN: F stat for adding " << name << " is " << f.f_stat() << std::endl;
      if(f.f_stat() > 0.000001) regr.add_predictors();
      x = x.cwiseProduct(tokenCount);
    }
    if (os)
      os << "Tokens " << regr.r_squared() << " " << regr.adj_r_squared() << " " << regr.residual_ss() << " " << regr.q() << " " << regr.aic_c() << std::endl;
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
    if(f.f_stat() > 0.00001) regr.add_predictors();
    else std::clog << "MAIN: F = " << f.f_stat() << " for word j=" << index << ", type=" << tv[index] << " is (near) singular in sequence r2 and skipped.\n";
    if (os)
      os << tv[index] << " " << regr.r_squared() << " " << regr.adj_r_squared() << " " << regr.residual_ss() << " " << regr.q() <<" "<< regr.aic_c() <<std::endl;
  }
  std::clog << "MAIN: Regression on W completed with results written to " << file << ".\n";
}
