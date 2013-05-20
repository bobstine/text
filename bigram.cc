/*
  R Code for check/comparison of SVD
 
  X <- read.table("/Users/bob/C/text/random_projection.txt")

  udv <- svd (X)

  udv$d

  for(j in 1:10) {
  cat(j,"  ",any(is.nan(X[,j])),"  \n")
  cat(which (is.nan(X[,j])), "  \n");
  }
*/

#include "bigram.h"
#include "eigen_utils.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SVD>

#include <iostream>
#include <sstream>

#include <time.h>

typedef Eigen::VectorXf                   Vector;
typedef Eigen::MatrixXf                   Matrix;

typedef Eigen::SparseMatrix<float,Eigen::RowMajor> SparseMatrix;


void
print_time(std::ostream& os, std::string const& s, clock_t const& start, clock_t const& stop)
{
  clock_t diff (stop-start);
  os << "TIME: " << s << " [" << ((float)diff)/CLOCKS_PER_SEC << " sec]\n";
}

std::ostream&
operator<< (std::ostream& os, Vector const& x)
{
  os << "{";
  for (int i=0; i<x.size(); ++i)
    os << " " << x[i];
  os << "}";
  return os;
}

  
int main(int, char **)
{
  const int nProjections = 30;

  using std::string;
  
  std::ostringstream ss;
  
  // read tokens as first string in file, one item per line
  std::list<string>    tokenList;
  std::list<string>    goldPOSList;
  // maps hold token assignments, counts
  std::map<string, std::map<string,int>> tokenPOSMap;
  std::map<string,int> tokenMap;
  
  clock_t startTime;

  startTime = clock();
  while (std::cin)
  { std::string token;
    std::cin >> token;
    if (token.empty()) break;   // end of file
    tokenList.push_back(token);
    ++tokenMap[token];
    std::string pos;
    std::cin >> pos;
    goldPOSList.push_back(pos);
    ++tokenPOSMap[token][pos];
  }
  ss << "Read " << tokenMap.size() << " unique tokens from input of length " << tokenList.size() << ".";
  print_time(std::clog, ss.str(), startTime, clock());
  ss.str("");
  

  const int nTokens ((int)tokenMap.size());


  // then assign unique integer id to tokens
  startTime = clock();
  TokenManager tokens (tokenMap);
  print_time(std::clog, "Sort tokens and assign IDs in TokenManager.", startTime, clock());

  // build the sparse bigram array using integer id for words
  startTime = clock();
  SparseMatrix B(nTokens,nTokens);
  {
    typedef std::pair<int,int> Index;
    std::map<Index,int> bgram;
    Index i = std::make_pair(0,tokens[tokenList.front()]);
    for(auto it = ++tokenList.begin(); it!=tokenList.end(); ++it)
    { i.first = i.second;
      i.second= tokens[*it];
      ++bgram[i];
    }
    typedef Eigen::Triplet<float> T;
    std::list<T> triplets (nTokens);
    for(auto it = bgram.begin(); it != bgram.end(); ++it)
      triplets.push_back(T(it->first.first, it->first.second, it->second));
    B.setFromTriplets(triplets.begin(), triplets.end());
  }
  ss << "Init sparse bigram B[" << B.rows() << "x" << B.cols() << "] from map.";
  print_time(std::clog, ss.str(), startTime, clock());
  ss.str("");
  
  // test by adding up all of the elements in B
  ss << "Sums of first rows of B are ";
  startTime = clock();
  {
    Vector one  (Vector::Constant(nTokens,1.0));
    Vector sums (Vector::Zero    (nTokens));
    sums = B * one;
    for(int i=0; i<5; ++i)
      ss << "[" << i << "]=" << sums[i] << "   ";
    ss << " with +/+/B = " << sums.sum() << ".";
  }
  print_time(std::clog, ss.str(), startTime, clock());
  ss.str("");
  
  // Random projection of the array
  startTime = clock();
  Matrix RP (B.rows(),nProjections);
  {
    Vector norms(B.rows());
    for (int i=0; i<B.rows(); ++i)
    { float ss = B.row(i).norm();
      if (0 == ss)
      { ss = 1;
	std::clog << "MAIN: Zero norm for B[" << i << "] with token ----->" << tokens[i] << "<----- with count " << tokenMap[tokens[i]] << std::endl;
      }
      else
        norms(i) = 1.0/ss;
    }
    Matrix R = Matrix::Random(B.cols(),nProjections);
    RP = norms.asDiagonal() * (B * R);
  }
  ss << "Compute scaled random projection RP[" << RP.rows() << "x" << RP.cols() << "].";
  print_time(std::clog, ss.str(), startTime, clock());
  ss.str("");
  std::clog << " RP matrix : \n"
	    << RP.topLeftCorner(5,nProjections) << "\n  ...\n"
	    << RP.bottomLeftCorner(5,nProjections) << std::endl;

  // cluster tokens using k-means
  const int nClusters = 10;
  KMeansClusters clusters (RP, 20);

  // write data and clusters
  {
    std::ios_base::openmode mode = (append) ? std::ios_base::app : std::ios_base::trunc;
    std::ofstream output (fileName.c_str(), mode);
    if (! output)
      std::cerr << "ERROR: Cannot open output text file for writing matrix to file " << fileName << std::endl;
    else
    { std::vector<int> group = cluster.tags();
      for (i=0; i<RP.rows(); ++i)
      { output << RP.row(i) << " " << group(i) << " " << 
  }




  // SVD of random projection array
  if (false)
  { startTime = clock();
    Eigen::JacobiSVD<Matrix, Eigen::HouseholderQRPreconditioner> svd(RP, Eigen::ComputeThinU);
    Matrix UD = svd.matrixU() * svd.singularValues().asDiagonal();
    // std::clog << " U matrix      : \n" << UD.topLeftCorner(10,nProjections) << std::endl;
    std::clog << "Singular values: \n" << svd.singularValues().transpose().head(nProjections) << std::endl;
    print_time(std::clog, "Compute SVD and extract U*D.", startTime, clock());
    
    // write u matrix to file
    // write_matrix_to_file("svd.ud", U);
  }
  
    // print first items in the map
  /*
    for(auto it = sortedTokenMap.begin(); it != sortedTokenMap.end(); ++it)
    std::cout << it->second << " " << -it->first << std::endl;  
  */
  return 0;
}

