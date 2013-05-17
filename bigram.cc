#include "eigen_utils.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SVD>

#include <list>
#include <map>
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
  const int nProjections = 10;

  std::ostringstream ss;
  
  // read tokens as first string in file, one item per line
  std::list<std::string>    tokenList;
  std::list<std::string>    goldPOSList;
  std::map<std::string,int> tokenMap;
  
  clock_t startTime;

  startTime = clock();
  while (std::cin)
  { std::string s;
    std::cin >> s;
    tokenList.push_back(s);
    ++tokenMap[s];
    std::cin >> s;
    goldPOSList.push_back(s);
  }
  ss << "Read " << tokenMap.size() << " unique tokens from input of length " << tokenList.size();
  print_time(std::clog, ss.str(), startTime, clock());
  ss.str("");
  

  const int nTokens ((int)tokenMap.size());


  {
    using std::clog;
    using std::endl;
    
    Matrix m = Matrix::Random(3,2);
    clog << "Here is the matrix m:" << endl << m << endl;
    Eigen::JacobiSVD<Matrix> svd(m, Eigen::ComputeThinU | Eigen::ComputeThinV);
    clog << "Its singular values are:" << endl << svd.singularValues() << endl;
    clog << "Its left singular vectors are the columns of the thin U matrix:" << endl << svd.matrixU() << endl;
    clog << "Its right singular vectors are the columns of the thin V matrix:" << endl << svd.matrixV() << endl;

  }
  
  // sort map by values by inserting into multimap; note sign reversed so descends
  // then assign unique integer id to tokens
  startTime = clock();
  std::multimap<int,std::string> sortedTokenMap;
  for (auto it = tokenMap.begin(); it != tokenMap.end(); ++it)
    sortedTokenMap.insert( std::make_pair(-it->second,it->first) );
  std::map<std::string, int> tokenID;
  {
    int i=0;
    for (auto it = sortedTokenMap.begin(); it != sortedTokenMap.end(); ++it)
      tokenID[it->second] = i++;
  }
  print_time(std::clog, "Sort tokens and assign IDs", startTime, clock());

  // build the sparse bigram array
  startTime = clock();
  SparseMatrix B(nTokens,nTokens);
  {
    typedef std::pair<int,int> Index;
    std::map<Index,int> bgram;
    Index i = std::make_pair(0,tokenID[tokenList.front()]);
    for(auto it = ++tokenList.begin(); it!=tokenList.end(); ++it)
    { i.first = i.second;
      i.second= tokenID[*it];
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
    for(int i=0; i<10; ++i)
      ss << "[" << i << "]=" << sums[i] << "   ";
    ss << " with +/+/B = " << sums.sum() << std::endl;
  }
  print_time(std::clog, ss.str(), startTime, clock());
  ss.str("");
  
  // Random projection of the array
  startTime = clock();
  Matrix RP (B.rows(),nProjections);
  {
    Vector norms(B.rows());
    for (int i=0; i<B.rows(); ++i)
      norms(i) = 1.0/B.row(i).norm();
    Matrix R = Matrix::Random(B.cols(),nProjections);
    RP = norms.asDiagonal() * (B * R);
  }
  print_time(std::clog, "Compute scaled random projection", startTime, clock());
  std::clog << " RP matrix : \n"
	    << RP.topLeftCorner(10,nProjections) << "\n  ...\n"
	    << RP.bottomLeftCorner(10,nProjections) << std::endl;
  write_matrix_to_file("random_projection.txt", RP);

  // SVD of random projection array
  startTime = clock();
  Eigen::JacobiSVD<Matrix, Eigen::HouseholderQRPreconditioner> svd(RP, Eigen::ComputeThinU);
  Matrix U = svd.matrixU(); // * svd.singularValues().asDiagonal();
  std::clog << " U matrix : \n" << U.topLeftCorner(10,10) << std::endl;
  std::clog << "Singular values: \n    " << svd.singularValues().transpose().head(10) << std::endl;
  print_time(std::clog, "Compute SVD and extract U*D", startTime, clock());

  // write u matrix to file
  write_matrix_to_file("svd.ud", U);

  // print first items in the map
  for(auto it = sortedTokenMap.begin(); it != sortedTokenMap.end(); ++it)
    std::cout << it->second << " " << -it->first << std::endl;  
  
  return 0;
}

