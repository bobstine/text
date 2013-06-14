#include "k_means.h"
#include "eigen_utils.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>

const int nRows = 100;
const int nCols =  12;

using std::cout;
using std::endl;
using std::string;

int main (void)
{
  // data vary from [-10,10] after randomize seed
  srand(time(NULL));
  Eigen::MatrixXf baseData = 10.0 * Eigen::MatrixXf::Random(nRows,nCols+1);

  
  const int nClusters = 4;
  const int maxIter = 10;
  
  if (false)   // l2 cluster
  {
    Eigen::MatrixXf data = baseData;
    bool useL2 = true;
    bool scale = true;
    
    for(int i=0; i<nRows; ++i)
      data.row(i).array() += (i % nClusters);
    Eigen::VectorXi wts = Eigen::VectorXi::Ones(nRows);
    KMeansClusters clusters(data.leftCols(nCols), wts, useL2, scale, nClusters, maxIter);
    clusters.print_to_stream(std::cout);
  }

  if (true)        // cosine centers
  {
    Eigen::MatrixXf data = baseData;
    bool useL2 = false;
    bool scale = false;
    
    float shift = 10;
    for(int i=0; i<nRows; ++i)
    { int offset = i % nClusters;
      for (int j=offset; j<nCols; j += nClusters)
        data(i,j) += shift;
    }
    //    std::cout << data << std::endl;
    Eigen::VectorXi wts = Eigen::VectorXi::Ones(nRows);
    const int nToFind = 7;
    KMeansClusters clusters(data.leftCols(nCols), wts, useL2, scale, nToFind, maxIter);
    clusters.print_to_stream(std::cout);
  }
}
