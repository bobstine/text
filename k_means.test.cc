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
  // data vary from [-sigma, sigma] after randomize seed
  srand(time(NULL));
  const float sigma = 10.0;
  Eigen::MatrixXf baseData = sigma * Eigen::MatrixXf::Random(nRows,nCols+1);
  
  const int nClusters = 4;
  const int maxIter = 10;

  { // test data in row i are shifted by (i mod k)
    // if sigma is 'small enough', see (k=4) clusters with 0 4 8 12 ...   1 5 9 13 ...  etc
    Eigen::MatrixXf data = baseData;
    for(int i=0; i<nRows; ++i)
      data.row(i).array() += (i % nClusters);
    
    Eigen::VectorXi wts = Eigen::VectorXi::Ones(nRows);
    bool biDir =  true;
    bool scale = false;
    KMeansClusters clusters(data.leftCols(nCols), biDir, wts, scale, nClusters, scale, maxIter);
    clusters.print_summary_stats(std::cout);
    std::cout << "TEST: Within cluster summary statistics (n, avg dist, max dist, norm, centroid) \n"
	      << clusters.within_cluster_summary_stats()
	      << std::endl;
    clusters.print_to_stream(std::cout);
  }

}
