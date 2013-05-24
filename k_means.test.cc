#include "k_means.h"
#include "eigen_utils.h"

#include <iostream>
#include <fstream>
#include <time.h>

const int nRows = 30;
const int nCols = 12;

using std::cout;
using std::endl;

int main (void)
{
  std::vector<int> tags;
  // data vary from [-10,10] after randomize seed
  srand(time(NULL));
  Eigen::MatrixXf baseData = 10.0 * Eigen::MatrixXf::Random(nRows,nCols+1);
  
  int nClusters = 4;
  if (false)   // l2 cluster
  {
    Eigen::MatrixXf data = baseData;
    bool useL2 = true;
    bool scale = true;
    
    for(int i=0; i<nRows; ++i)
      data.row(i).array() += (i % nClusters);
    Eigen::VectorXi wts = Eigen::VectorXi::Ones(nRows);
    KMeansClusters clusters(data.leftCols(nCols),wts,useL2,scale,nClusters,11);
    clusters.print_to_stream(std::cout);
    tags = clusters.cluster_tags();
  }

  if (true)        // cosine centers
  {
    Eigen::MatrixXf data = baseData;
    bool useL2 = false;
    bool scale = false;
    
    float shift = 100;
    for(int i=0; i<nRows; ++i)
    { int offset = i % nClusters;
      for (int j=offset; j<nCols; j += nClusters)
        data(i,j) += shift;
    }
    //    std::cout << data << std::endl;
    Eigen::VectorXi wts = Eigen::VectorXi::Ones(nRows);
    KMeansClusters clusters(data.leftCols(nCols),wts, useL2, scale, nClusters,11);
    clusters.print_to_stream(std::cout);
    tags = clusters.cluster_tags();
  }
  
  if (false) // write to file
  { for (int i=0; i<nRows; ++i)
      baseData(i,nCols) = tags[i];
    write_matrix_to_file("/Users/bob/Desktop/cluster_test.txt", baseData);
  }
}
