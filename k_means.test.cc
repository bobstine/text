#include "k_means.h"

#include <iostream>

const int nRows = 10;
const int nCols =  3;

using std::cout;
using std::endl;

int main (void)
{
  
  Eigen::MatrixXf data = Eigen::MatrixXf::Random(nRows,nCols);

  cout << "TEST: Random data is \n" << data << endl;

  int nClusters = 3;
  
  ClusterMap clusters = k_means_cluster_map(data,nClusters);


  for(int i=0; i<nClusters; ++i)
  { cout << "Cluster " << i << ": ";
    for(int j=0; j<(int)clusters[i].size(); ++j)
      cout << " " << clusters[i][j];
    cout << endl;
  }
}
