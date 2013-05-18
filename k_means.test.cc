#include "k_means.h"
#include "eigen_utils.h"

#include <iostream>
#include <fstream>

const int nRows = 100;
const int nCols =   3;

using std::cout;
using std::endl;

int main (void)
{
  Eigen::MatrixXf data = Eigen::MatrixXf::Random(nRows,nCols+1);

  cout << "TEST: Random data is \n" << data << endl;

  int nClusters = 3;
  KMeansClusters clusters(data.topLeftCorner(nRows,nCols),nClusters);
  clusters.print_to_stream(std::cout);

  std::vector<int> ids = clusters.cluster_tags();
  for (int i=0; i<nRows; ++i)
    data(i,nCols) = ids[i];

  write_matrix_to_file("/Users/bob/Desktop/cluster_test.txt", data);
}
