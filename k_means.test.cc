#include "k_means.h"
#include "eigen_utils.h"

#include <iostream>
#include <fstream>
#include <time.h>

const int nRows = 200;
const int nCols =  10;

using std::cout;
using std::endl;

int main (void)
{
  // data vary from [-10,10] after randomize seed
  srand(time(NULL));
  Eigen::MatrixXf data = Eigen::MatrixXf::Random(nRows,nCols+1);
  data *= 10.0;
  
  int nClusters = 4;
  for(int i=0; i<nRows; ++i)
    data.row(i).array() += (i % 4);
  
  //  cout << "TEST: Random data is \n" << data << endl;

  Eigen::VectorXi wts = Eigen::VectorXi::Ones(nRows);
  KMeansClusters clusters(data.leftCols(nCols),wts,k_means::l2_distance,k_means::identity,nClusters,11);
  clusters.print_to_stream(std::cout);

  std::vector<int> ids = clusters.cluster_tags();
  for (int i=0; i<nRows; ++i)
    data(i,nCols) = ids[i];

  write_matrix_to_file("/Users/bob/Desktop/cluster_test.txt", data);
}
