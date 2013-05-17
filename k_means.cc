#include "k_means.h"

typedef Eigen::MatrixXf Matrix;
typedef Eigen::VectorXf Vector;

int
closest_cluster (Eigen::RowVectorXf const&x, Matrix const& centers)
{
  float minD = (x - centers.row(0)).squaredNorm();
  int   minI = 0;
  for (int i=1; i<centers.rows(); ++i)
  { float d = (x - centers.row(i)).squaredNorm();
    if (d < minD) { minD = d; minI = i; }
  }
  return minI;
}


ClusterMap
k_means_cluster_map(Matrix const& data, int nClusters)
{
  // init cluster centers to top rows
  Matrix centers = data.topLeftCorner(nClusters, data.cols());

  std::vector<int> clusters (data.rows());

  // assign to clusters
  for(int i=0; i<data.rows(); ++i)
    clusters[i] = closest_cluster(data.row(i), centers);

  // find new cluster centers
  centers = Matrix::Zero(centers.rows(), centers.cols());
  std::vector<int> counts (nClusters);
  
  for(int i=0; i<data.rows(); ++i)
  { centers.row(clusters[i]) += data.row(i);
    ++counts[clusters[i]];
  }
  for(int i=0; i<centers.rows(); ++i)
    centers.row(i).array() /= counts[i];

  // repeat <<<<<<<<<<--------------

  // return as map of vectors
  std::map<int, std::vector<int>> clusterMap;
  for(int i=0; i<(int)clusters.size(); ++i)
    clusterMap[clusters[i]].push_back(i);
  return clusterMap;
}
  

