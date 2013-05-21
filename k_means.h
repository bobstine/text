#ifndef _K_MEANS_H_
#define _K_MEANS_H_

#include <Eigen/Core>

#include <vector>
#include <map>
#include <iostream>

class KMeansClusters
{
  typedef Eigen::MatrixXf                  Matrix;
  typedef std::map<int, std::vector<int>>  Map;
  typedef Eigen::RowVectorXf               RowVector;

  Matrix const&      mData;
  int                mNClusters;
  Matrix             mClusterCenters;
  std::vector<int>   mClusterTags;

 public:
  KMeansClusters (Matrix const& data, int nClusters, int maxIterations = 10)
    : mData(data), mNClusters(nClusters), mClusterCenters(Matrix::Zero(nClusters,data.cols())), mClusterTags(data.rows())
    { find_clusters(maxIterations); }

  Map              cluster_map  ()  const;
  std::vector<int> cluster_tags ()  const;

  void             find_clusters(int maxIterations);

  void             print_to_stream (std::ostream& os) const;

 private:
  int    closest_cluster (RowVector const& r, Matrix const& m) const;
  double relative_squared_distance (Matrix const& a, Matrix const& b) const;
};


inline
std::ostream&
operator<< (std::ostream& os, KMeansClusters const& x)
{
  x.print_to_stream(os);
  return os;
}

#endif
