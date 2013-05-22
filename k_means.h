#ifndef _K_MEANS_H_
#define _K_MEANS_H_

#include <Eigen/Core>

#include <functional>
#include <vector>
#include <map>
#include <iostream>


class L2Distance: public std::binary_function<Eigen::RowVectorXf const&, Eigen::RowVectorXf const&, float>
{
  typedef Eigen::RowVectorXf RowVector;
 public:
  float
    operator()(RowVector const& a, RowVector const& b) const
  { return (a.array() - b.array()).matrix().squaredNorm(); }
};

class CosineDistance: public std::binary_function<Eigen::RowVectorXf const&, Eigen::RowVectorXf const&, float>
{
  typedef Eigen::RowVectorXf RowVector;
 public:
  float
    operator()(RowVector const& a, RowVector const& b) const
  { return a.dot(b); }
};


template <class Distance>
class KMeansClusters
{
  typedef Eigen::MatrixXf                  Matrix;
  typedef Eigen::VectorXf                  Vector;
  typedef Eigen::VectorXi                  IntVector;
  typedef Eigen::RowVectorXf               RowVector;
  typedef std::map<int, std::vector<int>>  Map;

  Matrix const&      mData;
  IntVector const&   mWeights;
  Distance           mDist;
  int                mNClusters;
  Matrix             mClusterCenters;
  std::vector<int>   mClusterTags;

 public:
    KMeansClusters (Matrix const& data, IntVector const& wts, Distance const& f, int nClusters, int maxIterations = 10)
      : mData(data), mWeights(wts), mDist(f), mNClusters(nClusters), mClusterCenters(Matrix::Zero(nClusters,data.cols())), mClusterTags(data.rows())
    { find_clusters(maxIterations); }

  Map              cluster_map  ()  const;
  std::vector<int> cluster_tags ()  const { return mClusterTags; }

  void             find_clusters(int maxIterations);

  void             print_to_stream (std::ostream& os) const;

 private:
  int    closest_cluster           (RowVector const& r, Matrix const& m) const;
  double relative_squared_distance (Matrix const& a, Matrix const& b) const;
  void   renorm_center             (RowVector &pR) const;  // evil reference
};


template <class D>
inline
std::ostream&
operator<< (std::ostream& os, KMeansClusters<D> const& x)
{
  x.print_to_stream(os);
  return os;
}

#endif
