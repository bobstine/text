#ifndef _K_MEANS_H_
#define _K_MEANS_H_

#include <Eigen/Core>

#include <functional>
#include <vector>
#include <map>
#include <iostream>


namespace k_means
{
  typedef Eigen::RowVectorXf RowVector;
  typedef float (*Distance)(RowVector const&, RowVector const&);
  typedef RowVector (*Renorm) (RowVector const&);
  
  inline float
    l2_distance(RowVector const& a, RowVector const& b)
  {
    return (a.array() - b.array()).matrix().squaredNorm();
  }
  
  inline RowVector
    identity(RowVector const&x)
    {
      return x;
    }
  
  inline float
    cosine_distance(RowVector const& a, RowVector const& b)
  {
    float dp = a.dot(b);
    return (dp < 0.0) ? -dp : dp;
  }
  
  inline RowVector
    two_balls(RowVector const&x)
    {
      int n = x.size()/2;
      assert (x.size() == 2*n);
      RowVector y = RowVector::Zero(2*n);
      y.head(n) = x.head(n)/x.head(n).norm();
      y.tail(n) = x.tail(n)/x.tail(n).norm();
      return y;
    }
}

class KMeansClusters
{
  typedef Eigen::MatrixXf                  Matrix;
  typedef Eigen::VectorXf                  Vector;
  typedef Eigen::VectorXi                  IntVector;
  typedef Eigen::RowVectorXf               RowVector;
  typedef std::map<int, std::vector<int>>  Map;
  typedef float (*Distance)(RowVector const&, RowVector const&);
  typedef RowVector (*Renorm) (RowVector const&);
  
  
  Matrix const&      mData;
  IntVector const&   mWeights;
  Distance           mDist;
  Renorm             mRenorm;
  int                mNClusters;
  Matrix             mClusterCenters;
  std::vector<int>   mClusterTags;

 public:
  KMeansClusters (Matrix const& data, IntVector const& wts, Distance f, Renorm g, int nClusters, int maxIterations = 10)
    : mData(data), mWeights(wts), mDist(f), mRenorm(g), mNClusters(nClusters),
      mClusterCenters(Matrix::Zero(nClusters,data.cols())), mClusterTags(data.rows())
    { find_clusters(maxIterations); }

  Map              cluster_map  ()  const;
  std::vector<int> cluster_tags ()  const { return mClusterTags; }

  void             find_clusters(int maxIterations);

  void             print_to_stream (std::ostream& os) const;

 private:
  int    closest_cluster           (RowVector const& r, Matrix const& m) const;
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
