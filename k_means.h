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
  
  inline float
    l2_distance(RowVector const& a, RowVector const& b)
  {
    return (a.array() - b.array()).matrix().squaredNorm();
  }
  
  inline float
    cosine_distance(RowVector const& a, RowVector const& b)
  {
    return 1.0 - a.dot(b);
  }
}

class KMeansClusters
{
 public:
  typedef Eigen::MatrixXf                  Matrix;
  typedef Eigen::VectorXf                  Vector;
  typedef Eigen::VectorXi                  IntVector;
  typedef Eigen::RowVectorXf               RowVector;
  typedef std::map<int, std::vector<int>>  Map;
  typedef k_means::Distance                Distance;
  typedef std::vector<int>::const_iterator Iterator;

 private:
  Matrix             mData;
  IntVector const&   mWeights;
  bool               mUseL2;
  Distance           mDist;
  bool               mScaleData;
  int                mNClusters;
  Matrix             mClusterCenters;
  std::vector<int>   mClusterTags;

 public:
  
  KMeansClusters (Matrix const& data, IntVector const& wts, bool l2, bool scaleData, int nClusters, int maxIterations = 10)
    : mData(data), mWeights(wts),
      mUseL2(l2), mDist(l2 ? k_means::l2_distance : k_means::cosine_distance ), mScaleData(scaleData),
      mNClusters(nClusters), mClusterCenters(Matrix::Zero(nClusters,data.cols())),
      mClusterTags(data.rows())
    { prepare_data(&mData); find_clusters(maxIterations); }

  Map              cluster_map  ()                                       const;  // map of vectors for indices in each cluster
  
  Iterator         cluster_tags_begin ()                                 const { return mClusterTags.cbegin(); }
  Iterator         cluster_tags_end ()                                   const { return mClusterTags.cend(); }

  void             print_to_stream (std::ostream& os)                    const;

  std::vector<int> assign_to_clusters (Matrix *data)                     const;  // note modifies argument
  
 private:
  void   find_clusters(int maxIterations);
  void   prepare_data              (Matrix *data)                        const;  // note in place argument
  int    closest_cluster           (RowVector const& r, Matrix const& m) const;
  double relative_squared_distance (Matrix const& a, Matrix const& b)    const;
};


inline
std::ostream&
operator<< (std::ostream& os, KMeansClusters const& x)
{
  x.print_to_stream(os);
  return os;
}

#endif
