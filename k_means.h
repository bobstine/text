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
  typedef std::map<int, std::vector<int>>  ClusterMap;
  typedef k_means::Distance                Distance;
  typedef std::vector<int>                 IntegerVector; 
  typedef std::vector<int>::const_iterator Iterator;

 private:
  Matrix                   mData;
  IntVector const&         mWeights;
  bool                     mUseL2;
  Distance                 mDist;
  bool                     mScaleData;
  int                      mNClusters;
  Matrix                   mClusterCenters;
  IntegerVector            mDataClusterIndex;    // cluster index for each data row

 public:
  
  KMeansClusters (Matrix const& data, IntVector const& wts, bool l2, bool scaleData, int nClusters, int maxIterations = 10)
    : mData(data), mWeights(wts), mUseL2(l2), mDist(l2 ? k_means::l2_distance : k_means::cosine_distance ), mScaleData(scaleData),
      mNClusters(nClusters), mClusterCenters(Matrix::Zero(nClusters,data.cols())), mDataClusterIndex(data.rows())
    { prepare_data(&mData); find_clusters(maxIterations);  }

  int            number_of_cases()                                       const { return mData.rows(); }
  int            number_of_clusters()                                    const { return mNClusters; }
  
  ClusterMap     cluster_map  ()                                         const;  // map of vectors for indices in each cluster
  
  // these iterate over input data (types in text application)
  Iterator       data_cluster_index_begin ()                             const { return mDataClusterIndex.cbegin(); }
  Iterator       data_cluster_index_end ()                               const { return mDataClusterIndex.cend(); }

  IntegerVector  assign_cluster_indices (Matrix *data)                   const;  // note modifies argument

  void           print_to_stream (std::ostream& os)                      const;
  
 private:
  void   prepare_data              (Matrix *data)                        const;  // note in place argument
  void   find_clusters             (int maxIterations);
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
