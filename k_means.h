#ifndef _K_MEANS_H_
#define _K_MEANS_H_

#include <Eigen/Core>

#include <functional>
#include <vector>
#include <map>
#include <iostream>
#include <cmath>

class KMeansClusters
{
  
 public:
  typedef Eigen::MatrixXf                  Matrix;
  typedef Eigen::VectorXf                  Vector;
  typedef Eigen::VectorXi                  IntVector;
  typedef Eigen::RowVectorXf               RowVector;
  typedef std::map<int, std::vector<int>>  ClusterMap;
  typedef std::vector<int>                 IntegerVector; 
  typedef std::vector<int>::const_iterator Iterator;

 private:
  Matrix                   mData;                // need local copy of data since modified
  bool                     mBidirectional;       // effects normalization
  IntVector const&         mWeights;
  bool                     mScaleData;           // normalize to |x|=1 
  int                      mNClusters;
  bool                     mScaleCentroid;       // normalize to |c|=1  (if both normed, then cosine distance)
  Matrix                   mClusterCenters;
  IntegerVector            mDataClusterIndex;    // cluster index for each data row

 public:
  
  KMeansClusters (Matrix const& data, bool biDir, IntVector const& wts, bool scaleData, int nClusters, bool scaleCentroid, int maxIterations = 10)
    : mData(data), mBidirectional(biDir), mWeights(wts), mScaleData(scaleData), mNClusters(nClusters), mScaleCentroid(scaleCentroid),
      mClusterCenters(Matrix::Zero(nClusters,data.cols())), mDataClusterIndex(data.rows())
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
  void   prepare_data      (Matrix *m)                                   const;
  void   find_clusters     (int maxIterations);
  int    closest_cluster   (RowVector const& r, Matrix const& m)         const;
  float  distance          (RowVector const& a, RowVector const& b)      const  { return (a.array() - b.array()).matrix().squaredNorm(); }
  double relative_squared_distance (Matrix const& a, Matrix const& b)    const;

  float  average_within_cluster_variance(Matrix const& centroids)        const;


};


inline
std::ostream&
operator<< (std::ostream& os, KMeansClusters const& x)
{
  x.print_to_stream(os);
  return os;
}

#endif
