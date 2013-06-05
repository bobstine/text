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

using std::string;

template <class I, class C>
class IndexedStringIterator: public std::iterator<std::forward_iterator_tag, string>
{
  I         mIndexIt;
  C const&  mStrings;
  
public:
  IndexedStringIterator  (I it, C const& labels)  : mIndexIt (it), mStrings(labels) { }
  
  bool          operator==(IndexedStringIterator const& it) const  { return (mIndexIt == it.mIndexIt) && (mStrings == it.mStrings); }
  bool          operator!=(IndexedStringIterator const& it) const  { return (mIndexIt != it.mIndexIt) || (mStrings != it.mStrings); }
  
  IndexedStringIterator&  operator++()                             { ++mIndexIt; return *this; }
  string                  operator*()                       const  { return mStrings[*mIndexIt]; }
};


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
  typedef IndexedStringIterator<std::vector<int>::const_iterator, std::vector<string>> ISIterator;
  
 private:
  Matrix                   mData;
  IntVector const&         mWeights;
  bool                     mUseL2;
  Distance                 mDist;
  bool                     mScaleData;
  int                      mNClusters;
  Matrix                   mClusterCenters;
  std::vector<std::string> mClusterLabels;
  std::vector<int>         mDataClusterIndex;       // for the original cases which are types

 public:
  
  KMeansClusters (Matrix const& data, IntVector const& wts, std::vector<std::string> const& caseLabels, bool l2, bool scaleData, int nClusters, int maxIterations = 10)
    : mData(data), mWeights(wts),
      mUseL2(l2), mDist(l2 ? k_means::l2_distance : k_means::cosine_distance ), mScaleData(scaleData),
      mNClusters(nClusters), mClusterCenters(Matrix::Zero(nClusters,data.cols())), mClusterLabels(nClusters),
      mDataClusterIndex(data.rows())
	{ prepare_data(&mData); find_clusters(maxIterations); label_clusters(caseLabels); }

  Map              cluster_map  ()                                       const;  // map of vectors for indices in each cluster

  // these iterate over input data (types in text application)
  Iterator         cluster_index_begin ()                                const { return mDataClusterIndex.cbegin(); }
  Iterator         cluster_index_end ()                                  const { return mDataClusterIndex.cend(); }
  ISIterator       cluster_tag_begin()                                   const { return ISIterator(mDataClusterIndex.cbegin(), mClusterLabels); }
  ISIterator       cluster_tag_end()                                     const { return ISIterator(mDataClusterIndex.cend(), mClusterLabels); }
  

  void             print_to_stream (std::ostream& os)                    const;

  std::vector<int> assign_to_clusters (Matrix *data)                     const;  // note modifies argument
  
 private:
  void   prepare_data              (Matrix *data)                        const;  // note in place argument
  void   find_clusters             (int maxIterations);
  void   label_clusters            (std::vector<std::string> const& labels);
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
