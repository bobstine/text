#include "k_means.h"
#include <utility>
#include <iostream>


template <class D>
int
KMeansClusters<D>::closest_cluster (RowVector const& r, Matrix const& m) const
{
  typename D::result_type minD = mDist(r,m.row(0));
  int   minI = 0;
  for (int i=1; i<m.rows(); ++i)
  { typename D::result_type d = mDist(r,m.row(i));
    if (d < minD) { minD = d; minI = i; }
  }
  return minI;
}

template <class D>
void
KMeansClusters<D>::find_clusters(int maxIterations)
{
  // init new cluster centers to top rows
  Matrix newCenters = Matrix (mNClusters,mData.cols());
  newCenters = mData.topLeftCorner(mNClusters, mData.cols());
  // set 'old' center to zero so initial change is infinite
  Matrix oldCenters = Matrix::Zero       (mNClusters, mData.cols());

  std::vector<int> counts   (mNClusters);
  int itCount = 0;
  Matrix *pNew = &newCenters;
  Matrix *pOld = &oldCenters;
  while ((0.001 < relative_squared_distance(*pNew,*pOld)) && (++itCount < maxIterations))
  { // assign cases to clusters
    for(int i=0; i<mData.rows(); ++i)
      mClusterTags[i] = closest_cluster(mData.row(i), *pNew);
    std::swap(pNew,pOld);
    // calculate new centers of clusters
    *pNew = Matrix::Zero(mNClusters, mData.cols());
    for(int i=0; i<mNClusters; ++i) counts[i]=0;
    for(int i=0; i<mData.rows(); ++i)
    { (*pNew).row(mClusterTags[i]) += mData.row(i) * mWeights(i);
      counts[mClusterTags[i]] += mWeights(i);
    }
    for(int i=0; i<mNClusters; ++i)
    { (*pNew).row(i).array() /= counts[i];
      renorm_center((*pNew).row(i));
    }
    // std::clog << "Center counts   \n" ;
    // for(int i=0;i<mNClusters; ++i) std::clog << counts[i] << " "; std::clog << std::endl;
    // std::clog << "Centers at step " << itCount << ":\n" << *pNew << std::endl;
  }
  mClusterCenters = *pNew;
}


template<>
void
KMeansClusters<L2Distance>::renorm_center (RowVector &) const
{
  // do nothing
}

template<>
void
KMeansClusters<CosineDistance>::renorm_center (RowVector &r) const
{
  // weird renorm for split ball
  int j (r.cols()/2);
  r.leftCols(j)  /= r.leftCols(j).norm();
  r.rightCols(j) /= r.rightCols(j).norm();
}


template <class D>
typename KMeansClusters<D>::Map
KMeansClusters<D>::cluster_map() const
{
  KMeansClusters::Map m;
  for(int i=0; i<(int)mClusterTags.size(); ++i)
    m[mClusterTags[i]].push_back(i);
  return m;
}

template <class D>
double
KMeansClusters<D>::relative_squared_distance (Matrix const& newCenters, Matrix const& oldCenters) const
{
  // avg squared distance relative to avg squared size
  double dist = (newCenters.array() - oldCenters.array()).matrix().squaredNorm()/oldCenters.squaredNorm();
  std::clog << "KMCS: distance = " << dist << std::endl;
  return dist;
}
  
template <class D>
void
KMeansClusters<D>::print_to_stream (std::ostream& os) const
{
  Map m (cluster_map());
  for(int i=0; i<mNClusters; ++i)
  { os << "Cluster " << i << ": ";
    for(int j=0; j<(int)m[i].size(); ++j)
      os << " " << m[i][j];
    os << std::endl;
  }
}
