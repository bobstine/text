#include "k_means.h"
#include "order.h"

#include <utility>
#include <iostream>

const float sqrt2 = 1.41421356;

inline int min(int x, int y) { return (x < y) ? x : y; }


std::vector<int>
KMeansClusters::assign_cluster_indices (Matrix *data) const
{
  prepare_data(data);
  std::vector<int> intTags;
  for (int i=0; i<data->rows(); ++i)
    intTags.push_back(closest_cluster(data->row(i), mClusterCenters));
  return intTags;
}

void
KMeansClusters::prepare_data(Matrix *data) const
{
  if(mUseL2)                               // L2 norm on unit ball is optionally scaled
  { if (mScaleData)
    { for (int i=0; i<data->rows(); ++i)
	data->row(i) /= data->row(i).norm();
    }
  }
  else                                     // norm on split-ball, always scaled for cosine
  { const int n = data->cols()/2;
    assert (data->cols() == 2*n);
    for (int i=0; i<data->rows(); ++i)
    { data->row(i).head(n) /= sqrt2*data->row(i).head(n).norm();
      data->row(i).tail(n) /= sqrt2*data->row(i).tail(n).norm();
    }
  }
}


int
KMeansClusters::closest_cluster (RowVector const& r, Matrix const& m) const
{
  float minD = mDist(r,m.row(0));
  int   minI = 0;
  for (int i=1; i<m.rows(); ++i)
  { float d = mDist(r,m.row(i));
    if (d < minD) { minD = d; minI = i; }
  }
  return minI;
}

void
KMeansClusters::find_clusters(int maxIterations)
{
  // init new cluster centers to top rows
  Matrix newCenters = Matrix (mNClusters,mData.cols());
  newCenters = mData.topLeftCorner(mNClusters, mData.cols());
  // set 'old' center to zero so initial change is infinite
  Matrix oldCenters = Matrix::Zero(mNClusters, mData.cols());

  std::vector<int> counts (mNClusters);
  int itCount = 0;
  Matrix *pNew = &newCenters;
  Matrix *pOld = &oldCenters;
  const int n = mData.cols()/2;
  while ((0.001 < relative_squared_distance(*pNew,*pOld)) && (++itCount < maxIterations))
  { // assign cases to clusters
    for(int i=0; i<mData.rows(); ++i)
      mDataClusterIndex[i] = closest_cluster(mData.row(i), *pNew);
    std::swap(pNew,pOld);
    // calculate new centers of clusters
    *pNew = Matrix::Zero(mNClusters, mData.cols());
    for(int i=0; i<mNClusters; ++i) counts[i]=0;
    for(int i=0; i<mData.rows(); ++i)
    { (*pNew).row(mDataClusterIndex[i]) += mData.row(i) * mWeights(i);
      counts[mDataClusterIndex[i]] += mWeights(i);
    }
    for(int i=0; i<mNClusters; ++i)
    { pNew->row(i).array() /= counts[i];
      if (!mUseL2)
      { pNew->row(i).head(n) /= sqrt2*pNew->row(i).head(n).norm();
	pNew->row(i).tail(n) /= sqrt2*pNew->row(i).tail(n).norm();
      }
    }
    // std::clog << "Center counts   \n" ;
    // for(int i=0;i<mNClusters; ++i) std::clog << counts[i] << " "; std::clog << std::endl;
    // std::clog << "Centers at step " << itCount << ":\n" << *pNew << std::endl;
  }
  for(int i=0; i<mData.rows(); ++i)
    mDataClusterIndex[i] = closest_cluster(mData.row(i), *pNew);
  mClusterCenters = *pNew;
}


KMeansClusters::ClusterMap
KMeansClusters::cluster_map() const
{
  KMeansClusters::ClusterMap m;
  for(int i=0; i<(int)mDataClusterIndex.size(); ++i)
    m[mDataClusterIndex[i]].push_back(i);
  return m;
}


double
KMeansClusters::relative_squared_distance (Matrix const& newCenters, Matrix const& oldCenters) const
{
  // avg squared distance relative to avg squared size
  double dist = (newCenters.array() - oldCenters.array()).matrix().squaredNorm()/oldCenters.squaredNorm();
  std::clog << "KMCS: " << ((mUseL2) ? "L2 " : "Cos ") << "distance = " << dist << std::endl;
  return dist;
}


void
KMeansClusters::print_to_stream (std::ostream& os) const 
{
  os << "K-Means cluster analysis, with " << mNClusters << " clusters with n=" << mData.rows() << std::endl;
  ClusterMap m (cluster_map());
  const int nShown = 25;
  for(int i=0; i<mNClusters; ++i)
  { int clusterSize = m[i].size();
    os << "Cluster " << i << " (" << clusterSize << " items): ";
    { for(int j=0; j<min(clusterSize, nShown); ++j)  os << " " << m[i][j];
      if(clusterSize > nShown)                       os << " + " << clusterSize-nShown << " more.";
      os << std::endl;
    }
  }
}
