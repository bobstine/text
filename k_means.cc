#include "k_means.h"
#include <iostream>



KMeansClusters::Map
KMeansClusters::cluster_map() const
{
  KMeansClusters::Map m;
  for(int i=0; i<(int)mClusterTags.size(); ++i)
    m[mClusterTags[i]].push_back(i);
  return m;
}

std::vector<int>
KMeansClusters::cluster_tags () const
{
  return mClusterTags;
}


int
KMeansClusters::closest_cluster (int row) const
{
  float minD = (mData.row(row) - mClusterCenters.row(0)).squaredNorm();
  int   minI = 0;
  for (int i=1; i<mClusterCenters.rows(); ++i)
  { float d = (mData.row(row) - mClusterCenters.row(i)).squaredNorm();
    if (d < minD) { minD = d; minI = i; }
  }
  return minI;
}

double
KMeansClusters::relative_squared_distance (Matrix const& newCenters) const
{
  // avg squared distance relative to avg squared size
  double dist = (mClusterCenters.array() - newCenters.array()).matrix().squaredNorm()/mClusterCenters.squaredNorm();
  std::clog << "KMCS: disance = " << dist << std::endl;
  return dist;
}

void
KMeansClusters::find_clusters(int maxIterations)
{
  // init new cluster centers to top rows
  Matrix centersNew = mData.topLeftCorner(mNClusters, mData.cols());
  mClusterCenters   = Matrix::Zero       (mNClusters, mData.cols());

  std::vector<int> counts   (mNClusters);
  int itCount = 0;
  while ((relative_squared_distance(centersNew) > 0.01) && (++itCount < maxIterations))
  { mClusterCenters = centersNew;
    // assign cases to clusters
    for(int i=0; i<mData.rows(); ++i)
      mClusterTags[i] = closest_cluster(i);
    // calculate new centers of clusters
    centersNew = Matrix::Zero(centersNew.rows(), centersNew.cols());
    for(int i=0; i<mData.rows(); ++i)
    { centersNew.row(mClusterTags[i]) += mData.row(i);
      ++counts[mClusterTags[i]];
    }
    for(int i=0; i<centersNew.rows(); ++i)
      centersNew.row(i).array() /= counts[i];
  }
}
  

void
KMeansClusters::print_to_stream (std::ostream& os) const
{
  Map m (cluster_map());
  for(int i=0; i<mNClusters; ++i)
  { os << "Cluster " << i << ": ";
    for(int j=0; j<(int)m[i].size(); ++j)
      os << " " << m[i][j];
    os << std::endl;
  }
}
  
