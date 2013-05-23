#include "k_means.h"

#include <utility>
#include <iostream>


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
      (*pNew).row(i) = mRenorm( (*pNew).row(i) );
      if (i < 5)
	std::clog << "Center norm " << i << " = {" << (*pNew).row(i).head(5) << " ... }  = " << (*pNew).row(i).norm() << std::endl;
    }
    // std::clog << "Center counts   \n" ;
    // for(int i=0;i<mNClusters; ++i) std::clog << counts[i] << " "; std::clog << std::endl;
    // std::clog << "Centers at step " << itCount << ":\n" << *pNew << std::endl;
  }
  mClusterCenters = *pNew;
}


KMeansClusters::Map
KMeansClusters::cluster_map() const
{
  KMeansClusters::Map m;
  for(int i=0; i<(int)mClusterTags.size(); ++i)
    m[mClusterTags[i]].push_back(i);
  return m;
}


double
KMeansClusters::relative_squared_distance (Matrix const& newCenters, Matrix const& oldCenters) const
{
  // avg squared distance relative to avg squared size
  double dist = (newCenters.array() - oldCenters.array()).matrix().squaredNorm()/oldCenters.squaredNorm();
  std::clog << "KMCS: distance = " << dist << std::endl;
  return dist;
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
