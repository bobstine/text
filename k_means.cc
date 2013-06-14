#include "k_means.h"
#include "order.h"

#include <utility>
#include <iostream>
#include <fstream>  // debug temporary
#include <set>      // debug temporary

const float sqrt2 = 1.41421356;

static std::string messageTag = "KMCS: ";


float
check_norm(float norm, std::string msg)
{
  if ((norm <= 0) || (!std::isfinite(norm)))
  { std::cerr << messageTag << "*** ERROR *** Norm=" << norm << " " << msg << std::endl;
    return 1.0;
  }
  else
    return norm;
}

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
KMeansClusters::prepare_data(Matrix *mat) const
{
  if (mScaleData)
  { if (mBidirectional)
    { assert(0 == mat->cols() % 2);
      const int n = mat->cols()/2;
      for(int i=0; i<mat->rows(); ++i)
      { mat->row(i).head(n).normalize();
	mat->row(i).tail(n).normalize();
	mat->row(i).array() /= sqrt2;
      }
    }
    else
      for (int i=0; i<mat->rows(); ++i)
	mat->row(i).normalize();
  }
}
  
int
KMeansClusters::closest_cluster (RowVector const& r, Matrix const& m) const
{
  float minD = distance(r, m.row(0));
  int   minI = 0;
  for (int i=1; i<m.rows(); ++i)
  { float d = distance(r,m.row(i));
    if (d < minD) { minD = d; minI = i; }
  }
  return minI;
}

void
KMeansClusters::find_clusters(int maxIterations)
{

  // write data to file
  std::ofstream kos ("/Users/bob/Desktop/kmeans_data.txt");
  kos << mData << std::endl;

  // init cluster centers to top rows
  Matrix newCenters = Matrix (mNClusters,mData.cols());
  newCenters = mData.topLeftCorner(mNClusters, mData.cols());
  // set 'old' center to zero so initial rel change is infinite
  Matrix oldCenters = Matrix::Zero(mNClusters, mData.cols());

  std::vector<int> counts (mNClusters);
  int itCount = 0;
  Matrix *pNew = &newCenters;
  Matrix *pOld = &oldCenters;
  // const int n = mData.cols()/2;
  while ((0.001 < relative_squared_distance(*pNew,*pOld))  // prints message
	 && (++itCount < maxIterations))
  {
    std::set<int> clusterZero;
    for(int i=0; i<mData.rows(); ++i)
    { mDataClusterIndex[i] = closest_cluster(mData.row(i), *pNew);
      if (0 == mDataClusterIndex[i])
	clusterZero.insert(i);
    }
    // avg distance to cluster 0 centroid
    std::clog << messageTag << clusterZero.size() << " cases in Cluster 0.";
    if(clusterZero.size()>0)
    { float mean = 0.0;
      for(size_t i=0; i<clusterZero.size(); ++i)
	mean += distance(mData.row(i), pNew->row(0));
      mean /= clusterZero.size();
      std::clog << messageTag << "Average distance to center is " << mean;
    }
    else std::clog << std::endl;
     
    std::swap(pNew,pOld);
    // calculate new centers of clusters
    *pNew = Matrix::Zero(mNClusters, mData.cols());
    for(int i=0; i<mNClusters; ++i) counts[i]=0;
    for(int i=0; i<mData.rows(); ++i)
    { pNew->row(mDataClusterIndex[i]) += mData.row(i) * mWeights(i);
      counts[mDataClusterIndex[i]] += mWeights(i);
    }
    for(int c=0; c<mNClusters; ++c)
      if(counts[c])
      {	pNew->row(c).array() /= counts[c];
	if (mScaleCentroid)
	{ if (mBidirectional)
	  { int half = mData.cols()/2;
	    pNew->row(c).head(half).normalize();
	    pNew->row(c).tail(half).normalize();
	  }
	  else pNew->row(c).normalize();
	}
      }
      else std::clog << messageTag << "Count=0 in cluster " << c << std::endl;
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


std::vector<std::pair<int,float> >
KMeansClusters::average_centroid_dist(Matrix const& centroids) const
{
  std::vector<std::pair<int,float>> result;
  KMeansClusters::ClusterMap m = cluster_map();
  for(auto it=m.cbegin(); it != m.cend(); ++it)
  { int                     cluster = it->first;
    std::vector<int> const& cases   = it->second;
    float xbar=0;
    for (size_t i=0; i<cases.size(); ++i)
      xbar += distance(mData.row(cases[i]), centroids.row(cluster));
    if(cases.size()) xbar /=  cases.size();
    result.push_back( std::make_pair(cases.size(), xbar) );
  }
  return result;
}


double
KMeansClusters::relative_squared_distance (Matrix const& newCenters, Matrix const& oldCenters) const
{
  // avg squared distance relative to avg squared size
  double oldNorm  = oldCenters.squaredNorm();
  double diffNorm = (newCenters.array() - oldCenters.array()).matrix().squaredNorm();
  double dist = diffNorm/oldNorm;
  std::clog << messageTag << "distance = " << dist << std::endl;
  return dist;
}


void
KMeansClusters::print_to_stream (std::ostream& os) const 
{
  os << "K-Means cluster analysis, building " << mNClusters << " clusters from n=" << mData.rows() << " cases." << std::endl;
  ClusterMap m (cluster_map());
  const int maxClustersShown = 15;
  const int maxShownEach     = 15;
  for(int i=0; i<min(mNClusters, maxClustersShown); ++i)
  { int clusterSize = m[i].size();
    os << "Cluster " << i << " (" << clusterSize << " items): ";
    { for(int j=0; j<min(clusterSize, maxShownEach); ++j)  os << " " << m[i][j];
      if(clusterSize > maxShownEach)                       os << " + " << clusterSize-maxShownEach << " more.";
      os << std::endl;
    }
  }
}
