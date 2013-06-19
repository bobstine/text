#include "k_means.h"
#include "order.h"
#include "utils.h"

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
  while ((0.001 < relative_squared_distance(*pNew,*pOld))  // prints message without eol
	 && (++itCount < maxIterations))
  { for(int i=0; i<mData.rows(); ++i)
      mDataClusterIndex[i] = closest_cluster(mData.row(i), *pNew);
    // show concentrations within clusters
    std::clog << " with sqrt(avg cluster var)=" << sqrt(average_within_cluster_variance(*pNew)) << std::endl;
    // calculate new centers of clusters
    std::swap(pNew,pOld);
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
  std::clog << std::endl;
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


void
KMeansClusters::summarize_rare_cases(std::ostream &os) const
{
  std::map<int,int> w1ClusterCounts; int nW1=0;
  std::map<int,int> w2ClusterCounts; int nW2=0;
  for(int i=0; i<mData.rows(); ++i)
  { int cluster = mDataClusterIndex[i];
    if      (mWeights(i) == 1) { ++w1ClusterCounts[cluster]; ++nW1; }
    else if (mWeights(i) == 2) { ++w2ClusterCounts[cluster]; ++nW2; }
  }
  os << messageTag << nW1 << " cases with weight 1 appear in " << w1ClusterCounts.size() << " clusters:\n       ";
  for(auto p : utils::invert_map_descending(w1ClusterCounts))
    os << "(C" << p.second << ":" << p.first <<") ";
  os << std::endl
     << messageTag << nW2 << " cases with weight 2 appear in " << w2ClusterCounts.size() << " clusters:\n       ";
  for(auto p : utils::invert_map_descending(w2ClusterCounts))
    os << "(C" << p.second << ":" << p.first <<") ";
  os << std::endl;
}


float 
KMeansClusters::average_within_cluster_variance(Matrix const& centroids) const
{
  IntVector counts = IntVector::Zero(mNClusters);
  Vector    ss     = Vector::Zero(mNClusters);
  for(int i=0; i<mData.rows(); ++i)
  { counts(mDataClusterIndex[i]) += mWeights(i);
    ss(mDataClusterIndex[i]) += mWeights(i) * distance(mData.row(i), centroids.row(mDataClusterIndex[i]));
  }
  float avg=0;
  int n = 0;
  for(int c=0; c<mNClusters; ++c)
    if(counts[c]>0)
    { avg += ss[c]/counts[c];
      ++n;
    }      
  return avg/n;
}


double
KMeansClusters::relative_squared_distance (Matrix const& newCenters, Matrix const& oldCenters) const
{
  // avg squared distance relative to avg squared size
  double oldNorm  = oldCenters.squaredNorm();
  double diffNorm = (newCenters.array() - oldCenters.array()).matrix().squaredNorm();
  double dist = diffNorm/oldNorm;
  std::clog << messageTag << "distance = " << dist;
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
