#include "k_means.h"
#include "order.h"

#include <utility>
#include <iostream>

const float sqrt2 = 1.41421356;

std::vector<int>
KMeansClusters::assign_cluster_indices (Matrix *data) const
{
  prepare_data(data);
  std::vector<int> intTags;
  for (int i=0; i<data->rows(); ++i)
    intTags.push_back(closest_cluster(data->row(i), mClusterCenters));
  return intTags;
}

std::vector<string>
KMeansClusters::assign_cluster_labels (Matrix *data) const
{
  std::vector<int> indices (assign_cluster_indices(data));
  std::vector<string> labels (indices.size());
  for(size_t i=0; i<indices.size(); ++i)
    labels[i] = mClusterLabels[indices[i]];
  return labels;
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


KMeansClusters::Map
KMeansClusters::cluster_map() const
{
  KMeansClusters::Map m;
  for(int i=0; i<(int)mDataClusterIndex.size(); ++i)
    m[mDataClusterIndex[i]].push_back(i);
  return m;
}


void
KMeansClusters::label_clusters (std::vector<std::string> const& caseLabels)
{
  // lookup table for case labels
  std::map<std::string, int> labelMap;
  for(int i=0; i<(int)caseLabels.size(); ++i)
  { if (0 == labelMap.count(caseLabels[i]))
    { labelMap[caseLabels[i]]=(int)mUniqueLabels.size();
      mUniqueLabels.push_back(caseLabels[i]);
    }
    mDataLabelIndex[i] = labelMap[caseLabels[i]];
  }
  // identify indices of cases in various clusters
  KMeansClusters::Map m = cluster_map();
  for(int c=0; c<mNClusters; ++c)
  { std::vector<int> const& casesInCluster = m[c];
    std::map<std::string,int> counts;
    for (int i=0; i<(int)casesInCluster.size(); ++i)
      ++counts[caseLabels[casesInCluster[i]]];
    int max=0;
    std::string maxLabel = "";
    for (auto it=counts.cbegin(); it!=counts.cend(); ++it)
      if (it->second > max)
      { max = it->second;
	maxLabel = it->first;
      }
    mClusterLabels[c]=maxLabel;
  }
}
  

double
KMeansClusters::relative_squared_distance (Matrix const& newCenters, Matrix const& oldCenters) const
{
  // avg squared distance relative to avg squared size
  double dist = (newCenters.array() - oldCenters.array()).matrix().squaredNorm()/oldCenters.squaredNorm();
  std::clog << "KMCS: " << ((mUseL2) ? "L2 " : "Cos ") << "distance = " << dist << std::endl;
  return dist;
}
  
float
KMeansClusters::purity() const
{
  int nRight=0;
  for(int i=0; i<mData.rows(); ++i)
    if(mClusterLabels[mDataClusterIndex[i]] == mUniqueLabels[mDataLabelIndex[i]])
      ++nRight;
  return ((float)nRight)/mData.rows();
}

void
KMeansClusters::fill_with_fitted_cluster_labels(KMeansClusters::OutIter b, KMeansClusters::OutIter e)   const
{
  for(int i=0; i<mData.rows(); ++i)
  { assert (b != e);
    *b = mClusterLabels[mDataClusterIndex[i]];
    ++b;
  }
}

std::vector<string>
KMeansClusters::fitted_cluster_labels()   const
{
  StringVector labels(mData.rows());
  for(int i=0; i<mData.rows(); ++i)
    labels[i] = mClusterLabels[mDataClusterIndex[i]];
  return labels;
}


void
KMeansClusters::print_to_stream (std::ostream& os, bool showTag) const 
{
  os << "K-Means cluster analysis, with " << mNClusters << " clusters of purity " << purity()
     << " from observing " << mData.rows() << " items with the following "
     << mUniqueLabels.size() << " unique group labels:\n   -->  ";
  for(int i=0; i<(int)mUniqueLabels.size(); ++i)
    os << mUniqueLabels[i] << " ";
  os << std::endl;
  Map m (cluster_map());
  for(int i=0; i<mNClusters; ++i)
  { int clusterSize = m[i].size();
    os << "Cluster " << i << " <" << mClusterLabels[i] << "," << clusterSize << " items>: ";
    if (clusterSize < 15)  // show the whole thing
    { for(int j=0; j<clusterSize; ++j)
      { int caseIndex = m[i][j];
	if (showTag)
	  os << " (" << caseIndex << " " << mUniqueLabels[mDataLabelIndex[caseIndex]] << ")";
	else
	  os << " " << caseIndex;
      }
      os << std::endl;
    }
    else // show summary for cluster
    { std::vector<int> labelCounts (mUniqueLabels.size());
      for(int j=0; j<clusterSize; ++j)
      { int caseIndex = m[i][j];
	++labelCounts.at(mDataLabelIndex[caseIndex]);
      }
      const bool descending = true;
      const size_t maxPrint = 10;
      std::vector<size_t> sortIdx = order(labelCounts, descending);
      size_t nPrint=0;
      for(size_t j=0; j<labelCounts.size(); ++j)
      { size_t ii = sortIdx[j];
	if (labelCounts[ii]>0)
	{ ++nPrint;
	  os << " [" << mUniqueLabels.at(ii) << "," << labelCounts.at(ii) << "]";
	}
	if (maxPrint == nPrint)
	{ os << " ... ";
	  break;
	}
      }
      os << std::endl;
    }
  }
}
