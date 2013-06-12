#ifndef _CLASSIFIER_H_
#define _CLASSIFIER_H_

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




  
  // label the cases
  std::vector<string> caseLabels (nRows);
  { std::vector<string> labels(4);
    labels[0] = "AAA";
    labels[1] = "BBB";
    labels[2] = "CCC";
    labels[3] = "DDD";
    srand( time(NULL) );
    for (int i=0; i<nRows; ++i)
      caseLabels[i] = labels[ rand() % 4 ];  // rand() % 4   or i % 4
  }


    // check two ways to get group labels
    std::vector<string> labelA (nRows);
    clusters.fill_with_fitted_cluster_labels(labelA.begin(), labelA.end());
    std::vector<string> labelB;
    int matchCount = 0; 
    {
      int i=0;
      for(auto it = clusters.data_cluster_label_begin(); it!= clusters.data_cluster_label_end(); ++it)
      { labelB.push_back(*it);
	if (labelB[i] == caseLabels[i]) ++matchCount;
	++i;
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

// percent of items for which input label matches cluster label 
  float          purity()                                                 const;

#endif
