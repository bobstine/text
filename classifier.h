#ifndef _CLASSIFIER_H_
#define _CLASSIFIER_H_

#include "token_manager.h"
#include "k_means.h"
#include "confusion_matrix.h"

#include <functional>
#include <map>

using std::string;


class ClusterClassifier: public std::unary_function<Type,POS>
{
  typedef std::vector<POS>   POSVector;
  typedef std::map<Type,int> TypeToIntMap;
  typedef std::map<Type,POS> TypeToPOSMap;
  
  TypeToPOSMap mMap;
  
 public:
  ClusterClassifier (KMeansClusters const& km, TokenManager const& tm) : mMap() { fill_map(km, tm); }
  
  POS operator()(Type const& t)                 const;
  POS operator()(std::pair<Type,POS> const& p)  const { return operator()(p.first); }

 private:
  void         fill_map                 (KMeansClusters const& clusters,               TokenManager const& tm);

  POSVector    POS_vector_from_clusters (KMeansClusters::ClusterMap const& clusterMap, TokenManager const& tm) const;
  TypeToIntMap type_to_cluster_map      (KMeansClusters::ClusterMap const& clusterMap, TokenManager const& tm) const; 

};

ConfusionMatrix
make_confusion_matrix (ClusterClassifier const& classifier, TokenManager const& tm);


#endif
