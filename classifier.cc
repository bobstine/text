#include "classifier.h"

#include "iterators.h"

//     ClusterClassifier     ClusterClassifier     ClusterClassifier     ClusterClassifier     ClusterClassifier     

 void
 ClusterClassifier::fill_map (KMeansClusters const& clusters, TokenManager const& tm)
 {
   KMeansClusters::ClusterMap cMap = clusters.cluster_map();

   POSVector     clusterToPOSVec  = POS_vector_from_clusters(cMap, tm);
   TypeToIntMap  typeToClusterMap = type_to_cluster_map(cMap, tm);
   for(auto it = typeToClusterMap.cbegin(); it != typeToClusterMap.cend(); ++it)
   { Type t = it->first;
     int  i = it->second;
     mMap[t]=clusterToPOSVec[i];
   }
 }


 std::vector<POS>
 ClusterClassifier::POS_vector_from_clusters (KMeansClusters::ClusterMap const& clusterMap, TokenManager const& tm) const
 {
   std::vector<Type>     typeVec = tm.type_vector();
   std::vector<POS>      posVec;
   for(auto it=clusterMap.cbegin(); it!=clusterMap.cend(); ++it)
   { std::map<POS,int> clusterPOSCount;
     for(size_t i=0; i<it->second.size(); ++i)
       ++clusterPOSCount[ tm.POS_of_type(typeVec[it->second[i]]) ];
     string dummy("dummy");
     POS maxPOS(dummy);
     int maxCount = 0;
     for(auto iter=clusterPOSCount.cbegin(); iter!=clusterPOSCount.cend(); ++iter)
       if(iter->second > maxCount)
       { maxCount = iter->second;
	 maxPOS   = iter->first;
       }
     posVec.push_back(maxPOS);
   }
   return posVec;
 }


 ClusterClassifier::TypeToIntMap
 ClusterClassifier::type_to_cluster_map (KMeansClusters::ClusterMap const& clusterMap, TokenManager const& tm) const
 {
   std::vector<Type> typeVec = tm.type_vector();
   std::map<Type,int> mapTypeToCluster;
   for(auto it=clusterMap.cbegin(); it!=clusterMap.cend(); ++it)
     for(size_t i=0; i<it->second.size(); ++i)
       mapTypeToCluster[ typeVec[it->second[i]] ] = it->first;
   return mapTypeToCluster;
 }


//     ConfusionMatrix     ConfusionMatrix     ConfusionMatrix     ConfusionMatrix     ConfusionMatrix     ConfusionMatrix

//   Iterators

class ExtractPOSString: public std::unary_function<std::pair<Type,POS>, string>
{
public:
  string operator()(std::pair<Type,POS> const& p) const { return p.second; }
};


//    Builder

ConfusionMatrix
make_confusion_matrix (ClusterClassifier const& classifier, TokenManager const& tm)
{
  ConfusionMatrix cm(make_function_iterator(tm.token_list_begin(), ExtractPOSString()),
		     make_function_iterator(tm.token_list_end  (), ExtractPOSString()),
		     make_function_iterator(tm.token_list_begin(), classifier)
		     );
  return cm;
}


