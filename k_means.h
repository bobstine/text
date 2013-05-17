#ifndef _K_MEANS_H_
#define _K_MEANS_H_


#include <Eigen/Core>

#include <vector>
#include <map>

typedef std::map<int, std::vector<int>>  ClusterMap;

ClusterMap k_means_cluster_map(Eigen::MatrixXf const& data, int nClusters);

#endif

