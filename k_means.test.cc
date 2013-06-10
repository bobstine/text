#include "k_means.h"
#include "eigen_utils.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>

const int nRows = 30;
const int nCols = 12;

using std::cout;
using std::endl;
using std::string;

int main (void)
{
  // data vary from [-10,10] after randomize seed
  srand(time(NULL));
  Eigen::MatrixXf baseData = 10.0 * Eigen::MatrixXf::Random(nRows,nCols+1);

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
  
  const int nClusters = 4;
  const int maxIter = 10;
  
  if (false)   // l2 cluster
  {
    Eigen::MatrixXf data = baseData;
    bool useL2 = true;
    bool scale = true;
    
    for(int i=0; i<nRows; ++i)
      data.row(i).array() += (i % nClusters);
    Eigen::VectorXi wts = Eigen::VectorXi::Ones(nRows);
    KMeansClusters clusters(data.leftCols(nCols),wts,caseLabels,useL2,scale,nClusters,maxIter);
    clusters.print_to_stream(std::cout);
    std::vector<string> labels;
    int ct=0;
    for(auto it = clusters.item_cluster_tag_begin(); it!= clusters.item_cluster_tag_end(); ++it)
    { labels.push_back(*it); 
      std::cout << "TEST: label for item " << ct++ << " is " << *it << endl;
    }
  }

  if (true)        // cosine centers
  {
    Eigen::MatrixXf data = baseData;
    bool useL2 = false;
    bool scale = false;
    
    float shift = 100;
    for(int i=0; i<nRows; ++i)
    { int offset = i % nClusters;
      for (int j=offset; j<nCols; j += nClusters)
        data(i,j) += shift;
    }
    //    std::cout << data << std::endl;
    Eigen::VectorXi wts = Eigen::VectorXi::Ones(nRows);
    KMeansClusters clusters(data.leftCols(nCols), wts, caseLabels, useL2, scale, nClusters, maxIter);
    std::cout << "TEST: Cluster label purity is " << clusters.purity() << std::endl;
    clusters.print_to_stream(std::cout, true);
    // check two ways to get group labels
    std::vector<string> labelA (nRows);
    clusters.fill_with_fitted_cluster_tags(labelA.begin(), labelA.end());
    std::vector<string> labelB;
    for(auto it = clusters.item_cluster_tag_begin(); it!= clusters.item_cluster_tag_end(); ++it)
      labelB.push_back(*it);
    std::cout << "TEST:   Showing data in first column, then two copies of fitted labels\nData   A   B\n";
    for (int i=0; i<nRows; ++i)
      std::cout << caseLabels[i] << "   " << labelA[i] << " " << labelB[i] << std::endl;  // labelA should match labelB
  }
}
