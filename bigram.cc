#include "bigram.h"

/*
  R Code to check SVD
 
  X <- read.table("/Users/bob/C/text/random_projection.txt")
  udv <- svd (X);  udv$d
  for(j in 1:10) {
     cat(j,"  ",any(is.nan(X[,j])),"  \n")
     cat(which (is.nan(X[,j])), "  \n");  }
*/

typedef Eigen::VectorXf                   Vector;
typedef Eigen::VectorXi                   IntVector;
typedef Eigen::MatrixXf                   Matrix;

typedef Eigen::SparseMatrix<float,Eigen::RowMajor> SparseMatrix;


void
print_time(std::ostream& os, std::string const& s, clock_t const& start, clock_t const& stop)
{
  clock_t diff (stop-start);
  os << "TIME: " << s << " [" << ((float)diff)/CLOCKS_PER_SEC << " sec]\n";
}

void
parse_arguments(int argc, char** argv,		float &posThreshold, int &nProj, int &scaling, char &dist, int &weighting, int &nClus, int &nIter, int &nPrint);

  
int main(int argc, char **argv)
{
  // argument defaults
  float posThreshold = 0.0;        // accum all pos tags with relative frequency below this threshold
  int   nProjections = 30;
  int   scaling      =  0;
  char  distance     ='2';         // letter for type of distance to use
  int   weighting    =  1;
  int   nClusters    = 10;
  int   nIterations  = 10;
  int   nPrint       =  0;

  parse_arguments(argc, argv, posThreshold, nProjections, scaling, distance, weighting, nClusters, nIterations, nPrint);
  std::clog << "MAIN: Arguments threshold=" << posThreshold << " nProjections=" << nProjections << " scaling=" << scaling
	    << " distance=" << distance << " weighting=" << weighting << " nClusters=" << nClusters << " nIterations=" << nIterations
	    << " nPrint=" << nPrint << std::endl;
  
  using std::string;
  std::ostringstream ss;
  
  // read parsed tokens and pos from input
  clock_t startTime = clock();
  TokenManager tokenManager (std::cin, posThreshold);
  print_time(std::clog, "Read tokens from cin, sort and assign IDs in TokenManager.", startTime, clock());
  if (nPrint) tokenManager.print_tags(nPrint);
  int nAmbiguous = tokenManager.n_ambiguous();
  std::clog << "MAIN: Source has " << nAmbiguous << " ambiguous tokens among " << tokenManager.input_length()
	    << "  (" << 100.0*(1.0 - ((float)nAmbiguous)/tokenManager.input_length()) << "% pure)" << std::endl;
  
  // build the sparse bigram array using integer id for words
  const int nTokens (tokenManager.n_unique_tokens());
  SparseMatrix B(nTokens,nTokens);
  {
    startTime = clock();
    std::map<std::pair<int,int>,int> bgramMap;
    tokenManager.fill_bigram_map(bgramMap);
    {
      typedef Eigen::Triplet<float> T;
      std::list<T> triplets (nTokens);
      for(auto it = bgramMap.cbegin(); it != bgramMap.cend(); ++it)
	triplets.push_back(T(it->first.first, it->first.second, it->second));
      B.setFromTriplets(triplets.begin(), triplets.end());
    }
    ss << "Init sparse bigram B[" << B.rows() << "x" << B.cols() << "] from map.";
    print_time(std::clog, ss.str(), startTime, clock());
    ss.str("");
  }

  // count number of times each type appears; check for an empty row or col in B
  IntVector typeCounts = IntVector::Zero(B.rows());
  {
    startTime = clock();
    std::vector<int> iZero;
    for(int i=0; i<B.rows(); ++i)
    { typeCounts[i] = B.row(i).sum();
      if (0 == typeCounts[i]) iZero.push_back(i);
      // if (B.col(i).sum() == 0) iZero.push_back(-i);
    }
    ss << "Zero check finds " << iZero.size() << " empty rows and columns.";
    if(!iZero.empty())
    { ss << std::endl;
      for (auto it=iZero.begin(); it!=iZero.end(); ++it) ss << *it << " ";
    }
    print_time(std::clog, ss.str(), startTime, clock());
    ss.str("");
  }
  
  // Random projections of row and column spaces 
  startTime = clock();
  Matrix RP (B.rows(), 2*nProjections);
  RP.leftCols (nProjections) = B             * Matrix::Random(B.cols(), nProjections);
  RP.rightCols(nProjections) = B.transpose() * Matrix::Random(B.rows(), nProjections);
  ss << "Compute random projection RP[" << RP.rows() << "x" << RP.cols() << "]";
  print_time(std::clog, ss.str(), startTime, clock());
  ss.str("");

  // optional scaling of projection rows
  if(scaling)
  { for(int i = 0; i<RP.rows(); ++i)
    { float norm = RP.leftCols(nProjections).row(i).norm();
      RP.leftCols(nProjections).row(i) /= norm;
    }
    for(int i = 0; i<RP.rows(); ++i)
    { float norm = RP.rightCols(nProjections).row(i).norm();
      RP.rightCols(nProjections).row(i) /= norm;
    }
    for (int i = 0; i<2; ++i)
      std::clog << " Norm of random projection row " << i << " is " << RP.row(i).norm() << std::endl;
  }
  
  // cluster tokens using k-means
  startTime = clock();
  std::vector<int> tags;
  {
    IntVector wts;
    if (weighting)
      wts = typeCounts;
    else
      wts = IntVector::Ones(RP.rows());
    k_means::Distance d;
    k_means::Renorm   g;
    if ('2' == distance)
    { d = k_means::l2_distance;      g = k_means::identity;
    } else
    { d = k_means::cosine_distance;  g = k_means::two_balls;
    }
    KMeansClusters clusters (RP, wts, d, g, nClusters, nIterations);
    tags = clusters.cluster_tags();
  }
  ss << "Compute " << nClusters << " cluster centers.";
  print_time(std::clog, ss.str(), startTime, clock());
  ss.str("");

  // Measure classifier error rate (cross-classify, then count number not at max)
  {
    std::map<string,int> posMap = tokenManager.POS_map();
    std::map<string,int> posIndex;
    std::vector<string>  indexPos (posMap.size());
    int index = 0;
    for(auto it=posMap.cbegin(); it != posMap.cend(); ++it)
    { posIndex[it->first] = index;
      indexPos[index]=it->first;
      ++index;
    }
    Eigen::MatrixXi crossTab = Eigen::MatrixXi::Zero(nClusters, posMap.size());
    index = 0;
    for(auto it=tokenManager.token_list_begin(); it != tokenManager.token_list_end(); ++it)
    { crossTab(tags[tokenManager[it->first]], posIndex[it->second]) += 1;
      ++index;
    }
    int width = 2+log10(crossTab.array().maxCoeff());
    for(int col=0; col<(int)indexPos.size(); ++col)
      std::clog << std::setw(width) << indexPos[col];
    std::clog << std::endl << crossTab << std::endl;
    /*
    for(int col=0; col<crossTab.cols(); ++col)
      std::clog << "Column margin in crosstab: " << indexPos[col] << "  " << posMap[indexPos[col]] << " == " << crossTab.col(col).sum() << std::endl;
    */
    int nCorrect=0;
    for(int i=0; i<crossTab.rows(); ++i)
      nCorrect += crossTab.row(i).maxCoeff();
    std::clog << "Classify " << nCorrect << " correctly out of " << tokenManager.input_length()
	      << "(" << 100.0 * ((float) nCorrect)/tokenManager.input_length() << "%)\n";
  }
  
  // Write original tags and cluster ids to file
  {
    std::ios_base::openmode mode = std::ios_base::trunc;
    std::string fileName ("/Users/bob/Desktop/tags.txt");
    std::ofstream file (fileName.c_str(), mode);
    file << "Token\tPOS\tCluster" << std::endl;
    for(auto it = tokenManager.token_list_begin(); it != tokenManager.token_list_end(); ++it)
      file << it->first << "\t" << it->second << "\t" << tags[tokenManager[it->first]] << std::endl;
  }
  
  // SVD of random projection array
  if (false)
  { startTime = clock();
    Eigen::JacobiSVD<Matrix, Eigen::HouseholderQRPreconditioner> svd(RP, Eigen::ComputeThinU);
    Matrix UD = svd.matrixU() * svd.singularValues().asDiagonal();
    // std::clog << " U matrix      : \n" << UD.topLeftCorner(10,nProjections) << std::endl;
    std::clog << "Singular values: \n" << svd.singularValues().transpose().head(nProjections) << std::endl;
    print_time(std::clog, "Compute SVD and extract U*D.", startTime, clock());
  }
  
  // print first items in the map
  /*
    for(auto it = sortedTokenMap.begin(); it != sortedTokenMap.end(); ++it)
    std::cout << it->second << " " << -it->first << std::endl;  
  */
  return 0;
}



void
parse_arguments(int argc, char** argv,		float &threshold, int &nProj, int &scaling, char &distance, int &weighting, int &nClus, int &nIter, int &nPrint)
{
  static struct option long_options[] = {
    // bigram prep options
    {"threshold",    required_argument, 0, 't'},
    {"projections",  required_argument, 0, 'r'},
    {"scaling",      required_argument, 0, 's'},
    // clustering options
    {"distance",     required_argument, 0, 'd'},
    {"weighting",    required_argument, 0, 'w'},
    {"clusters",     required_argument, 0, 'c'},
    {"iterations",   required_argument, 0, 'n'},
    // misc
    {"print",        required_argument, 0, 'p'},
    {0, 0, 0, 0}                             // terminator 
  };
  int key;
  int option_index = 0;
  while (-1 !=(key = getopt_long (argc, argv, "t:r:s:w:c:n:p:", long_options, &option_index))) // colon means has argument
  {
    // std::cout << "Option key " << char(key) << " for option " << long_options[option_index].name << ", option_index=" << option_index << std::endl;
    switch (key)
    {
    case 't' :
      {
	threshold = read_utils::lexical_cast<float>(optarg);
	break;
      }
    case 'n' :
      {
	nIter = read_utils::lexical_cast<int>(optarg);
	break;
      }
    case 'r' : 
      {
	nProj = read_utils::lexical_cast<int>(optarg);
	break;
      }
    case 's' : 
      {
	scaling = read_utils::lexical_cast<int>(optarg);
	break;
      }
    case 'd' : 
      {
	distance = read_utils::lexical_cast<char>(optarg);
	break;
      }
    case 'w' : 
      {
	weighting = read_utils::lexical_cast<int>(optarg);
	break;
      }
    case 'c' : 
      {
	nClus= read_utils::lexical_cast<int>(optarg);
	break;
      }
    case 'p' : 
      {
	nPrint = read_utils::lexical_cast<int>(optarg);
	break;
      }
    default:
      {
	std::cout << "PARSE: Option not recognized; returning.\n";
      }
    } // switch
  } // while
}
