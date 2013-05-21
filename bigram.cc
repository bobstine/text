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
typedef Eigen::MatrixXf                   Matrix;

typedef Eigen::SparseMatrix<float,Eigen::RowMajor> SparseMatrix;


void
print_time(std::ostream& os, std::string const& s, clock_t const& start, clock_t const& stop)
{
  clock_t diff (stop-start);
  os << "TIME: " << s << " [" << ((float)diff)/CLOCKS_PER_SEC << " sec]\n";
}

void
parse_arguments(int argc, char** argv,		float &posThreshold, int &nProj, int &scaling, int &nClus, int &nIter, int &nPrint);

  
int main(int argc, char **argv)
{
  // argument defaults
  float posThreshold = 0.0;        // accum all pos tags with relative frequency below this threshold
  int   nProjections = 30;
  int   scaling      =  0;
  int   nClusters    = 10;
  int   nIterations  = 10;
  int   nPrint       =  0;

  parse_arguments(argc, argv, posThreshold, nProjections, scaling, nClusters, nIterations, nPrint);
  std::clog << "MAIN: Echo arguments, threshold=" << posThreshold << " nProjections=" << nProjections << " scaling=" << scaling
	    << " nClusters=" << nClusters << " nIterations=" << nIterations << " nPrint=" << nPrint << std::endl;
  
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

  // test by summing elements in B
  if(false)
  { ss << "Sums of first rows of B are ";
    startTime = clock();
    {
      Vector one  (Vector::Constant(nTokens,1.0));
      Vector sums (Vector::Zero    (nTokens));
      sums = B * one;
      for(int i=0; i<5; ++i)
	ss << "[" << i << "]=" << sums[i] << "   ";
      ss << " with +/+/B = " << sums.sum() << ".";
    }
    print_time(std::clog, ss.str(), startTime, clock());
    ss.str("");
  }
  
  // Random projection of the array
  startTime = clock();
  Matrix RP = B * Matrix::Random(B.cols(),nProjections);;
  if (scaling)
  { Vector recipNorms(B.rows());
    for (int i=0; i<B.rows(); ++i)
    { float ss = B.row(i).norm();
      if (0 == ss)
      { ss = 1;
	std::clog << "MAIN: Zero norm for B[" << i << "] with token ----->" << tokenManager[i]
		  << "<----- with count " << tokenManager.token_freq(tokenManager[i]) << std::endl;
      }
      else
	recipNorms(i) = 1.0/ss;
    }
    RP = recipNorms.asDiagonal() * RP;
  }
  ss << "Compute random projection RP[" << RP.rows() << "x" << RP.cols() << "].";
  print_time(std::clog, ss.str(), startTime, clock());
  ss.str("");
  /*  std::clog << " First 10 cols of first and last 5 rows of the RP matrix : \n"
	    << RP.topLeftCorner(5,10) << "\n  ...\n"
	    << RP.bottomLeftCorner(5,10) << std::endl;
      write_matrix_to_file("/Users/bob/Desktop/rand_projection.txt",RP);
  */
  // cluster tokens using k-means
  KMeansClusters clusters (RP, nClusters, nIterations);
  std::vector<int> tags = clusters.cluster_tags();

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
    /*
    for(int col=0; col<crossTab.cols(); ++col)
      std::clog << "Column margin in crosstab: " << indexPos[col] << "  " << posMap[indexPos[col]] << " == " << crossTab.col(col).sum() << std::endl;
    */
    for(int col=0; col<(int)indexPos.size(); ++col)
      std::clog << std::setw(6) << indexPos[col];
    std::clog << std::endl << crossTab << std::endl;
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
parse_arguments(int argc, char** argv,		float &threshold, int &nProj, int &scaling, int &nClus, int &nIter, int &nPrint)
{
  static struct option long_options[] = {
    {"threshold",    required_argument, 0, 't'},
    {"projections",  required_argument, 0, 'd'},
    {"scaling",      required_argument, 0, 's'},
    {"clusters",     required_argument, 0, 'c'},
    {"iterations",   required_argument, 0, 'n'},
    {"print",        required_argument, 0, 'p'},
    {0, 0, 0, 0}                             // terminator 
  };
  int key;
  int option_index = 0;
  while (-1 !=(key = getopt_long (argc, argv, "t:d:s:c:n:p:", long_options, &option_index))) // colon means has argument
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
    case 'd' : 
      {
	nProj = read_utils::lexical_cast<int>(optarg);
	break;
      }
    case 's' : 
      {
	scaling = read_utils::lexical_cast<int>(optarg);
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
