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
print_time(std::string const& s, clock_t const& start, clock_t const& stop)
{
  clock_t diff (stop-start);
  std::cout << "TIME: " << s << " [" << ((float)diff)/CLOCKS_PER_SEC << " sec]\n";
  std::clog << "TIME: " << s << " [" << ((float)diff)/CLOCKS_PER_SEC << " sec]\n";
}

void
parse_arguments(int argc, char** argv,
		int &nSkip, float &posThreshold,
		int &nProj, int &scaling, char &dist, int &weighting, int &nClus, int &nIter, int &nPrint,
		string &validationFileName);

  
int main(int argc, char **argv)
{
  using std::string;
  
  // argument defaults
  int     nSkip        =  0;         // skipped words in bigram
  float   posThreshold = 0.0;        // accum all pos tags with relative frequency below this threshold
  int     nProjections = 30;
  int     scaling      =  0;
  char    distance     ='2';         // identify distance to use  (2 for L2, c for cosine)
  int     weighting    =  1;
  int     nClusters    = 10;
  int     nIterations  = 10;
  int     nPrint       =  0;
  string  vFileName    = "";

  parse_arguments(argc, argv, nSkip, posThreshold, nProjections, scaling, distance, weighting, nClusters, nIterations, nPrint, vFileName);
  std::cout << "MAIN: Arguments skip=" << nSkip << " threshold=" << posThreshold << " nProjections=" << nProjections << " scaling=" << scaling
	    << " distance=" << distance << " weighting=" << weighting << " nClusters=" << nClusters << " nIterations=" << nIterations
	    << " nPrint=" << nPrint << " validation=" << vFileName << std::endl;
  
  using std::string;
  std::ostringstream ss;
  
  // read tagged tokens and pos from input
  clock_t startTime = clock();
  TokenManager tokenManager (std::cin, posThreshold);
  {
    print_time("Read tokens from cin, sort and assign IDs in TokenManager.", startTime, clock());
    if (nPrint) tokenManager.print_tags(nPrint);
    int nAmbiguous = tokenManager.n_ambiguous();
    ss << "MAIN: Source has " << nAmbiguous << " ambiguous tokens among " << tokenManager.input_length()
       << "  (" << 100.0*(1.0 - ((float)nAmbiguous)/tokenManager.input_length()) << "% pure)" << std::endl;
    std::cout << ss.str();
    std::clog << ss.str();
    ss.str("");
  }
  
  // build the sparse bigram array using integer id for words
  const int nTypes (tokenManager.n_types());
  SparseMatrix B(nTypes,nTypes);
  {
    startTime = clock();
    std::map<std::pair<int,int>,int> bgramMap;
    tokenManager.fill_bigram_map(bgramMap, nSkip);
    {
      typedef Eigen::Triplet<float> T;
      std::list<T> triplets (nTypes);
      for(auto it = bgramMap.cbegin(); it != bgramMap.cend(); ++it)
	triplets.push_back(T(it->first.first, it->first.second, it->second));
      B.setFromTriplets(triplets.begin(), triplets.end());
    }
    ss << "Init sparse bigram B[" << B.rows() << "x" << B.cols() << "] from map.";
    print_time(ss.str(), startTime, clock());
    ss.str("");
  }

  // optional validation which is huge and very sparse
  SparseMatrix V(0,0);
  if (0 < vFileName.size())
  { startTime = clock();
    TokenManager validationTM(vFileName, posThreshold);
    int nOOV = validationTM.n_types_oov(tokenManager);
    std::map<std::pair<int,int>,int> bgMap;
    int dim = validationTM.fill_bigram_map(bgMap, nSkip, tokenManager);
    std::clog << "MAIN: Validation dim = " << dim << " with oov count = " << nOOV << "  Resize to " << nTypes + nOOV << std::endl;
    typedef Eigen::Triplet<float> T;
    std::list<T> triplets (validationTM.n_types());
    for(auto it = bgMap.cbegin(); it != bgMap.cend(); ++it)
      triplets.push_back(T(it->first.first, it->first.second, it->second));
    V.resize(nTypes + nOOV, nTypes + nOOV);
    V.setFromTriplets(triplets.begin(), triplets.end());
    ss << "Init validation sparse bigram V[" << V.rows() << "x" << V.cols() << "] from map." << std::endl;
    print_time(startTime, clock());
    ss.str("");
    }
    
  // count number of times each type appears; check for an empty row or col in B; print margins
  IntVector typeCounts = IntVector::Zero(B.rows());
  {
    std::vector<int> iZero;
    for(int i=0; i<B.rows(); ++i)
    { typeCounts[i] = B.row(i).sum();
      if (0 == typeCounts[i]) iZero.push_back(i);
    }
    if(!iZero.empty())
    { std::clog << "Zero check finds " << iZero.size() << " empty rows and columns." << std::endl;
      for (auto it=iZero.begin(); it!=iZero.end(); ++it) std::clog << *it << " ";
    }
  }

  // write the type counts and frequencies to file
  {
    std::map<string, int> posMap (tokenManager.POS_map());
    IntVector posCounts (posMap.size());
    int i = 0;
    for(auto it=posMap.cbegin(); it != posMap.cend(); ++it)
    { posCounts[i++] = it->second; }
    std::ofstream output ("results/margins.txt");
    output << "Types\n" << typeCounts.transpose() << std::endl;
    output << "POS\n" << posCounts.transpose()  << std::endl;
    output.close();
  }
  
  // Random projections of row and column spaces 
  startTime = clock();
  Matrix RP (B.rows(), 2*nProjections);
  Matrix rightR = Matrix::Random(B.cols(), nProjections);
  Matrix leftR  = Matrix::Random(B.rows(), nProjections);
  RP.leftCols (nProjections) = B             * rightR;
  RP.rightCols(nProjections) = B.transpose() * leftR;
  ss << "Compute random projection RP[" << RP.rows() << "x" << RP.cols() << "]";
  print_time(ss.str(), startTime, clock());
  ss.str("");

  // Repeat for validation data if present
  Matrix vRP (V.rows(), 2*nProjections);
  if(v.rows() > 0)
  { startTime = clock();
    vRP.leftCols (nProjections) = V.leftCols(B.cols())            * rightR;
    vRP.rightCols(nProjections) = V.transpose().leftCols(B.cols())* leftR;
    ss << "Compute validation random projection RP[" << vRP.rows() << "x" << vRP.cols() << "]";
    print_time(ss.str(), startTime, clock());
    ss.str("");
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
    bool useL2      ('2' == distance);
    bool useScaling (scaling != 0);
    KMeansClusters clusters (RP, wts, useL2, useScaling, nClusters, nIterations);
    tags = clusters.cluster_tags();
  }
  ss << "Compute " << nClusters << " cluster centers.";
  print_time(ss.str(), startTime, clock());
  ss.str("");

  // Cluster validation data if present
  
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
      std::cout << std::setw(width) << indexPos[col];
    std::cout << std::endl << crossTab << std::endl;
    /*
    for(int col=0; col<crossTab.cols(); ++col)
      std::clog << "Column margin in crosstab: " << indexPos[col] << "  " << posMap[indexPos[col]] << " == " << crossTab.col(col).sum() << std::endl;
    */
    int nCorrect = crossTab.rowwise().maxCoeff().sum();
    ss << "Classify " << nCorrect << " correctly out of " << tokenManager.input_length()
       << "(" << 100.0 * ((float) nCorrect)/tokenManager.input_length() << "%)\n";
    std::clog << ss.str();
    std::cout << ss.str();
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
parse_arguments(int argc, char** argv,
		int &nSkip, float &threshold,
		int &nProj, int &scaling, char &distance, int &weighting, int &nClus, int &nIter, int &nPrint,
		std::string &vFileName)
{
  static struct option long_options[] = {
    // bigram prep options
    {"skip",         required_argument, 0, 'k'},
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
    {"validation",   required_argument, 0, 'v'},
    {0, 0, 0, 0}                             // terminator 
  };
  int key;
  int option_index = 0;
  while (-1 !=(key = getopt_long (argc, argv, "k:t:r:s:w:c:n:p:v:", long_options, &option_index))) // colon means has argument
  {
    // std::cout << "Option key " << char(key) << " for option " << long_options[option_index].name << ", option_index=" << option_index << std::endl;
    switch (key)
    {
    case 'k' :
      {
	nSkip = read_utils::lexical_cast<int>(optarg);
	break;
      }
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
    case 'v' : 
      {
	vFileName = optarg;
	break;
      }
    default:
      {
	std::cout << "PARSE: Option not recognized; returning.\n";
      }
    } // switch
  } // while
}
