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


class CrossTabulator
{
  using std::string;
  
  int mN;
  std::vector<string> mColLabels;
  std::vector<string> mRowLabels;
  Eigen::MatrixXi mTable;

public:
  
  CrossTabulator (int nRows, int nCols) : mN(0), mColLabels(), mRowLabels(), mTable(Eigen::MatrixXi::Zero(nRows,nCols))  {}
  
  template <class ItR, class ItC>
  CrossTabulator (ItR rBegin, ItR rEnd, ItC cBegin, ItC cEnd) : mN(0), mColLabels(), mRowLabels(), mTable(Eigen::MatrixXi::Zero(nRows,nCols))
    { int nR = 0; int nC = 0;
      for(ItR it=rBegin; it != rEnd; ++it)
      { ++nR;
	mRowLabels.push_back(*it);
      }
      for(ItC it=cBegin; it != cEnd; ++it)
      { 

  void            increment(int i, int j)        { ++mN; ++mTable(i,j); }
  
  int             element(int i, int j)    const { return mTable(i,j); }
  Eigen::VectorXi row_margins()            const { return mTable.rowwise().sum(); }
  Eigen::VectorXi col_margins()            const { return mTable.colwise().sum(); }
  int             total()                  const { return mN; }

  float           accuracy()               const { return 100.0*((float)sum_row_max())/mN; }
  int             sum_row_max()            const { return mTable.rowwise().maxCoeff().sum(); }

  void            print_to_stream (std::ostream &os) const
    { 
      os << "Classify " << sum_row_max() << " correctly out of " << mN << "(" << accuracy() << "%)\n";
      int width = 2+log10(mTable.array().maxCoeff());
      Eigen::VectorXi colMargins = col_margins();
      for(int col=0; col<mTable.cols(); ++col)
	os << std::setw(width) << colMargins(col);
      os << std::endl << mTable << std::endl;
    }
};



void
print_time(std::string const& s, clock_t const& start, clock_t const& stop)
{
  clock_t diff (stop-start);
  std::cout << "TIME: " << s << " [" << ((float)diff)/CLOCKS_PER_SEC << " sec]\n";
  std::clog << "TIME: " << s << " [" << ((float)diff)/CLOCKS_PER_SEC << " sec]\n";
}


void
fill_sparse_bigram_matrix(SparseMatrix &B, int skip, TokenManager const& tmRow, TokenManager const& tmCol)
{
  std::map<std::pair<int,int>,int> bgramMap;
  if(&tmRow == &tmCol)
    tmRow.fill_bigram_map(bgramMap, skip);
  else
    tmRow.fill_bigram_map(bgramMap, skip, tmCol);    
  fill_sparse_matrix(B, bgramMap);
}


void
parse_arguments(int argc, char** argv,
		int &nSkip, float &posThreshold,
		int &nProj, int &scaling, char &dist, int &weighting, int &nClus, int &nIter, int &nPrint,
		string &validationFileName);

  
int main(int argc, char **argv)
{
  using std::string;
  std::ostringstream ss;

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

  const bool useValidation (vFileName.size() > 0);
  
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
    tokenManager.write_frequencies_to_file("results/margins.txt");
  }
  
  // build the sparse bigram array using integer id for words
  const int nTypes (tokenManager.n_types());
  SparseMatrix B(nTypes,nTypes);
  {
    startTime = clock();
    fill_sparse_bigram_matrix(B, nSkip, tokenManager, tokenManager);
    ss << "Init sparse bigram B[" << B.rows() << "x" << B.cols() << "] from map; sum +/,B= " << B.sum() << "  +/B[452,]=" << B.row(452).sum();
    print_time(ss.str(), startTime, clock());
    ss.str("");
  }

  // optional validation which has same column indices as B
  SparseMatrix V(0,0);
  TokenManager validationTM;
  if (useValidation)
  { startTime = clock();
    validationTM = TokenManager(vFileName, posThreshold);
    int nOOV = validationTM.n_types_oov(tokenManager);
    std::clog << "MAIN: Validation has " << validationTM.n_types() << " with OOV count " << nOOV << " types." << std::endl;
    V.resize(validationTM.n_types(), tokenManager.n_types());
    fill_sparse_bigram_matrix(V, nSkip, validationTM, tokenManager);
    ss << "Init validation sparse bigram V[" << V.rows() << "x" << V.cols() << "] from map; sum +/,V=" << V.sum();
    print_time(ss.str(), startTime, clock());
    ss.str("");
  }
  
  // Random projections of row and column spaces 
  startTime = clock();
  Matrix rightR = Matrix::Random(B.cols(), nProjections);
  Matrix RP = B * rightR;
  // RP.rightCols(nProjections) = B.transpose() * rightR;
  // Matrix leftR  = Matrix::Random(B.rows(), nProjections);
  // RP.leftCols (nProjections) = B             * leftR;
  ss << "Compute random projection RP[" << RP.rows() << "x" << RP.cols() << "]";
  print_time(ss.str(), startTime, clock());
  ss.str("");

  // Repeat for validation data if present
  Matrix vRP (V.rows(), nProjections); 
  if(useValidation)
  { startTime = clock();
    vRP = V * rightR;
    ss << "Compute validation random projection vRP[" << vRP.rows() << "x" << vRP.cols() << "]";
    print_time(ss.str(), startTime, clock());
    ss.str("");
  }
  
  // cluster tokens using k-means
  startTime = clock();
  std::vector<int> tags, vTags;
  {
    IntVector wts = IntVector::Ones(B.rows());
    if (weighting)
    { wts.resize(B.rows());
      for(int i=0; i<B.rows(); ++i)
      { wts(i) = B.row(i).sum();
	if (0 == wts(i))
	  std::clog << "MAIN: row " << i << " of B sums to zero.\n";
      }
      std::clog << "MAIN: First five weights are " << wts.head(5).transpose() << std::endl;  // should match next line
      std::clog << "MAIN: First five tag freqs   "; for(int i=0; i<5; ++i) std::clog << tokenManager.type_freq(i) << " "; std::clog << std::endl;
    }
    bool useL2      ('2' == distance);
    bool useScaling (scaling != 0);
    KMeansClusters clusters (RP, wts, useL2, useScaling, nClusters, nIterations);
    tags = clusters.cluster_tags();
    if(useValidation)
      vTags = clusters.assign_to_clusters(&vRP);
  }
  ss << "Compute " << nClusters << " cluster centers and optionally assign validation clusters.";
  print_time(ss.str(), startTime, clock());
  ss.str("");
  
  
  CrossTabulator crossTab(nClusters, tokenManager.n_POS());
  for(auto it=tokenManager.token_list_begin(); it != tokenManager.token_list_end(); ++it)
    crossTab.increment(tags[tokenManager.index_of_type(it->first)], tokenManager.index_of_POS(it->second));
  crossTab.print_to_stream(std::clog);

  if(useValidation)
  { CrossTabulator vCrossTab(nClusters, tokenManager.n_POS());
    for(auto it=validationTM.token_list_begin(); it != validationTM.token_list_end(); ++it)
      vCrossTab.increment(tags[tokenManager.index_of_type(it->first)], tokenManager.index_of_POS(it->second));
    vCrossTab.print_to_stream(std::clog);
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
    print_time("Compute SVD and extract U*D.", startTime, clock());
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
