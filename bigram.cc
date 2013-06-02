#include "bigram.h"

typedef Eigen::VectorXf                   Vector;
typedef Eigen::VectorXi                   IntVector;
typedef Eigen::MatrixXf                   Matrix;

typedef Eigen::SparseMatrix<float,Eigen::RowMajor> SparseMatrix;



CrossTab
build_cross_tab (KMeansClusters::Iterator begin, int nClusters, TokenManager const& tm)
{
  std::vector<string> labels;
  std::ostringstream ss;
  for (int i=1; i<=nClusters; ++i)
  { ss << "Clus " << i; labels.push_back(ss.str()); ss.str(""); }
  CrossTab table(labels.cbegin(), labels.cend(), tm.POS_begin(), tm.POS_end());
  for(auto it=tm.token_list_begin(); it != tm.token_list_end(); ++it)
    table.increment(*(begin+tm.index_of_type(it->first)), tm.index_of_POS(it->second));
  return table;
}


void
print_time(std::string const& s, clock_t const& start, clock_t const& stop)
{
  clock_t diff (stop-start);
  std::cout << "TIME: " << s << " [" << ((float)diff)/CLOCKS_PER_SEC << " sec]\n";
  std::clog << "TIME: " << s << " [" << ((float)diff)/CLOCKS_PER_SEC << " sec]\n";
}


void
fill_sparse_bigram_matrix(SparseMatrix &B, int skip, TokenManager const& tmRow, TokenManager const& tmCol, bool transpose=false)
{
  std::map<std::pair<int,int>,int> bgramMap;
  tmRow.fill_bigram_map(bgramMap, skip, tmCol, transpose);    
  fill_sparse_matrix(B, bgramMap);
}


void
parse_arguments(int argc, char** argv,
		int &nSkip, float &posThreshold, bool &bidirectional,
		int &nProj, int &scaling, char &dist, int &weighting, int &nClus, int &nIter, int &nPrint,
		string &validationFileName);

  
int main(int argc, char **argv)
{
  using std::string;
  std::ostringstream ss;

  // argument defaults
  int     nSkip         =  0;         // skipped words in bigram
  float   posThreshold  = 0.0;        // accum all pos tags with relative frequency below this threshold
  bool    bidirectional = false;
  int     nProjections  = 30;
  int     scaling       =  0;
  char    distance      ='2';         // identify distance to use  (2 for L2, c for cosine)
  int     weighting     =  1;
  int     nClusters     = 10;
  int     nIterations   = 10;
  int     nPrint        =  0;
  string  vFileName     = "";

  parse_arguments(argc, argv, nSkip, posThreshold, bidirectional, nProjections, scaling, distance, weighting, nClusters, nIterations, nPrint, vFileName);
  std::cout << "MAIN: Arguments skip=" << nSkip << " threshold=" << posThreshold << " bidirectional=" << bidirectional
	    << " nProjections=" << nProjections << " scaling=" << scaling
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
  
  // Random projections of row and column spaces 
  startTime = clock();
  Matrix RP, rightR, leftR;
  if (!bidirectional)
  { rightR = Matrix::Random(B.cols(), nProjections);
    RP = B * rightR;
  }
  else
  { assert (nProjections % 2 == 0);
    RP.resize(B.rows(),nProjections);
    rightR = Matrix::Random(B.cols(), nProjections/2);
    RP.rightCols(nProjections/2) = B * rightR;
    leftR = Matrix::Random(B.cols(), nProjections/2);
    RP.leftCols(nProjections/2) = B.transpose() * leftR;
  }
  ss << "Compute random projection RP[" << RP.rows() << "x" << RP.cols() << "]";
  print_time(ss.str(), startTime, clock());
  ss.str("");

  // cluster types using k-means
  startTime = clock();
  IntVector wts = IntVector::Ones(B.rows());
  if (weighting)  
  { for(int i=0; i<B.rows(); ++i)
    { wts(i) = B.row(i).sum();  // equiv to tokenManager.type_freq(i) so use as check
      if (0 == wts(i)) std::clog << "MAIN: row " << i << " of B sums to zero.\n";
    }
  }
  bool useL2      ('2' == distance);
  bool useScaling (scaling != 0);
  KMeansClusters clusters(RP, wts, useL2, useScaling, nClusters, nIterations);
  Eigen::VectorXi estClusterPOS;
  {
    ss << "Compute " << nClusters << " cluster centers.";
    print_time(ss.str(), startTime, clock());
    ss.str("");
    CrossTab table = build_cross_tab (clusters.cluster_tags_begin(), nClusters, tokenManager);
    estClusterPOS = table.most_common_col_in_each_row();
    table.print_accuracy_to_stream(std::clog);
    table.print_accuracy_to_stream(std::cout);
    table.print_to_stream(std::cout);
  }

  // optional validation which has same column indices as B
  if (useValidation)
  { std::clog << "\n\nMAIN: Running validation.\n";
    TokenManager  validationTM(vFileName, posThreshold);
    int nOOV = validationTM.n_types_oov(tokenManager);
    std::clog << "MAIN: Validation data has " << validationTM.n_types() << " types (" << nOOV << " OOV) and " << validationTM.n_POS() << " POS." << std::endl;

    SparseMatrix V, Vt;
    startTime = clock();
    V.resize(validationTM.n_types(), tokenManager.n_types());
    fill_sparse_bigram_matrix(V , nSkip, validationTM, tokenManager);
    ss << "Init validation sparse bigram V[" << V.rows() << "x" << V.cols() << "] from map; sum +/,V=" << V.sum();
    if (bidirectional)
    { Vt.resize(validationTM.n_types(), tokenManager.n_types());
      fill_sparse_bigram_matrix(Vt, nSkip, validationTM, tokenManager, bidirectional);
      ss << "  Vt[" << Vt.rows() << "x" << Vt.cols() << "] from map; sum +/,Vt=" << Vt.sum();
    }
    print_time(ss.str(), startTime, clock());
    ss.str("");

    Matrix vRP (V.rows(), nProjections); 
    startTime = clock();
    if (!bidirectional)
      vRP = V * rightR;
    else
    { vRP.resize(V.rows(), nProjections);
      vRP.rightCols(nProjections/2) = V  * rightR;
      vRP.leftCols (nProjections/2) = Vt * leftR;
    }
    ss << "Compute validation random projection vRP[" << vRP.rows() << "x" << vRP.cols() << "]";
    print_time(ss.str(), startTime, clock());
    ss.str("");
    std::clog << "MAIN: Assigning validation data to clusters." << std::endl; 
    std::vector<int> vTags = clusters.assign_to_clusters(&vRP);
    std::clog << "MAIN: Tabulating validation data." << std::endl;
    CrossTab vCrossTab = build_cross_tab (vTags.cbegin(), nClusters, validationTM);
    
    // check max count by direct computation from input tokens
    int nRightPOS = 0;
    for(auto it=validationTM.token_list_begin(); it != validationTM.token_list_end(); ++it)
    { int cluster = vTags[validationTM.index_of_type(it->first)];
      int truPOS = tokenManager.index_of_POS(it->second);
      if(estClusterPOS[cluster] == truPOS) ++nRightPOS;
    }
    std::clog << "MAIN: Direct token calculation finds correct tags assigned to " << nRightPOS << " tokens.\n";
    vCrossTab.print_accuracy_to_stream(std::clog);   // number correct agrees with above if most common label is right label
    vCrossTab.print_accuracy_to_stream(std::cout);
    vCrossTab.print_to_stream(std::cout);

    // repeat for just oov
    CrossTab oovCrossTab(nClusters, validationTM.n_POS());
    for(auto it=validationTM.token_list_begin(); it != validationTM.token_list_end(); ++it)
      if(!tokenManager.known_type(it->first))
	oovCrossTab.increment(vTags[validationTM.index_of_type(it->first)], validationTM.index_of_POS(it->second));
    oovCrossTab.print_accuracy_to_stream(std::clog);
    oovCrossTab.print_accuracy_to_stream(std::cout);
    oovCrossTab.print_to_stream(std::cout);
  }
  
  // Write original tags and cluster ids to file
  if (false)
  { std::ios_base::openmode mode = std::ios_base::trunc;
    std::string fileName ("/Users/bob/Desktop/tags.txt");
    std::ofstream file (fileName.c_str(), mode);
    file << "Token\tPOS\tCluster" << std::endl;
    KMeansClusters::Iterator clusterIndex = clusters.cluster_tags_begin();
    for(auto it = tokenManager.token_list_begin(); it != tokenManager.token_list_end(); ++it)
      file << it->first << "\t" << it->second << "\t" << *(clusterIndex+tokenManager[it->first]) << std::endl;
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
		int &nSkip, float &threshold, bool &bidirectional,
		int &nProj, int &scaling, char &distance, int &weighting, int &nClus, int &nIter, int &nPrint,
		std::string &vFileName)
{
  static struct option long_options[] = {
    // bigram prep options
    {"skip",         required_argument, 0, 'k'},
    {"threshold",    required_argument, 0, 't'},
    {"bidirectional",      no_argument, 0, 'b'},
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
  while (-1 !=(key = getopt_long (argc, argv, "k:t:br:s:w:c:n:p:v:", long_options, &option_index))) // colon means has argument
  {
    // std::cout << "Option key " << char(key) << " for option " << long_options[option_index].name << ", option_index=" << option_index << std::endl;
    switch (key)
    {
    case 'k' :
      {
	nSkip = read_utils::lexical_cast<int>(optarg);
	break;
      }
    case 'b' :
      {
	bidirectional = true;
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
  } 
}
