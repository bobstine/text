#include "bigram.h"

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
fill_sparse_bigram_matrix(SparseMatrix &B, int skip, TokenManager const& tmRow, TokenManager const& tmCol, bool transpose=false)
{
  std::map<std::pair<int,int>,int> bgramMap;
  tmRow.fill_bigram_map(bgramMap, skip, tmCol, transpose);    
  fill_sparse_matrix(B, bgramMap);
}



void
parse_arguments(int argc, char** argv,
		int &nSkip, float &threshold, bool &bidirectional, int &nProj, bool &scaleData,
		bool &weightCentroid, int &nClus, bool &scaleCentroid, int &nIter, int &nPrint,
		std::string &vFileName);
  
int main(int argc, char **argv)
{
  using std::string;
  std::ostringstream ss;

  // argument defaults
  int     nSkip          =     0;       // skipped words in bigram; 0 is usual bigram
  float   posThreshold   =   0.0;       // accum all pos tags with relative frequency below this threshold
  bool    bidirectional  = false;
  int     nProjections   =    30;
  bool    scaleData      = false;
  bool    weightCentroid = false;
  int     nClusters      =   100;
  bool    scaleCentroid  = false; 
  int     nIterations    =    10;
  int     nPrint         =     0;
  string  vFileName      =    "";

  parse_arguments(argc, argv, nSkip, posThreshold, bidirectional, nProjections, scaleData, weightCentroid, nClusters, scaleCentroid, nIterations, nPrint, vFileName);
  std::cout << "MAIN: Arguments skip=" << nSkip << " threshold=" << posThreshold << " bidirectional=" << bidirectional
	    << " projections=" << nProjections << " scale_data=" << scaleData << " weight_centroid=" << weightCentroid
	    << " clusters=" << nClusters << " iterations=" << nIterations
	    << " print=" << nPrint << " validation=" << vFileName << std::endl;

  const bool useValidation (vFileName.size() > 0);
  
  // read tagged tokens and pos from input
  clock_t startTime = clock();

  // use for debug running gdb
  //   std::ifstream is("/Users/bob/C/text/tagged/train.tagged");
  //   TokenManager tokenManager (is, posThreshold);

  TokenManager tokenManager = TokenManager(std::cin, posThreshold);
  
  print_time("Read tokens from cin, sort and assign IDs in TokenManager.", startTime, clock());
  if (nPrint) tokenManager.print_type_tags(nPrint);
  int nAmbiguous = tokenManager.n_ambiguous();
  ss << "MAIN: Source has " << nAmbiguous << " ambiguous tokens among " << tokenManager.input_length()
     << "  (" << 100.0*(1.0 - ((float)nAmbiguous)/tokenManager.input_length()) << "% pure)" << std::endl;
  std::cout << ss.str();
  std::clog << ss.str();
  ss.str("");
  tokenManager.write_frequencies_to_file("results/margins.txt");
  
  // build the sparse bigram array
  const int nTypes (tokenManager.n_types());
  SparseMatrix B(nTypes,nTypes);
  startTime = clock();
  fill_sparse_bigram_matrix(B, nSkip, tokenManager, tokenManager);
  std::clog << "MAIN: Initial bigram row sums are "
	    << B.row(0).sum() << " "<< B.row(1).sum() << " " << B.row(2).sum() << " "
	    << B.row(3).sum() << " " << B.row(4).sum() << " " << std::endl;
  ss << "Init sparse bigram B[" << B.rows() << "x" << B.cols() << "] from map; sum +/,B= " << B.sum();
  print_time(ss.str(), startTime, clock());
  ss.str("");
  
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
  if (weightCentroid)  
  { for(int i=0; i<B.rows(); ++i)
    { wts(i) = B.row(i).sum();
      if (0 == wts(i)) std::clog << "MAIN: row " << i << " of B sums to zero.\n";
    }
  }  
  KMeansClusters clusters(RP, bidirectional, wts, scaleData, nClusters, scaleCentroid, nIterations);
  clusters.print_summary_stats(std::clog);
  clusters.print_to_stream(std::clog);
  clusters.summarize_rare_cases(std::clog);
  {
    std::ofstream deskTopFile ("/Users/bob/Desktop/centers.txt");
    deskTopFile << "    n    Mean    Max    Norm ";  // space before name for JMP
    for(int i=0; i<clusters.number_of_clusters(); ++i)
      deskTopFile << " D" << i;
    deskTopFile << std::endl << clusters.within_cluster_summary_stats() << std::endl;
  }
  
  ClusterClassifier classifier(clusters, tokenManager);
  {
    ConfusionMatrix table = make_confusion_matrix (classifier, tokenManager);
    table.print_to_stream(std::cout);
    table.print_to_stream(std::clog);
  }

  // optional validation which has same column indices as B
  if (false) // (useValidation)
  { std::clog << "\n\nMAIN: Running validation.\n";
    TokenManager  validationTM(vFileName, posThreshold);
    int nOOV = validationTM.n_types_oov(tokenManager);
    std::clog << "MAIN: Validation data has " << validationTM.n_types() << " types (" << nOOV << " OOV) and " << validationTM.n_POS() << " POS." << std::endl;

    SparseMatrix V, Vt;
    startTime = clock();
    V.resize(validationTM.n_types(), tokenManager.n_types());
    fill_sparse_bigram_matrix(V, nSkip, validationTM, tokenManager);
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
    {
      ConfusionMatrix table = make_confusion_matrix (classifier, validationTM);
      table.print_to_stream(std::cout);
      table.print_to_stream(std::clog);
    }
    // repeat for just oov
    /*
      CrossTab oovCrossTab(nClusters, validationTM.n_POS());
      for(auto it=validationTM.token_list_begin(); it != validationTM.token_list_end(); ++it)
      if(!tokenManager.known_type(it->first))
      oovCrossTab.increment(vTags[validationTM.index_of_type(it->first)], validationTM.index_of_POS(it->second));
      oovCrossTab.print_accuracy_to_stream(std::clog);
      oovCrossTab.print_accuracy_to_stream(std::cout);
      oovCrossTab.print_to_stream(std::cout);
    */
  }
  
  return 0;

}



void
parse_arguments(int argc, char** argv,
		int &nSkip, float &threshold, bool &bidirectional, int &nProj, bool &scaleData,
		bool &weightCentroid, int &nClus, bool &scaleCentroid, int &nIter, int &nPrint,
		std::string &vFileName)
{
  static struct option long_options[] = {
    // bigram prep options
    {"skip",         required_argument, 0, 'k'},
    {"threshold",    required_argument, 0, 't'},
    {"bidirectional",      no_argument, 0, 'b'},
    {"projections",  required_argument, 0, 'r'},
    {"scale_data",         no_argument, 0, 's'},
    // clustering options
    {"weight_centroid",    no_argument, 0, 'w'},
    {"clusters",     required_argument, 0, 'c'},
    {"scale_centroid",     no_argument, 0, 'C'},
    {"iterations",   required_argument, 0, 'n'},
    // misc
    {"print",        required_argument, 0, 'p'},
    {"validation",   required_argument, 0, 'v'},
    {0, 0, 0, 0}                             // terminator 
  };
  int key;
  int option_index = 0;
  while (-1 !=(key = getopt_long (argc, argv, "k:t:br:swc:Cn:p:v:", long_options, &option_index))) // colon means has argument
  {
    // std::cout << "Option key " << char(key) << " for option " << long_options[option_index].name << ", option_index=" << option_index << std::endl;
    switch (key)
    {
    case 'k' : { nSkip          = read_utils::lexical_cast<int>(optarg);   break; }
    case 't' : { threshold      = read_utils::lexical_cast<float>(optarg); break; }
    case 'b' : { bidirectional  = true;                                    break; }
    case 'r' : { nProj          = read_utils::lexical_cast<int>(optarg);   break; }
    case 's' : { scaleData      = true;                                    break; }
    case 'w' : { weightCentroid = true;                                    break; }
    case 'c' : { nClus          = read_utils::lexical_cast<int>(optarg);   break; }
    case 'C' : { scaleCentroid  = true;                                    break; }
    case 'n' : { nIter          = read_utils::lexical_cast<int>(optarg);   break; }
    case 'p' : { nPrint         = read_utils::lexical_cast<int>(optarg);   break; }
    case 'v' : { vFileName      = optarg;                                  break; }
    default  : { std::cout << "PARSE: Option not recognized; returning.\n";       }
    } // switch
  } 
}
