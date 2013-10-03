#include "vocabulary.h"
#include "eigen_utils.h"
#include "regex.h"            // parse_numeric_string

#include <assert.h>
#include <fstream>
#include <sstream>

//     globals     globals     globals     globals     globals     globals     globals     globals     globals

static std::string messageTag ("VOCB: ");

static std::string oovStr ("OOV");
Type Vocabulary::OOV = Type(oovStr);

static std::string eolStr ("EOL");     // marks end of document/line
Type Vocabulary::EOL = Type(eolStr);


//     print     print     print     print     print     print     print     print     print     print     print     

void
Vocabulary::print_to_stream(std::ostream & os) const
{
  os << "Vocabulary has " << n_types() << " types from " << mNTokens << " tokens, with "
     << mFreqMap.at(OOV) << " OOV.  Most frequent 10 types are:\n      ";
  for (int i=0; i<10; ++i)
  { assert (i == mIndexMap.at(mTypeVector[i]));
    os << " " << mTypeVector[i] << "->" << mFreqMap.at(mTypeVector[i]);
  }
}

void
Vocabulary::write_type_freq(std::ostream & os) const
{
  for(auto x : mFreqMap)
    os << x.second << " ";
}


//     init     init     init     init     init     init     init     init     init     init     init     init     

void
Vocabulary::init_from_file(std::string fileName)
{
  std::ifstream inputFile (fileName);
  if(!inputFile)
    std::cerr << messageTag << " ******  ERROR  ****** Input file " << fileName << " not found; exiting." << std::endl;
  else
    init_from_stream (inputFile);
}

void
Vocabulary::parse_line(std::string const& line, std::map<Type,int> &vocab)
{
   std::stringstream ss(line);
   std::string token;
   for(int i=0; i<mSkipInitial; ++i)
     ss >> token;
   while (ss >> token)
   { if (can_parse_as_numeric_string(token))   // converts digits, dollar amts, phone nos
       token = parse_numeric_string(token);
     ++mNTokens;
     mTokens.push_back(Type(token));
     ++vocab[Type(token)];
   }
   if (mMarkEOL)
   { ++mNTokens;
     mTokens.push_back(EOL);
     ++vocab[Type(EOL)];
   }
}


void
Vocabulary::init_from_stream(std::istream& is)
{
  std::map<Type,int> fullVocab;
  fullVocab[OOV] = 0;                         // insert oov

  int nLines = 0;
  std::string line;
  while(getline(is, line))
  { ++nLines;
    parse_line(line, fullVocab);
  }
  /*
  std::string token;
  while (is >> token)
  { if (can_parse_as_numeric_string(token))   // converts digits, dollar amts, phone nos
      token = parse_numeric_string(token);
    ++mNTokens;
    mTokens.push_back(Type(token));
    ++fullVocab[Type(token)];
  }
   */
  std::clog << messageTag << "Full vocabulary has " << fullVocab.size() << " types from input of " << mNTokens << " tokens on "
	    << nLines << " lines. (Skip " << mSkipInitial << " initial tokens with EOL mark " << mMarkEOL << ")." << std::endl;
  if(false) // write counts to file
  { std::clog << messageTag << "Writing word counts to file\n";
    std::ofstream os("/Users/bob/Desktop/word_counts.txt");
    for(auto x: fullVocab)
      os << x.second << std::endl;
  }
  for (auto x : fullVocab)
  { if (x.second < mMinFrequency)             // mark as oov
    { mOOVMap.insert(x);
      mFreqMap[OOV] += x.second;
    }
    else mFreqMap.insert(x);                  // transfer
  }
  int checkSum (0);
  for (auto x: mFreqMap)
    checkSum += x.second;
  std::clog << messageTag
	    << "Thresholded vocabulary of " << mFreqMap.size() << " types with token count " << checkSum << " tokens." << std::endl;
  std::multimap<int,Type> sortTypesMap;                               // sort types in frequency order
  for(auto x : mFreqMap)
    sortTypesMap.insert( std::make_pair(x.second, x.first) );
  int index = 0;
  for(auto it=sortTypesMap.crbegin(); it != sortTypesMap.crend(); ++it, ++index) // reverse iter
  { mTypeVector.push_back(it->second);
    mIndexMap[it->second] = index;
  }
  assert(n_types() == (int) mTypeVector.size());
  mOOVIndex = mIndexMap[OOV];
  std::clog << messageTag << "Position of OOV type is " << mOOVIndex << " with frequency " << mFreqMap[OOV] << std::endl;
}
  

//     access     access     access     access     access     access     access     access     access     access     access

int
Vocabulary::type_index(Type const& type) const
{
  auto it = mIndexMap.find(type);
  if(it != mIndexMap.end())
    return it->second;
  else
    return mOOVIndex;
}

//     bigram     bigram     bigram     bigram     bigram     bigram     bigram     bigram     bigram     bigram     bigram

void
Vocabulary::fill_sparse_bigram_matrix(Vocabulary::SparseMatrix &B, int skip) const
{
  BigramMap bgramMap;
  fill_bigram_map(bgramMap, skip); 
  EigenUtils::fill_sparse_matrix(B, bgramMap);
}


void
Vocabulary::fill_bigram_map(BigramMap &bm, int skip) const
{
  std::pair<int,int> ij;                // beyond those indices for tokens found in tm
  auto itBack = mTokens.cbegin();
  for (int i=0; i<skip+1; ++i) ++itBack;
  for (auto it=mTokens.cbegin(); itBack != mTokens.cend(); ++it, ++itBack)
  { int i = type_index(*it);       
    assert (i >= 0);
    int j = type_index(*itBack);   
    assert (j >= 0);
    ij.first = i;
    ij.second = j;
    ++bm[ij];
  }
}


//     build count vector     build count vector     build count vector     build count vector     build count vector

void
Vocabulary::fill_sparse_regr_design (Vector &Y, Vocabulary::SparseMatrix &X, std::istream &is) const
{
  typedef Eigen::Triplet<int> T;
  std::vector<T> triplets;         // list worked in Eigen 3.1.3
  std::string line;
  std::string token;
  assert (Y.size() == X.rows());
  for (int i=0; i<X.rows(); ++i)
  { getline(is, line);
    std::istringstream is(line);
    std::map<Type,int> counts;
    is >> Y(i);                     // read y-value at head of line
    while (is >> token)
      ++counts[ Type(token) ];      // map of remaining text tokens
    for (auto x : counts)
      triplets.push_back(T(i, type_index(x.first), x.second)); // convert to indices
  }
  X.setFromTriplets(triplets.begin(), triplets.end());
}

/*  // Extra code for debugging; all agreed when testing
      {
      IntVector counts = IntVector::Zero(n_types());
      std::istringstream is(line);
      std::string token;
      std::clog << messageTag << "Read response y(" << i << ") = " << Y(i) << std::endl;
      while (is >> token)
      ++counts[ type_index(Type(token)) ];
      std::clog << messageTag << "Size of map is " << counts.size() << std::endl;
      std::clog << "Direct calc for first line has sum " << counts.sum() << " with leading 100 elements " << counts.transpose().head(100) << std::endl;
      }
*/


