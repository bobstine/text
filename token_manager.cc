#include "token_manager.h"

#include <algorithm>
#include <set>
#include <iostream>
#include <fstream>
#include <ios>


std::string messageTag = "TKMG: ";


void
TokenManager::init_from_file(std::string &fileName, float posThreshold)
{
  std::ifstream input (fileName.c_str());
  if (! input)
    std::cerr << messageTag << "Cannot open input token file " << fileName << std::endl;
  else
    init_from_stream(input, posThreshold);
}


void
TokenManager::init_from_stream(std::istream &input, float posThreshold)
{
  std::map<POS,int> posMap;   // initial pos prior to collapse rare pos to OTH
  while (input)
  { string typeStr;
    input >> typeStr;
    if (typeStr.empty()) break;   // end of file
    Type type(typeStr);
    ++mTypeMap[type];
    string posStr;
    input >> posStr;
    POS pos(posStr);
    ++posMap[pos];
    mTokens.push_back(std::make_pair(type,pos));
  }
  std::clog << messageTag << "Read " << n_types() << " types from input of " << mTokens.size() << " strings." << std::endl;
  std::clog << messageTag << "Identified " << posMap.size() << " distinct POS." << std::endl;
  if(posThreshold > 0.0)                                             // remove rare pos tags
  { std::map<POS,bool> reduce;
    int reduceCount = 0;
    for(auto it=posMap.cbegin(); it != posMap.cend(); ++it)
    { float freq = ((float)it->second)/mTokens.size();
      if (freq < posThreshold)
      { reduce[it->first] = true; ++reduceCount; }
      else reduce[it->first]=false;
    }
    if (reduceCount > 0)
    { std::clog << messageTag << "Mapped " << reduceCount << " POS categories to OTH." << std::endl;
      for(auto it=mTokens.begin(); it != mTokens.end(); ++it)
	if (reduce[it->second]) it->second = POS("OTH");
    }
  }
  for(auto it=mTokens.cbegin(); it != mTokens.cend(); ++it)          // build pos map for each type
  { ++mPOSMap[it->second];
    ++mTypePOSMap[it->first][it->second];
  }
  int posIndex = 0;
  std::map<POS,int> newMap;
  for(std::map<POS,int>::iterator it=mPOSMap.begin(); it != mPOSMap.end(); ++it)          // define index for each pos
  { POS pos(it->first.pos_as_string());
    pos.assign_index(posIndex);   // compiler would not allow assignment of index (const)
    newMap[pos] = it->second; 
    mPOSVec.push_back(it->first);
  }
  mPOSMap = newMap;
  std::multimap<int,Type> sortedTypeMap;                               // sort types by frequency by inserting into multimap
  for (auto it = mTypeMap.begin(); it != mTypeMap.end(); ++it)
    sortedTypeMap.insert( std::make_pair(-it->second,it->first) );    // negate so decreasing
  int tokenID = 0;
  mTypeVec.resize(n_types(), Type("dummy"));
  for(auto it = sortedTypeMap.begin(); it != sortedTypeMap.end(); ++it)
  { it->second.assign_index(tokenID);
    mTypeIndexMap[it->second] = tokenID;
    mTypeVec[tokenID] = it->second;
    ++tokenID;
  }
}

int
TokenManager::index_of_type (Type const& t) const
{
  if(mTypeIndexMap.count(t))
    return mTypeIndexMap.at(t);
  else
  { assert(false);  // type not found
    return -1;
  }
}


TokenManager::TypeVector
TokenManager::type_vector()        const
{
  TypeVector sv;

  for(auto it=mTypeMap.cbegin(); it != mTypeMap.cend(); ++it)
    sv.push_back(it->first);
  return sv;
}


TokenManager::POSVector
TokenManager::POS_vector()        const
{
  POSVector sv;

  for(auto it=mPOSMap.cbegin(); it != mPOSMap.cend(); ++it)
    sv.push_back(it->first);
  return sv;
}


int
TokenManager::n_ambiguous ()       const
{
  int nAmbiguous = 0;
  for(auto it=mTypePOSMap.cbegin(); it!=mTypePOSMap.cend(); ++it)
  { Type type = it->first;
    int max=0;
    for(auto it2=it->second.cbegin(); it2 != it->second.cend(); ++it2)
      if (it2->second > max) max = it2->second;
    nAmbiguous += mTypeMap.at(type)-max;
  }
  return nAmbiguous;
}


POS
TokenManager::POS_of_type (Type const& t)      const
{
  const std::map<POS,int> posMap = mTypePOSMap.at(t);
  if(posMap.empty())
  { std::cerr << messageTag << "Type " << t << " does not have POS recorded. " << std::endl;
    return POS("");
  }
  else
  { POS pos   = posMap.cbegin()->first;
    int count = posMap.cbegin()->second;
    for(auto it = ++posMap.cbegin(); it != posMap.cend(); ++it)
      if(it->second > count)
      { count = it->second;
	pos   = it->first;
      }
    return pos;
  }
}


std::vector< std::pair<POS,int> >
TokenManager::type_POS_tags (Type const& type, bool sort) const
{
  typedef std::pair<POS,int> PSI;
  std::vector<PSI> counts;
  const std::map<POS,int> posMap = mTypePOSMap.at(type);
  if(posMap.empty())
    std::cerr << messageTag << "Token " << type << " does not have POS recorded. " << std::endl;
  else
  { for(auto it = posMap.cbegin(); it !=posMap.cend(); ++it)
      counts.push_back(*it);
    if(sort)
      std::sort(counts.begin(), counts.end(), [](PSI const& a, PSI const& b) { return a.second < b.second; });
  }
  return counts;
}


int
TokenManager::n_types_oov(TokenManager const& tm) const
{
  int nOOV = 0;
  for (auto it = mTypeMap.cbegin(); it != mTypeMap.cend(); ++it)
    if(tm.type_freq(it->first) == 0) ++ nOOV;
  return nOOV;
}



void
TokenManager::fill_bigram_map(BigramMap &bm, int skip, TokenManager const& tm, bool transpose) const
{
  std::pair<int,int> ij;                 // beyond those indices for tokens found in tm
  auto itBack = tm.token_list_begin();
  for (int i=0; i<skip+1; ++i) ++itBack;
  int nNotFound = 0;
  if (!transpose)
  { for (auto it=mTokens.cbegin(); itBack != mTokens.cend(); ++it, ++itBack)
    { int i = it->first.index();             // i from this
      int j = itBack->first.index();         // j from index set used in tm
      if (0 <= j)                            // valid
      { ij.first = i;
	ij.second = j;
	++bm[ij];
      }
      else ++nNotFound;
    }
  }
  else 
  { for (auto it=mTokens.cbegin(); itBack != mTokens.cend(); ++it, ++itBack)
    { int i = itBack->first.index();         // i from tm
      int j = it->first.index();             // j from this
      if (0 <= j)                            // valid
      { ij.first = i;
	ij.second = j;
	++bm[ij];
      }
      else ++nNotFound;
    }
  }
  if (nNotFound)
    std::clog << messageTag << " *** WARNING *** " << nNotFound << " tokens not found in bigram lookup." << std::endl;
}
    
void
TokenManager::print_tags(int num) const
{
  for(int k=0; k<num; ++k)
  { POSCountVector v = type_POS_tags(mTypeVec[k],true);
    std::clog << "k=" << k << "    Type   '" << mTypeVec[k] << "'" << std::endl;
    for (auto it=v.begin(); it != v.end(); ++it)
      std::clog << "     " << it->first << " " << it->second << std::endl;
  }
}

void
TokenManager::write_frequencies_to_file(std::string fileName) const   // write type and counts 
{
  std::map<POS, int> posMap (POS_map());
  std::vector<int> posCounts (posMap.size());
  int i = 0;
  for(auto it=posMap.cbegin(); it != posMap.cend(); ++it)
  { posCounts[i++] = it->second; }
  std::ofstream output (fileName);
  output << "Types\n";
  for (auto it=mTypeMap.cbegin(); it != mTypeMap.cend(); ++it)
    output << it->second << ", ";
  output << "\nPOS\n";
  for (auto it=posCounts.cbegin(); it != posCounts.cend(); ++it)
    output << *it << ", ";
  output << std::endl;
  output.close();
}

  
void
TokenManager::print_to_stream (std::ostream&os) const
{
  os << messageTag << "Token manager read " << input_length() << " <token,POS> pairs, finding "
     << n_POS() << " parts of speech among " << n_types() << " unique types ("
     << n_ambiguous() << " ambiguous).\n";
  os << " POS Counts among input tokens are : \n";
  for(auto it=mPOSMap.cbegin(); it != mPOSMap.cend(); ++it)
    std::cout << "       " << it->first << " " << it->second << std::endl;
}
