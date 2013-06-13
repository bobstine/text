#include "token_manager.h"

#include <algorithm>
#include <set>
#include <iostream>
#include <fstream>
#include <ios>


static std::string messageTag = "TKMG: ";

inline size_t min(size_t a, size_t b) { return (a < b) ? a : b; }


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
    ++mTypeCountMap[type];
    string posStr;
    input >> posStr;
    POS pos(posStr);
    ++posMap[pos];
    mTokens.push_back(std::make_pair(type,pos));
  }
  std::clog << messageTag << "Read " << mTypeCountMap.size() << " types from input of " << mTokens.size() << " strings." << std::endl;
  std::clog << messageTag << "Identified " << posMap.size() << " distinct POS." << std::endl;
  if(posThreshold > 0.0) // remove rare pos tags
  { std::map<POS,bool> reduce;
    int reduceCount = 0;
    for(auto it=posMap.cbegin(); it != posMap.cend(); ++it)
    { float freq = ((float)it->second)/mTokens.size();
      if (freq < posThreshold)
      { reduce[it->first] = true;
	++reduceCount;
      }
      else reduce[it->first]=false;
    }
    if (reduceCount > 0)
    { const string other("OTH");
      std::clog << messageTag << "Mapped " << reduceCount << " POS categories to OTH." << std::endl;
      for(auto it=mTokens.begin(); it != mTokens.end(); ++it)
	if (reduce[it->second]) it->second = POS(other);
    }
  }
  for(auto it=mTokens.cbegin(); it != mTokens.cend(); ++it) // build pos map for each type
  { ++mPOSCountMap[it->second];
    ++mTypePOSMap[it->first][it->second];
  }
  std::multimap<int,Type> sortTypesMap;                          // sort types in frequency order
  for(auto it = mTypeCountMap.cbegin(); it != mTypeCountMap.cend(); ++it)
    sortTypesMap.insert( std::make_pair(it->second, it->first) );
  int index = 0;
  for(auto it=sortTypesMap.crbegin(); it != sortTypesMap.crend(); ++it, ++index) // reverse iter
  { mTypeIndexMap[it->second] = index;
    mTypeVector.push_back(it->second);
  }
  assert(mTypeVector.size()   == mTypeCountMap.size());
  assert(mTypeIndexMap.size() == mTypeCountMap.size());
  std::clog << messageTag << "Ten most common types are ";
  for(int i=0; i<(int)min(mTypeVector.size(),10); ++i)
  { assert(i == mTypeIndexMap[mTypeVector[i]]);                //  quick check of just a few
    std::clog << mTypeVector[i] << " " << mTypeCountMap[mTypeVector[i]] << "   ";
  }
  std::clog << std::endl;
}


TokenManager::POSVector
TokenManager::POS_vector()        const
{
  POSVector sv;

  for(auto it=mPOSCountMap.cbegin(); it != mPOSCountMap.cend(); ++it)
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
    nAmbiguous += mTypeCountMap.at(type)-max;
  }
  return nAmbiguous;
}


POS
TokenManager::POS_of_type (Type const& t)      const
{
  const std::map<POS,int> posMap = mTypePOSMap.at(t);
  if(posMap.empty())
  { std::cerr << messageTag << "Type " << t << " does not have POS recorded. " << std::endl;
    const string empty ("");
    return POS(empty);
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
TokenManager::POS_tags_of_type (Type const& type, bool sort) const
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
  for (auto it = mTypeCountMap.cbegin(); it != mTypeCountMap.cend(); ++it)
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
    { int i = type_index(it->first);                // i from this
      int j = tm.type_index(itBack->first);         // j from index set used in tm
      if (0 <= j)                                   // valid
      { ij.first = i;
	ij.second = j;
	++bm[ij];
      }
      else ++nNotFound;
    }
  }
  else 
  { for (auto it=mTokens.cbegin(); itBack != mTokens.cend(); ++it, ++itBack)
    { int i = tm.type_index(itBack->first);         // i from tm
      int j = type_index(it->first);                // j from this
      if (0 <= j)                                   // valid
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
TokenManager::print_type_tags(int maxToPrint) const
{
  const bool sorted (true);
  int nPrinted = 0;
  for(auto it = mTypeCountMap.cbegin(); it != mTypeCountMap.cend(); ++it, ++nPrinted)
  { if (nPrinted == maxToPrint) break;
    POSCountVector v = POS_tags_of_type(it->first, sorted);
    std::clog << "    Type   '" << it->first << "'" << std::endl;
    for (auto vit=v.begin(); vit != v.end(); ++it)
      std::clog << "     " << vit->first << " " << vit->second << std::endl;
  }
}

void
TokenManager::write_frequencies_to_file(std::string fileName) const   // write type and counts 
{
  std::vector<int> posCounts (mPOSCountMap.size());
  int i = 0;
  for(auto it=mPOSCountMap.cbegin(); it != mPOSCountMap.cend(); ++it)
    posCounts[i++] = it->second; 
  std::ofstream output (fileName);
  output << "Types\n";
  for (auto it=mTypeCountMap.cbegin(); it != mTypeCountMap.cend(); ++it)
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
  os << messageTag << "Token manager read " << input_length() << " <word,POS> pairs, finding "
     << n_POS() << " parts of speech among " << n_types() << " unique types ("
     << n_ambiguous() << " ambiguous).\n";
  os << " POS Counts among input tokens are : \n";
  for(auto it=mPOSCountMap.cbegin(); it != mPOSCountMap.cend(); ++it)
    std::cout << "       " << it->first << " " << it->second << std::endl;
}
