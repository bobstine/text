#include "token_manager.h"

#include <algorithm>
#include <set>
#include <iostream>
#include <fstream>
#include <ios>


std::string outputTag = "TKMG: ";

int
TokenManager::operator[](string s)      const
{
  std::map<string,int>::const_iterator it = mStrToIntMap.find(s);
  if (it != mStrToIntMap.cend())
    return it->second;
  else
    return -1;
}

int
TokenManager::index_of_POS (string pos)  const
{
  auto it = mPOSIndex.find(pos);
  if (it != mPOSIndex.cend())
    return it->second;
  else
    return -1;
} 

void
TokenManager::init_from_file(std::string &fileName, float posThreshold)
{
  std::ifstream input (fileName.c_str());
  if (! input)
    std::cerr << outputTag << "Cannot open input token file " << fileName << std::endl;
  else
    init_from_stream(input, posThreshold);
}


void
TokenManager::init_from_stream(std::istream &input, float posThreshold)
{
  std::map<string,int>posMap;   // initial pos prior to collapse rare pos to OTH
  while (input)
  { string token;
    input >> token;
    if (token.empty()) break;   // end of file
    ++mTypeFreqMap[token];
    string pos;
    input >> pos;
    ++posMap[pos];
    mTokens.push_back(std::make_pair(token,pos));
  }
  std::clog << outputTag << "Read " << n_types() << " types from input of " << mTokens.size() << " strings." << std::endl;
  std::clog << outputTag << "Identified " << posMap.size() << " distinct POS." << std::endl;
  if(posThreshold > 0.0)                                             // remove rare pos tags
  { std::map<string,bool> reduce;
    int reduceCount = 0;
    for(auto it=posMap.cbegin(); it != posMap.cend(); ++it)
    { float freq = ((float)it->second)/mTokens.size();
      if (freq < posThreshold)
      { reduce[it->first] = true; ++reduceCount; }
      else reduce[it->first]=false;
    }
    if (reduceCount > 0)
    { std::clog << outputTag << "Mapped " << reduceCount << " POS categories to OTH." << std::endl;
      for(auto it=mTokens.begin(); it != mTokens.end(); ++it)
	if (reduce[it->second]) it->second = "OTH";
    }
  }
  for(auto it=mTokens.cbegin(); it != mTokens.cend(); ++it)          // build pos map for each type
  { ++mPOSMap[it->second];
    ++mTypePOSMap[it->first][it->second];
  }
  int posIndex = 0;
  for(auto it=mPOSMap.cbegin(); it != mPOSMap.cend(); ++it)          // define index for each pos
  { mPOSIndex[it->first] = posIndex++;
    mIntToPOSVec.push_back(it->first);
  }
  std::multimap<int,string> sortedTokenMap;                          // sort types by frequency by inserting into multimap
  for (auto it = mTypeFreqMap.begin(); it != mTypeFreqMap.end(); ++it)
    sortedTokenMap.insert( std::make_pair(-it->second,it->first) );  // negate so decreasing
  int tokenID = 0;
  mIntToStrVec.resize(n_types());
  for(auto it = sortedTokenMap.begin(); it != sortedTokenMap.end(); ++it)
  { mStrToIntMap[(*it).second] = tokenID;
    mIntToStrVec[tokenID] = (*it).second;
    ++tokenID;
  }
}

TokenManager::StrVector
TokenManager::type_labels()        const
{
  StrVector sv;

  for(auto it=mTypeFreqMap.cbegin(); it != mTypeFreqMap.cend(); ++it)
    sv.push_back(it->first);
  return sv;
}


TokenManager::StrVector
TokenManager::type_POS_labels()        const
{
  StrVector sv;

  for(auto it=mTypeFreqMap.cbegin(); it != mTypeFreqMap.cend(); ++it)
    sv.push_back(POS_of_type(it->first));
  return sv;
}


int
TokenManager::n_ambiguous ()       const
{
  int nAmbiguous = 0;
  for(auto it=mTypePOSMap.cbegin(); it!=mTypePOSMap.cend(); ++it)
  { string type = it->first;
    int max=0;
    for(auto it2=it->second.cbegin(); it2 != it->second.cend(); ++it2)
      if (it2->second > max) max = it2->second;
    nAmbiguous += mTypeFreqMap.at(type)-max;
  }
  return nAmbiguous;
}


string
TokenManager::POS_of_type (string const& s)      const
{
  const std::map<string,int> posMap = mTypePOSMap.at(s);
  if(posMap.empty())
  { std::cerr << outputTag << "Type " << s << " does not have POS recorded. " << std::endl;
    return "";
  }
  else
  { string pos   = posMap.cbegin()->first;
    int    count = posMap.cbegin()->second;
    for(auto it = ++posMap.cbegin(); it != posMap.cend(); ++it)
      if(it->second > count)
      { count = it->second;
	pos   = it->first;
      }
    return pos;
  }
}


TokenManager::CountVector
TokenManager::type_POS_tags (string const& s, bool sort) const
{
  typedef std::pair<string,int> PSI;
  CountVector v;
  const std::map<string,int> posMap = mTypePOSMap.at(s);
  if(posMap.empty())
    std::cerr << outputTag << "Token " << s << " does not have POS recorded. " << std::endl;
  else
  { for(auto it = posMap.cbegin(); it !=posMap.cend(); ++it)
      v.push_back(*it);
    if(sort)
      std::sort(v.begin(), v.end(), [](PSI const& a, PSI const& b) { return a.second < b.second; });
  }
  return v;
}


int
TokenManager::n_types_oov(TokenManager const& tm) const
{
  int nOOV = 0;
  for (auto it = mTypeFreqMap.cbegin(); it != mTypeFreqMap.cend(); ++it)
    if(tm[it->first] < 0) ++ nOOV;
  return nOOV;
}



void
TokenManager::fill_bigram_map(BigramMap &bm, int skip, TokenManager const& tm, bool transpose) const
{
  std::pair<int,int> ij;                 // beyond those indices for tokens found in tm
  auto itBack = mTokens.cbegin();
  for (int i=0; i<skip+1; ++i) ++itBack;
  int nNotFound = 0;
  if (!transpose)
  { for (auto it=mTokens.cbegin(); itBack != mTokens.cend(); ++it, ++itBack)
    { int i = mStrToIntMap.at(it->first);  // i from this
      int j = tm[itBack->first];           // j from index set used in tm
      if (0 <= j)                          // valid
      { ij.first = i;
	ij.second = j;
	++bm[ij];
      }
      else ++nNotFound;
    }
  }
  else 
  { for (auto it=mTokens.cbegin(); itBack != mTokens.cend(); ++it, ++itBack)
    { int i = mStrToIntMap.at(itBack->first);  // i from tm
      int j = tm[it->first];                   // j from this
      if (0 <= j)                              // valid
      { ij.first = i;
	ij.second = j;
	++bm[ij];
      }
      else ++nNotFound;
    }
  }
  if (nNotFound)
    std::clog << outputTag << " *** WARNING *** " << nNotFound << " tokens not found in bigram lookup." << std::endl;
}
    
void
TokenManager::print_tags(int num) const
{
  for(int k=0; k<num; ++k)
  { std::vector<std::pair<string,int> > v = type_POS_tags(mIntToStrVec[k],true);
    std::clog << "k=" << k << "    Type   '" << mIntToStrVec[k] << "'" << std::endl;
    for (auto it=v.begin(); it != v.end(); ++it)
      std::clog << "     " << it->first << " " << it->second << std::endl;
  }
}

void
TokenManager::write_frequencies_to_file(std::string fileName) const   // write type and counts 
{
  std::map<string, int> posMap (POS_map());
  std::vector<int> posCounts (posMap.size());
  int i = 0;
  for(auto it=posMap.cbegin(); it != posMap.cend(); ++it)
  { posCounts[i++] = it->second; }
  std::ofstream output (fileName);
  output << "Types\n";
  for (auto it=mTypeFreqMap.cbegin(); it != mTypeFreqMap.cend(); ++it)
    output << it->second << ", ";
  output << "\nPOS\n";
  for (auto it=posCounts.cbegin(); it != posCounts.cend(); ++it)
    output << *it << ", ";
  output << std::endl;
  output.close();
}

  
