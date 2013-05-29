#include "token_manager.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <ios>


int
TokenManager::operator[](string const& s)      const
{
  std::map<string,int>::const_iterator it = mStrToIntMap.find(s);
  if (it != mStrToIntMap.cend())
    return it->second;
  else
    return -1;
}



void
TokenManager::init_from_file(std::string &fileName, float posThreshold)
{
  std::ifstream input (fileName.c_str());
  if (! input)
    std::cerr << "TKMG: Cannot open input token file " << fileName << std::endl;
  else
    init_from_stream(input, posThreshold);
}


void
TokenManager::init_from_stream(std::istream &input, float posThreshold)
{
  std::map <string, int> posMap;
  mTokenLength = 0;
  while (input)
  { string token;
    input >> token;
    if (token.empty()) break;   // end of file
    ++mTypeFreqMap[token];
    string pos;
    input >> pos;
    ++posMap[pos];
    mTokens.push_back(std::make_pair(token,pos));
    ++mTokenLength;
  }
  std::clog << "TKMN: Read " << n_types() << " types from input of " << mTokenLength << " strings." << std::endl;
  std::clog << "TKMN: Identified " << posMap.size() << " distinct POS." << std::endl;
  if(posThreshold > 0.0)                                             // remove rare pos tags
  { std::map<string,bool> reduce;
    int reduceCount = 0;
    for(auto it=posMap.cbegin(); it != posMap.cend(); ++it)
    { float freq = ((float)it->second)/mTokenLength;
      if (freq < posThreshold)
      { reduce[it->first] = true; ++reduceCount; }
      else reduce[it->first]=false;
    }
    if (reduceCount > 0)
    { std::clog << "TKMN: Mapped " << reduceCount << " POS categories to OTH." << std::endl;
      for(auto it=mTokens.begin(); it != mTokens.end(); ++it)
	if (reduce[it->second]) it->second = "OTH";
    }
  }
  for(auto it=mTokens.cbegin(); it != mTokens.cend(); ++it)          // build pos map for each type
  { ++mPOSMap[it->second];
    ++mTypePOSMap[it->first][it->second];
  }
  std::multimap<int,string> sortedTokenMap;                          // sort by values by inserting into multimap
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
TokenManager::type_POS (string const& s)      const
{
  const std::map<string,int> posMap = mTypePOSMap.at(s);
  if(posMap.empty())
  { std::cerr << "TKMN: Token " << s << " does not have POS recorded. " << std::endl;
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
    std::cerr << "TKMN: Token " << s << " does not have POS recorded. " << std::endl;
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

int
TokenManager::fill_bigram_map(BigramMap &bm, int skip, TokenManager const& tm) const
{
  int oovIndex = tm.n_types();
  std::pair<int,int> ij;
  if (0 == skip)
  { int j = tm[mTokens.cbegin()->first];
    ij = std::make_pair(0, (j < 0) ? oovIndex++ : j);
    for(auto it=++mTokens.cbegin(); it != mTokens.cend(); ++it)
    { ij.first = ij.second;
      j = tm[it->first];
      ij.second = (j < 0) ? oovIndex++ : j;
      ++bm[ij];
    }
  }
  else
    std::cerr << "TKMG: *** ERROR *** The function fill_bigram_map not defined with skip > 0\n" ;
  return oovIndex;  // max 'rows' of map
}

    
int
TokenManager::fill_bigram_map (BigramMap &bm, int skip) const
{
  std::pair<int,int> i;
  if (0 == skip)
  { i = std::make_pair(0, mStrToIntMap.at(mTokens.cbegin()->first));
    for(auto it = ++mTokens.cbegin(); it!=mTokens.cend(); ++it)
    { i.first = i.second;
      i.second= mStrToIntMap.at(it->first);
      ++bm[i];
    }
  }
  else // track two pointers over the list
  { auto itBack = mTokens.cbegin();
    for (int i=0; i<skip+1; ++i) ++itBack;
    for(auto it = mTokens.cbegin(); itBack!=mTokens.cend(); ++it, ++itBack)
    { i.first = mStrToIntMap.at(it    ->first);
      i.second= mStrToIntMap.at(itBack->first);
      ++bm[i];
    }
  }
  return n_types();
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
