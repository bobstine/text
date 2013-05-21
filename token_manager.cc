#include "token_manager.h"

#include <algorithm>
#include <iostream>


void
TokenManager::init_from_stream(std::istream &input, float posThreshold)
{
  std::map <string, int> posMap;
  mTokenLength = 0;
  while (input)
  { string token;
    input >> token;
    if (token.empty()) break;   // end of file
    ++mTokenFreqMap[token];
    string pos;
    input >> pos;
    ++posMap[pos];
    mTokens.push_back(std::make_pair(token,pos));
    ++mTokenLength;
  }
  std::clog << "TKMN: Read " << n_unique_tokens() << " unique tokens from input of " << mTokenLength << " strings." << std::endl;
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
  for(auto it=mTokens.cbegin(); it != mTokens.cend(); ++it)          // build pos map for each token
  { ++mPOSMap[it->second];
    ++mTokenPOSMap[it->first][it->second];
  }
  std::multimap<int,string> sortedTokenMap;                          // sort by values by inserting into multimap
  for (auto it = mTokenFreqMap.begin(); it != mTokenFreqMap.end(); ++it)
    sortedTokenMap.insert( std::make_pair(-it->second,it->first) );  // negate so decreasing
  int tokenID = 0;
  mIntToStrVec.resize(n_unique_tokens());
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
  for(auto it=mTokenPOSMap.cbegin(); it!=mTokenPOSMap.cend(); ++it)
  { string token = it->first;
    int max=0;
    for(auto it2=it->second.cbegin(); it2 != it->second.cend(); ++it2)
      if (it2->second > max) max = it2->second;
    nAmbiguous += mTokenFreqMap.at(token)-max;
  }
  return nAmbiguous;
}


string
TokenManager::token_POS (string const& s)      const
{
  const std::map<string,int> posMap = mTokenPOSMap.at(s);
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
TokenManager::token_POS_tags (string const& s, bool sort) const
{
  typedef std::pair<string,int> PSI;
  CountVector v;
  const std::map<string,int> posMap = mTokenPOSMap.at(s);
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


void
TokenManager::fill_bigram_map (BigramMap &bm) const
{
  std::pair<int,int> i = std::make_pair(0,mStrToIntMap.at(mTokens.front().first));
  for(auto it = ++mTokens.cbegin(); it!=mTokens.cend(); ++it)
  { i.first = i.second;
    i.second= mStrToIntMap.at(it->first);
    ++bm[i];
  }
}


void
TokenManager::print_tags(int num) const
{
  for(int k=0; k<num; ++k)
  { std::vector<std::pair<string,int> > v = token_POS_tags(mIntToStrVec[k],true);
    std::clog << "k=" << k << "    Token   '" << mIntToStrVec[k] << "'" << std::endl;
    for (auto it=v.begin(); it != v.end(); ++it)
      std::clog << "     " << it->first << " " << it->second << std::endl;
  }
}
