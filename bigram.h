#ifndef _BIGRAM_H_
#define _BIGRAM_H_

#include <string>
#include <vector>
#include <list>
#include <map>

class TokenManager
{
  std::vector<std::string>  mIntToStrVec;
  std::map<std::string,int> mStrToIntMap;

 public:

 TokenManager(std::map<std::string,int> const &tokenMap)
   : mIntToStrVec(tokenMap.size()), mStrToIntMap()         { initialize(tokenMap); }
  
  int          operator[](std::string const& s) const { return mStrToIntMap.at(s); }
  std::string  operator[](int i)                const { return mIntToStrVec[i]; }

 private:
  // sort map by values by inserting into multimap
  void initialize(std::map<std::string,int> const& tokenMap)
  {
    std::multimap<int,std::string> sortedTokenMap;
    for (auto it = tokenMap.begin(); it != tokenMap.end(); ++it)
      sortedTokenMap.insert( std::make_pair(-it->second,it->first) );  // negate so decreasing
    int tokenID = 0;
    for(auto it = sortedTokenMap.begin(); it != sortedTokenMap.end(); ++it)
    {  mStrToIntMap[(*it).second] = tokenID;
       mIntToStrVec[tokenID] = (*it).second;
       ++tokenID;
    }
  }
};

#endif
