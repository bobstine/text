#ifndef _TOKEN_MANAGER_H_
#define _TOKEN_MANAGER_H_

#include <string>
#include <vector>
#include <list>
#include <map>
#include <iostream>


using std::string;

class TokenManager
{
  typedef  std::map<string, std::map<string,int>>              POSMap;
  typedef  std::map<std::pair<int,int>,int>                    BigramMap;
  typedef  std::list<std::pair<string,string>>::const_iterator Iter;
  typedef  std::vector<std::pair<string,int>>                  CountVector;
  
  std::list<std::pair<string,string>>    mTokens;          // source data
  int                                    mTokenLength;
  std::map<string, int>                  mTokenFreqMap;    // frequency of tokens

  std::vector<string>                    mIntToStrVec;     // sorted order of the tokens, inverse of str to int map
  std::map<string,int>                   mStrToIntMap;     // assigns integer to a token based on freq
  std::map<string,int>                   mPOSMap;          // count of POS in input tokens
  POSMap                                 mTokenPOSMap;     // identifies POS's for a token

 public:
  
  TokenManager (TokenManager const& tm)
    : mTokens(tm.mTokens), mTokenLength(tm.mTokenLength), mTokenFreqMap(tm.mTokenFreqMap),
      mIntToStrVec(tm.mIntToStrVec), mStrToIntMap(tm.mStrToIntMap), mPOSMap(tm.mPOSMap), mTokenPOSMap(tm.mTokenPOSMap)
    { std::cerr << "WARNING: Copy Construct token map\n"; }
  
  TokenManager (std::istream &input, float posThreshold = 0.0)
    : mTokens(), mTokenLength(0), mTokenFreqMap(), mIntToStrVec(), mStrToIntMap(), mPOSMap(), mTokenPOSMap()
    { init_from_stream(input, posThreshold); }

  int         input_length()                   const { return mTokenLength; }
  int         n_unique_tokens ()               const { return (int) mTokenFreqMap.size(); }
  Iter        token_list_begin()               const { return mTokens.cbegin(); }
  Iter        token_list_end()                 const { return mTokens.cend(); }
  
  int         operator[](string const& s)      const { return mStrToIntMap.at(s); }
  string      operator[](int i)                const { return mIntToStrVec[i]; }

  int         token_freq (string const& s)     const { return mTokenFreqMap.at(s); }
  int         n_ambiguous ()                   const;      // number tokens with ambiguous category
  string      token_POS (string const& s)      const;      // most common POS for this token
  CountVector token_POS_tags (string const& s, bool sort=false) const;

  std::map<string,int> POS_map()               const { return mPOSMap; }

  void        fill_bigram_map(BigramMap &bm)   const;
  
  void        print_tags(int k)                const;
  
 private:
  void init_from_stream(std::istream &input, float posThreshold);
};

#endif
