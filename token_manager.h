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
  std::map<string, int>                  mTypeFreqMap;

  std::vector<string>                    mIntToStrVec;     // sorted order of the tokens, inverse of str to int map
  std::map<string,int>                   mStrToIntMap;     // assigns integer to a token based on freq
  std::map<string,int>                   mPOSMap;          // count of POS in input tokens
  POSMap                                 mTypePOSMap;      // identifies POS's for a token

 public:
  
  TokenManager (TokenManager const& tm)
    : mTokens(tm.mTokens), mTokenLength(tm.mTokenLength), mTypeFreqMap(tm.mTypeFreqMap),
      mIntToStrVec(tm.mIntToStrVec), mStrToIntMap(tm.mStrToIntMap), mPOSMap(tm.mPOSMap), mTypePOSMap(tm.mTypePOSMap)
    { std::cerr << "WARNING: Copy Construct token map\n"; }
  
  TokenManager (std::istream &input, float posThreshold = 0.0)
    : mTokens(), mTokenLength(0), mTypeFreqMap(), mIntToStrVec(), mStrToIntMap(), mPOSMap(), mTypePOSMap()
    { init_from_stream(input, posThreshold); }

  TokenManager (std::string fileName, float posThreshold = 0.0)
    : mTokens(), mTokenLength(0), mTypeFreqMap(), mIntToStrVec(), mStrToIntMap(), mPOSMap(), mTypePOSMap()
    { init_from_file(fileName, posThreshold); }

  int         input_length()                      const { return mTokenLength; }
  int         n_types ()                          const { return (int) mTypeFreqMap.size(); }
  int         n_types_oov(TokenManager const& tm) const;
  Iter        token_list_begin()                  const { return mTokens.cbegin(); }
  Iter        token_list_end()                    const { return mTokens.cend(); }
  
  int         operator[](string const& s)         const;                          // returns -1 if not found
  string      operator[](int i)                   const { return mIntToStrVec[i]; }

  int         type_freq (string const& s)         const { return mTypeFreqMap.at(s); }
  int         n_ambiguous ()                      const;                          // number with ambiguous category
  string      type_POS (string const& s)          const;                          // most common POS for this type
  CountVector type_POS_tags (string const& s, bool sort=false) const;

  std::map<string,int> POS_map()                               const { return mPOSMap; }

  int         fill_bigram_map(BigramMap &bm, int skip)         const;
  int         fill_bigram_map(BigramMap &bm, int skip, TokenManager const& tm) const; // uses ints from other TM
			      
  void        print_tags(int k)                const;
  
 private:
  void init_from_file(std::string &fileName, float posThreshold);
  void init_from_stream(std::istream &input, float posThreshold);
};

#endif
