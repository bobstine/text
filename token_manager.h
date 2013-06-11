#ifndef _TOKEN_MANAGER_H_
#define _TOKEN_MANAGER_H_

#include "iterators.h"

#include <string>
#include <vector>
#include <list>
#include <map>
#include <iostream>


using std::string;

// In POS and Type objects, the index is assigned at the time that the token manager
// inserts the objects into its vectors

class TokenManager;

class POS
{
  friend TokenManager;

  string mStr;
  int    mIndex;
  
 public:
  explicit POS(string s) : mStr(s), mIndex(-1) {}

  int      index()                  const { assert(mIndex>=0); return mIndex; }
  
  string   underlying_string()      const { return mStr; }
  
  bool     operator==(POS const& p) const { return mStr == p.mStr; };
  bool     operator< (POS const& p) const { return mStr < p.mStr; };
 protected:
  void     assign_index(int i)            { mIndex = i; }

};


class Type
{
  friend TokenManager;

  string mStr;
  int    mIndex;
  
 public:
  explicit Type(string s) : mStr(s), mIndex(-1) {}

  int      index ()      const              { assert(mIndex>=0); return mIndex; }
  string   underlying_string_representation() const { return mStr; }
  
  bool     operator==(Type const& p) const { return mStr == p.mStr; };
  bool     operator< (Type const& p) const { return mStr  < p.mStr; };

  protected:
   void     assign_index (int i)             { mIndex = i; }

 };


 inline
 std::ostream&
 operator<<(std::ostream &os, POS const& p) { os <<  p.underlying_string(); return os; }


 inline
 std::ostream&
 operator<<(std::ostream &os, Type const& t) { os <<  t.underlying_string_representation(); return os; }


 class TokenManager
 {
   typedef  std::map<Type, std::map<POS,int>>                POSMap;
   typedef  std::vector<POS>                                 POSVector;
   typedef  std::vector<Type>                                TypeVector;
   typedef  std::map<std::pair<int,int>,int>                 BigramMap;
   typedef  std::list<std::pair<Type,POS>>::const_iterator   Iter;
   typedef  std::vector< std::pair<POS,int> >                POSCountVector;

   std::list<std::pair<Type,POS>> mTokens;           // source data
   std::map<Type, int>            mTypeMap;          // count of Types among input tokens
   std::map<POS,  int>            mPOSMap;           //           POS

   std::vector<Type>              mTypeVec;          // sorted order of the types, inverse of str to int map
   std::vector<POS>               mPOSVec;           // pos labels for integers
   POSMap                         mTypePOSMap;       // identifies all POSs observed for a type

  public:

   TokenManager (TokenManager const& tm)
     : mTokens(tm.mTokens), mTypeMap(tm.mTypeMap), mPOSMap(tm.mPOSMap), mTypeVec(tm.mTypeVec),
     mPOSVec(tm.mPOSVec),mTypePOSMap(tm.mTypePOSMap)
     { std::cerr << "WARNING: Copy Construct token map\n"; }

   TokenManager ()
     : mTokens(), mTypeMap(), mPOSMap(), mTypeVec(), mPOSVec(),mTypePOSMap()
     { std::cerr << "WARNING: Construct empty token map\n"; }

   TokenManager (std::istream &input, float posThreshold = 0.0)
     : mTokens(), mTypeMap(), mPOSMap(), mTypeVec(), mPOSVec(),mTypePOSMap()
     { init_from_stream(input, posThreshold); }

   TokenManager (std::string fileName, float posThreshold = 0.0)
     : mTokens(), mTypeMap(), mPOSMap(), mTypeVec(), mPOSVec(),mTypePOSMap()
     { init_from_file(fileName, posThreshold); }
   
  int            input_length()                      const { return (int) mTokens.size(); }  
  Iter           token_list_begin()                  const { return mTokens.cbegin(); }
  Iter           token_list_end()                    const { return mTokens.cend(); }
  
  int            n_types_oov(TokenManager const& tm) const;
  bool           known_type(std::string type)        const { return (mTypeMap.find(Type(type)) != mTypeMap.end()); }
  int            n_types ()                          const { return (int) mTypeMap.size(); }
  TypeVector     type_vector()                       const;
  
  int            n_POS()                             const { return (int) mPOSMap.size(); }
  POSVector      POS_vector()                        const;   
  
  Type           type_of_index (int i)               const { Type t = mTypeVec[i]; assert(t.index()==i); return t; }
  POS            POS_of_index (int i)                const { POS p = mPOSVec[i];   assert(p.index()==i); return p; }
  
  int            type_freq (Type const& type)        const { return mTypeMap.at(type); }
  int            type_freq (int i)                   const { return mTypeMap.at(mTypeVec[i]); }
  
  int            POS_freq (POS const& pos)           const { return mPOSMap.at(pos); }
  int            POS_freq (int i)                    const { return mPOSMap.at(mPOSVec[i]); }
  int            n_ambiguous ()                      const;                          // number with ambiguous category
  POS            POS_of_type  (Type const& type)  const;                          // most common POS for this type
  POSCountVector type_POS_tags(Type const& type, bool sort=false) const;

  
  std::map<POS,int> POS_map()                     const { return mPOSMap; }

  void        fill_bigram_map(BigramMap &bm, int skip)         const;
  void        fill_bigram_map(BigramMap &bm, int skip, TokenManager const& tm, bool transpose=false) const; // use col index from other TM

  void        print_to_stream(std::ostream &os)                const;
  void        print_tags(int k)                                const;
  void        write_frequencies_to_file (string filename)      const;
  
 private:
  void init_from_file(std::string &fileName, float posThreshold);
  void init_from_stream(std::istream &input, float posThreshold);
};

#endif
