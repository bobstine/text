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


class POS: public string
{  
 public:
  
  explicit POS()         : string("") { }
  explicit POS(string s) : string(s) { }
  explicit POS(char *s)  : string(s) { }
  
  void  print_to_stream(std::ostream &os) const   { os << "|" << static_cast<string>(*this) << "|"; }
};


class Type: public string
{
 public:
  explicit Type()         : string("") { }
  explicit Type(string s) : string(s)       { }
  explicit Type(char *s)  : string(s)       { }

  void  print_to_stream(std::ostream &os) const   { os << "`" << static_cast<string>(*this) << "'" ; }
 };


 inline
 std::ostream&
 operator<<(std::ostream &os, POS const& p) { p.print_to_stream(os); return os; }


 inline
 std::ostream&
 operator<<(std::ostream &os, Type const& t) { t.print_to_stream(os); return os; }


 class TokenManager
 {
   typedef  std::map<Type, std::map<POS,int>>                POSMap;
   typedef  std::vector<POS>                                 POSVector;
   typedef  std::vector<Type>                                TypeVector;
   typedef  std::map<std::pair<int,int>,int>                 BigramMap;
   typedef  std::list<std::pair<Type,POS>>::const_iterator   Iter;
   typedef  std::vector< std::pair<POS,int> >                POSCountVector;

   std::list<std::pair<Type,POS>> mTokens;           // source data
   std::map<Type, int>            mTypeCountMap;     // count of Types among input tokens
   std::map<Type, int>            mTypeIndexMap;     // map from type to index in following type vector
   TypeVector                     mTypeVector;       // holds all types, in order of decreasing frequency
   std::map<POS,  int>            mPOSCountMap;      //           POS
   POSMap                         mTypePOSMap;       // identifies all POSs observed for a type

  public:

   TokenManager (TokenManager const& tm)
     : mTokens(tm.mTokens), mTypeCountMap(tm.mTypeCountMap), mTypeIndexMap(tm.mTypeIndexMap), mTypeVector(tm.mTypeVector),
     mPOSCountMap(tm.mPOSCountMap), mTypePOSMap(tm.mTypePOSMap)
     { std::cerr << "WARNING: Copy Construct token map\n"; }

   TokenManager ()
     : mTokens(), mTypeCountMap(), mTypeIndexMap(), mTypeVector(), mPOSCountMap(), mTypePOSMap()
     { std::cerr << "WARNING: Construct empty token map\n"; }

   TokenManager (std::istream &input, float posThreshold = 0.0)
     : mTokens(), mTypeCountMap(), mTypeIndexMap(), mTypeVector(), mPOSCountMap(), mTypePOSMap()
     { init_from_stream(input, posThreshold); }

   TokenManager (std::string fileName, float posThreshold = 0.0)
     : mTokens(), mTypeCountMap(), mTypeIndexMap(), mTypeVector(), mPOSCountMap(), mTypePOSMap()
     { init_from_file(fileName, posThreshold); }
   
  int            input_length()                      const { return (int) mTokens.size(); }  
  Iter           token_list_begin()                  const { return mTokens.cbegin(); }
  Iter           token_list_end()                    const { return mTokens.cend(); }

  int            n_ambiguous ()                      const;                             // number types with ambiguous POS
  bool           known_type(std::string type)        const { return (mTypeCountMap.count(Type(type)) > 0); }
  
  int            n_types ()                          const; 
 
  int            n_types_oov(TokenManager const& tm) const;
  TypeVector     type_vector()                       const { return mTypeVector; }     // assert type_index(typeVector[i]) == i
  int            type_index(Type const& type)        const;
  int            type_freq (Type const& type)        const;
  
  int            n_POS()                             const { return (int) mPOSCountMap.size(); }
  POSVector      POS_vector()                        const;   
  int            POS_freq (POS const& pos)           const;
  
  POS            POS_of_type     (Type const& type)  const;                             // most common POS for this type
  POSCountVector POS_tags_of_type(Type const& type, bool sort=false) const;
                                                   
  void           fill_bigram_map(BigramMap &bm, int skip, TokenManager const& tm, bool transpose=false) const; // skip=0 for standard bigram

  void           print_to_stream(std::ostream &os)                const;
  void           print_type_tags(int maxToPrint)                  const;
  void           write_frequencies_to_file (string filename)      const;
  
 private:
  void init_from_file(std::string &fileName, float posThreshold);
  void init_from_stream(std::istream &input, float posThreshold);
};

#endif
