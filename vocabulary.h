#ifndef _VOCABULARY_H_
#define _VOCABULARY_H_

/*
  Notes:
          * parses numeric values as special tokens (parse_numeric_string)
	  * line oriented input, so optionally can mark line breaks (marked with EOL)
	  * optionally can skip initial tokens in a line (as in regression application)
*/

#include "base_classes.h"

#include <iostream>
#include <list>
#include <vector>
#include <map>

#include <Eigen/Core>
#include <Eigen/SparseCore>

using std::string;

class Vocabulary
{
 public:
  typedef std::vector<Type>                          TypeVector;
  typedef std::list<Type>                            TypeList;
  typedef std::map<Type,int>                         TypeMap;
  typedef std::map<std::pair<int,int>,int>           BigramMap;            // (row,col) positions of frequencies
  typedef Eigen::VectorXf                            Vector;
  typedef Eigen::VectorXi                            IntVector;
  typedef Eigen::SparseMatrix<float,Eigen::RowMajor> SparseMatrix;

  static Type OOV;
  static Type EOL;
  
 private:
  int  const      mSkipInitial;    // skip initial tokens in a line
  bool const      mMarkEOL;        // add end of line token
  int             mNTokens;
  int  const      mMinFrequency;
  int             mOOVIndex;
  TypeList        mTokens;
  TypeVector      mTypeVector;     // types and counts of types, collected into vectors (OOV collected), in freq order
  Vector          mTypeFreqVector; 
  TypeMap         mFreqMap;        // count of non-oov (OOV collected, sorted by types order)
  TypeMap         mOOVMap;         // count of oov types
  TypeMap         mIndexMap;       // position of types (oov mapped to OOV)

 public:

  Vocabulary ()
    : mSkipInitial(0), mMarkEOL(false), mNTokens(0), mMinFrequency(0), mOOVIndex(0), mTokens(),
      mTypeVector(), mTypeFreqVector(), mFreqMap(), mOOVMap(), mIndexMap() { };

 Vocabulary(std::string fileName, int skipInitial, bool markEOL, int minFrequency)
   : mSkipInitial(skipInitial), mMarkEOL(markEOL), mNTokens(0), mMinFrequency(minFrequency), mOOVIndex(0),
     mTokens(), mTypeVector(), mTypeFreqVector(), mFreqMap(), mOOVMap(), mIndexMap() { init_from_file(fileName); }

 Vocabulary(std::istream &is, int skipInitial, bool markEOL, int minFrequency)
   : mSkipInitial(skipInitial), mMarkEOL(markEOL), mNTokens(0), mMinFrequency(minFrequency), mOOVIndex(0),
     mTokens(), mTypeVector(), mTypeFreqVector(), mFreqMap(), mOOVMap(), mIndexMap() { init_from_stream(is); }

  int        n_types()                          const { return mFreqMap.size(); }
  int        type_index(Type const& type)       const;                                // OOV position if not found
  TypeVector types ()                           const { return mTypeVector; }         // types in Frequency order
  Vector     type_frequency_vector()            const { return mTypeFreqVector; }
  TypeMap    oov_map()                          const { return mOOVMap; }
  
  int        n_tokens()                         const { return mNTokens; }
  
  void       fill_sparse_bigram_matrix (SparseMatrix &B, int skip) const;
  void       fill_sparse_regr_design_from_stream (Vector &Y, SparseMatrix &X, Vector &nTokens, std::istream &is, bool normalize=false) const;
  
  void       print_to_stream    (std::ostream & os) const;
  void       write_type_freq    (std::ostream & os) const;
  void       write_oov_to_stream(std::ostream &os, int maxNumberToWrite) const;
  
 private:
  void init_from_file  (std::string fileName);
  void init_from_stream(std::istream &is);
  void parse_line (std::string const& line, std::map<Type,int> &vocab);
  void fill_bigram_map (BigramMap &bm, int skip) const;
};

inline
std::ostream&
operator<<(std::ostream& os, Vocabulary const& vocab)
{
  vocab.print_to_stream(os);
  return os;
}

#endif
