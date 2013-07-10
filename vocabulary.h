#ifndef _VOCABULARY_H_
#define _VOCABULARY_H_

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
  
 private:
  int             mNTokens;
  int             mMinFrequency;
  int             mOOVIndex;
  TypeList        mTokens;
  TypeVector      mTypeVector;
  TypeMap         mFreqMap;        // count of non-oov (OOV collected)
  TypeMap         mOOVMap;         // count of oov types
  TypeMap         mIndexMap;       // position of types (oov mapped to OOV)

 public:

  Vocabulary (): mNTokens(0), mMinFrequency(0), mOOVIndex(0), mTokens(), mTypeVector(), mFreqMap(), mOOVMap(), mIndexMap() { };

 Vocabulary(std::string fileName, int minFrequency)
   : mNTokens(0), mMinFrequency(minFrequency), mOOVIndex(0), mTokens(), mTypeVector(), mFreqMap(), mOOVMap(), mIndexMap() { init_from_file(fileName); }

 Vocabulary(std::istream &is, int minFrequency)
   : mNTokens(0), mMinFrequency(minFrequency), mOOVIndex(0), mTokens(), mTypeVector(), mFreqMap(), mOOVMap(), mIndexMap() { init_from_stream(is); }

  int        n_types()                          const { return mFreqMap.size(); }
  int        type_index(Type const& type)       const;                                // OOV position if not found
  TypeVector types ()                           const { return mTypeVector; }         // types in Frequency order 
  TypeMap    oov_map()                          const { return mOOVMap; }
  
  int        n_tokens()                         const { return mNTokens; }

  void       fill_sparse_bigram_matrix (SparseMatrix &B, int skip) const;

  void       fill_sparse_regr_design (Vector &Y, SparseMatrix &X, std::istream &is) const;
  void       print_to_stream(std::ostream & os) const;
  
 private:
  void init_from_file  (std::string fileName);
  void init_from_stream(std::istream &is);
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
