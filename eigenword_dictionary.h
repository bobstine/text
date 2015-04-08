#ifndef _EIGENWORD_DICTIONARY_H_
#define _EIGENWORD_DICTIONARY_H_

#include "text_base_classes.h"

#include <string>
#include <vector>
#include <map>
#include <iostream>

#include <Eigen/Core>

using std::string;

class EigenwordDictionary
{
 public:
  typedef Eigen::VectorXf         Vector;
  typedef Eigen::MatrixXf         Matrix;
  typedef std::map<Text::Type,int>TypeMap;
  typedef std::vector<Text::Type> TypeVector;

 private:
  TypeVector mTypes;
  TypeMap    mIndexMap;
  Matrix     mEigenwords;

 public:
  EigenwordDictionary (string fileName, int nTypes, int dim)                         // file needs header
    : mTypes(), mIndexMap(), mEigenwords(nTypes,dim) { init_from_file(fileName); }

  int        n_types()                           const  { return (int)mEigenwords.rows(); }
  int        eigen_dim()                         const  { return (int)mEigenwords.cols(); }
 
  TypeVector types()                             const  { return mTypes; }
  Matrix     eigenwords()                        const  { return mEigenwords; }
  Matrix     eigenwords(TypeVector const& types) const;

  int        type_index(Text::Type const& t)     const; // -1 if not found
  Vector     operator[](int row)                 const  { return mEigenwords.row(row); }
  Vector     operator[](Text::Type const& t)     const;

  void       print_to_stream(std::ostream &os)   const;
  
 private:
  int init_from_file(string fileName); 

};

inline
std::ostream&
operator<<(std::ostream &os, EigenwordDictionary const& ed)
{
  ed.print_to_stream(os);
  return os;
}

#endif
