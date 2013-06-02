#ifndef _CROSS_TAB_H_
#define _CROSS_TAB_H_

#include <Eigen/Core>

#include <vector>
#include <string>

using std::string;

class CrossTab
{
  typedef std::vector<std::string> StrVec;

  int mN;
  StrVec  mColLabels;
  StrVec  mRowLabels;
  Eigen::MatrixXi mTable;

public:
  
 CrossTab(int nRows, int nCols)
   : mN(0), mColLabels(), mRowLabels(), mTable(Eigen::MatrixXi::Zero(nRows,nCols)) { init_labels(); }
  
  template <class ItR, class ItC>
    CrossTab (ItR rBegin, ItR rEnd, ItC cBegin, ItC cEnd) : mN(0), mColLabels(), mRowLabels(), mTable()
    { int nR = 0; int nC = 0;
      for(ItR it=rBegin; it != rEnd; ++it)
      { ++nR;	mRowLabels.push_back(*it); }
      for(ItC it=cBegin; it != cEnd; ++it)
      { ++nC;   mColLabels.push_back(*it); }
      mTable = Eigen::MatrixXi::Zero(nR,nC);
    }

  void            increment(int i, int j);
  void            clear_cells()                           { mTable = Eigen::MatrixXi::Zero(mTable.rows(), mTable.cols()); }
  
  int             element(int i, int j)             const { return mTable(i,j); }
  Eigen::VectorXi row_margins()                     const { return mTable.rowwise().sum(); }
  Eigen::VectorXi col_margins()                     const { return mTable.colwise().sum(); }
  int             total()                           const { return mN; }

  Eigen::VectorXi most_common_col_in_each_row()     const; 
  StrVec          most_common_label_in_each_row()   const;

  float           accuracy()                        const { return 100.0*((float)sum_row_max())/mN; }
  int             sum_row_max()                     const { return mTable.rowwise().maxCoeff().sum(); }

  void            print_accuracy_to_stream (std::ostream &os) const;
  void            print_to_stream(std::ostream &os) const;

 private:
  void            init_labels();
  
};

#endif

