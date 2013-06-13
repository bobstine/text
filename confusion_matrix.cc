#include "confusion_matrix.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <assert.h>


static std::string messageTag ("CFMX: ");


inline
int
max(int a, int b) { return (a < b) ? b : a; }


//  accuracy     accuracy     accuracy     accuracy     accuracy     accuracy     accuracy

float
ConfusionMatrix::accuracy() const         // diagonal sum divided by total sum
{
  return ((float) mTable.diagonal().sum())/mN;
}


float
ConfusionMatrix::purity() const
{
  return ((float) mTable.colwise().maxCoeff().sum())/mN;
}


//  print     print     print     print     print     print     print     print     print     print


void
ConfusionMatrix::print_to_stream (std::ostream &os) const
{
  using std::setw;
  os << messageTag << " Accuracy=" << accuracy() << "  Purity=" << purity() << "      (rows=truth, cols=estimated)" << std::endl;
  int width = 2 + max(log10((float)mMaxFreq), mMaxLabelLen);
  
  os << std::setw(width) << messageTag << "  " << setw(width) << " ";
  for (auto it = mLabels.cbegin(); it != mLabels.cend(); ++it)
    os << setw(width) << *it;
  os << std::endl << "  " << setw(width) << " " << setw(width) << " ";
  Eigen::VectorXi colMargins = mTable.colwise().sum();
  for(int col=0; col<mTable.cols(); ++col)
    os << setw(width) << colMargins(col);
  os << std::endl << "        ------------------------------------------" << std::endl;
  Eigen::VectorXi rowMargins = mTable.rowwise().sum();
  for(int row=0; row<mTable.rows(); ++row)
  { os << setw(width) << mLabels[row] << setw(width) << rowMargins[row] << " | ";
    for (int c=0; c<mTable.cols(); ++c)
      os << setw(width) << mTable(row,c);
    os << std::endl;
  }
}

/////////////  Initialize     Initialize     Initialize     Initialize     Initialize     Initialize     Initialize     

void
ConfusionMatrix::set_labels(ConfusionMatrix::StrSet const& s)
{
  for(auto it=s.cbegin(); it!=s.cend(); ++it)
  { int len (it->size());
    if(len > mMaxLabelLen)
      mMaxLabelLen = len;
    mLabels.push_back(*it);
  }
}


void
ConfusionMatrix::convert_to_matrix(ConfusionMatrix::Map const& m)
{
  mTable = Matrix::Zero(mLabels.size(), mLabels.size());
  int row=0;
  for (auto it1=mLabels.cbegin(); it1 != mLabels.cend(); ++it1,++row)
  { int col=0;
    for (auto it2=mLabels.cbegin(); it2 != mLabels.cend(); ++it2,++col)
    { std::pair<string, string> cell = std::make_pair(*it1, *it2);
      if(m.count(cell)>0)
	mTable(row,col) = m.at(cell);
    }
  }
}


void
ConfusionMatrix::set_max_freq()
{
  mMaxFreq = mTable.maxCoeff();
  mN = mTable.sum();
  std::clog << messageTag << " Defined confusion matrix with " << mLabels.size()
	    << " rows & cols, with total count " << mN << "." << std::endl;
}
  
