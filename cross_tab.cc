#include "cross_tab.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <assert.h>

std::string tag ("CRTB: ");


void
CrossTab::increment(int i, int j)
{
  if ((i<0) || (i>=mTable.rows()))
    std::cerr << tag << "Row " << i << " is out of bounds." << std::endl;
  else if ((j<0) || (j>=mTable.cols()))
    std::cerr << tag << "Col " << j << " is out of bounds." << std::endl;
  else
  { ++mN; ++mTable(i,j); }
}


Eigen::VectorXi
CrossTab::most_common_col_in_each_row()     const
{
  Eigen::VectorXi cols (mTable.rows());

  for(int i=0; i<mTable.rows(); ++i)
  { int max=0;
    int maxIn = 0;
    for(int j=0; j<mTable.cols(); ++j)
    { if (mTable(i,j)>max)
      { max = mTable(i,j);
	maxIn=j;
      }
    }
    cols[i] = maxIn;
  }
  return cols;
}


CrossTab::StrVec
CrossTab::most_common_label_in_each_row() const
{
  Eigen::VectorXi cols = most_common_col_in_each_row();
  std::vector<string> labels;
  for (int i=0; i<cols.size(); ++i)
    labels.push_back(mColLabels[i]);
  return labels;
}

void
CrossTab::print_accuracy_to_stream (std::ostream &os) const
{ 
  os << tag << "Classify " << sum_row_max() << " correctly out of " << mN << " (" << accuracy() << "%)\n";
}


void
CrossTab::print_to_stream (std::ostream &os) const
{ 
  int width = 2+log10(mTable.array().maxCoeff());
  os << std::setw(width) << "CrT " << std::setw(width) << " ";
  for (auto it = mColLabels.cbegin(); it != mColLabels.cend(); ++it)
    os << std::setw(width) << *it;
  os << std::endl << std::setw(width) << " " << std::setw(width) << " ";
  Eigen::VectorXi colMargins = col_margins();
  for(int col=0; col<mTable.cols(); ++col)
    os << std::setw(width) << colMargins(col);
  os << std::endl;
  Eigen::VectorXi rowMargins = row_margins();
  for(int row=0; row<mTable.rows(); ++row)
    os << std::setw(width) << mRowLabels[row] << std::setw(width) << rowMargins[row] << " " << mTable.row(row) << std::endl;
}


void
CrossTab::init_labels ()
{
  std::ostringstream ss;

  for(int i=0; i<mTable.rows(); ++i)
  { ss << "R" << i << " ";
    mRowLabels.push_back(ss.str());
    ss.str("");
  }
  for(int i=0; i<mTable.cols(); ++i)
  { ss << "C" << i << " ";
    mColLabels.push_back(ss.str());
    ss.str("");
  }
}
