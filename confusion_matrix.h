#ifndef _CONFUSION_MATRIX_H_
#define _CONFUSION_MATRIX_H_

#include <Eigen/Core>

#include <vector>
#include <map>
#include <set>
#include <string>
#include <iostream>


using std::string;


class ConfusionMatrix
{
  typedef std::map<std::pair<string,string>,int> Map;
  typedef std::set<string>                       StrSet;
  typedef std::vector<string>                    StrVec;
  typedef Eigen::MatrixXi                        Matrix;

  int      mN;
  int      mMaxLabelLen;
  int      mMaxFreq;
  StrVec   mLabels;
  Matrix   mTable;

public:
  
 ConfusionMatrix()
   : mN(0), mMaxLabelLen(0), mMaxFreq(0), mLabels(), mTable() { }
  
  template <class It1, class It2>
    ConfusionMatrix (It1 trueLabelsBegin, It1 trueLabelsEnd, It2 estLabels)
    : mN(0), mMaxLabelLen(0), mMaxFreq(0), mLabels(), mTable()
    { Map m;
      StrSet s;
      for(It1 it=trueLabelsBegin; it != trueLabelsEnd; ++it, ++estLabels)
      { s.insert(*it);
	s.insert(*estLabels);
	++m[std::make_pair(*it, *estLabels)];
      }
      set_labels(s);
      convert_to_matrix(m);
      set_max_freq();
    }

  int         total_count()                     const { return mN; }

  float       accuracy()                        const;                   // (sum of diagonal)/total
  float       purity()                          const;                   // (sum of max in each col)/total
  
  void        print_to_stream(std::ostream &os) const;

 private:
  void        set_labels(StrSet const& s);
  void        convert_to_matrix(Map const& m);
  void        set_max_freq();    
  
};

#endif

