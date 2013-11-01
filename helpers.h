#ifndef _HELPERS_H_
#define _HELPERS_H_

#include "vocabulary.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace Helper
{
  
  typedef Eigen::VectorXf Vector;
  typedef Eigen::MatrixXf Matrix;
  
  
  void
    scan_google_vocabulary_for_oov (Vocabulary const& vocab);
  
  
  void
    write_eigenwords_to_file (string fileName, Matrix const& M, Vocabulary const& vocab);
  
  void
    write_word_counts_to_file(std::string fileName, Vocabulary::SparseMatrix const& W, int nCols, Vocabulary const& vocab);
  
  void
    calculate_sequence_r2 (Eigen::VectorXd const& Y, Eigen::VectorXd tokenCount, std::string xLabel, Eigen::MatrixXf X, std::string file);
  
  void
    calculate_sequence_r2 (Eigen::VectorXd const& Y, Eigen::VectorXd tokenCount, bool reverse, Vocabulary::SparseMatrix const& W,
			   Vocabulary const& vocab, int nToFit, std::string file);
  
  void
    fill_random_projection(Matrix &P, Vocabulary::SparseMatrix const& M, Vector const& wts, int power);
  
  void
    fill_left_SVD(Matrix &P, Vocabulary::SparseMatrix const& M);
}

#endif
  
