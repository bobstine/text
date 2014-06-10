#ifndef _HELPERS_H_
#define _HELPERS_H_

#include "vocabulary.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace Helper
{
  typedef float ScalarType;
  typedef Eigen::VectorXf Vector;
  typedef Eigen::MatrixXf Matrix;
  
  
  void
    scan_google_vocabulary_for_oov (Vocabulary const& vocab);

  float
    entropy(Vector const& prob);                     // assumes sum p(i) = 1 and p(i) >= 0
  float
    entropy(Vector const& x, float xSum);            // assumes x(i) >= 0

  void
    sort_columns_using_tfidf (Vocabulary::SparseMatrix *W, Vocabulary const& vocabulary);

  void
    scale_doc_term_matrix(char adjust, Vocabulary::SparseMatrix *W, Vocabulary const& vocabulary, Vector const& nTokens);
    
  
  Vector
    document_frequency_vector (Vocabulary::SparseMatrix const& W);    // count number of docs in which terms appear

  
  void
    write_eigenwords_to_file (string fileName, Matrix const& M, Vocabulary const& vocab);
  
  void
    write_word_counts_to_file(std::string fileName, Vocabulary::SparseMatrix const& W, int nCols, Vocabulary const& vocab);
  
  template <class Mat>
    void
    calculate_sequence_r2 (Eigen::VectorXd const& Y, Eigen::VectorXd tokenCount, int degree, bool reverse, Mat const& W,
				   Vocabulary const& vocab, int nToFit, string file);


  void
    fill_random_projection_svd(Matrix* U, Vector* d, Matrix *V, Vocabulary::SparseMatrix const& B, int powerIterations);
  
  // old version
  void
    fill_random_projection(Matrix &P, Vocabulary::SparseMatrix const& B, Vector const& leftWts, Vector const& rightWts, int power);

  
  void
    write_exact_svd_to_path(Vocabulary::SparseMatrix const& B, int nProjections, std::string path, std::string tag);

  void
    write_matrix_to_file(Matrix const& A, std::string fileName, std::string columnPrefix);
  
  void
    fill_left_SVD(Matrix &P, Vocabulary::SparseMatrix const& M);
}


  
#endif
  
