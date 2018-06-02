// C++ incQRdes
#include <iostream>
#include <cmath>
#include <climits>
#include <vector>
#include <queue>
#include <algorithm>

/** @brief Defines the QR decomposition algorithm of a square matrix.
 *
 *	   The matrices Q and R are such that \f$A = QR\f$ with
 *	    - \f$Q\f$ an orthogonal matrix;
 *	    - \f$R\f$ an upper triangular matrix.
 */
template<typename T> class DecompQR {
 public:
  /** @brief Default constructor, initilizes an empty matrix of dimension 0x0 */
  DecompQR();

  /** @brief Std constructor, sets the underlying std matrix */
  DecompQR(const std::vector<std::vector<T>>& mat);

  /** @brief Copy constructor */
  DecompQR(const DecompQR& decQR);

  /** @brief Equality operator */
  DecompQR& operator=(const DecompQR& decQR);

  /** @brief Destructor */
  ~DecompQR();

  /** @brief Decomposes a square matrix into a parity matrix P, an orthogonal
   *	     matrix Q and a right triangular matrix R
   *
   *  @param	mat	Input square matrix
   *
   *  @return		True if successful
   */
  void DecomposeGramSchmidt(const std::vector<std::vector<T>>& mat);

  /** @brief Returns the determinant of the decomposed matrix
   *  @return 		Determinant of (P^{-1}QR)
   */
  T Determinant();

  /** @brief Returns a matrix of m rows and n columns
   *
   *  @param	m	Number of rows
   *  @param	n	Number of coQRmns
   *
   *  @return 		True if successful
   */
  std::vector<std::vector<T>> NewMatrix(size_t m, size_t n);

  /** @brief Returns the dot produt of two vectors
   *
   *  @param	a	First vector
   *  @param	b	Second vector
   *
   *  @return 		Dot product of the two vectors
   */
  T DotProduct(const std::vector<T>& a, const std::vector<T>& b);

  /** @brief Gets and returns the lower triangular matrix */
  const std::vector<std::vector<T>>& Q()	{ return _Q; }

  /** @brief Gets and returns the upper triangular matrix */
  const std::vector<std::vector<T>>& R()	{ return _R; }

 private:

  std::vector<std::vector<T>>	_Q;	///< Q Matrix
  std::vector<std::vector<T>>	_R;	///< R Matrix
};

#include "DecompQR-inl.hh"
