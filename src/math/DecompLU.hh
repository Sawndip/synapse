// C++ includes
#include <iostream>
#include <cmath>
#include <climits>
#include <vector>
#include <queue>
#include <algorithm>

// Other includes
#include "Assert.hh"

/** @brief Defines the LU decomposition algorithm of a square matrix.
 *
 *	   The matrices L and U are such that \f$PA = LU\f$ or \f$A = P^{-1}LU\f$ with
 *	    - \f$P\f$ is the permutation applied to the matrix to decompose;
 *	    - \f$L\f$ is a lower triangular matrix of diagonal elements unity;
 *	    - \f$U\f$ is an upper triangular matrix.
 */
template<typename T> class DecompLU {
 public:
  /** @brief Default constructor, initilizes an empty matrix of dimension 0x0 */
  DecompLU();

  /** @brief Std constructor, sets the underlying std matrix */
  DecompLU(std::vector<std::vector<T>> mat);

  /** @brief Copy constructor */
  DecompLU(const DecompLU& decLU);

  /** @brief Equality operator */
  DecompLU& operator=(const DecompLU& decLU);

  /** @brief Destructor */
  ~DecompLU();

  /** @brief Decomposes a square matrix into a parity matrix P, a lower
   *	     triangular matrix L and an upper triangular matrix U
   *
   *  @param	mat	Input square matrix
   *
   *  @return		True if successful
   */
  void DecomposeCrout(std::vector<std::vector<T>> mat);

  /** @brief Returns the determinant of the decomposed matrix
   *  @return 		Determinant of (P*L*U)
   */
  T Determinant();

  /** @brief Returns a matrix of m rows and n columns
   *
   *  @param	m	Number of rows
   *  @param	n	Number of columns
   *
   *  @return 		True if successful
   */
  std::vector<std::vector<T>> NewMatrix(size_t m, size_t n);

  /** @brief Gets and returns the parity of the permutations */
  const T& Factor()				{ return _fact; }

  /** @brief Gets and returns whether the LU decomposition was successful or not */
  const bool& Complete()			{ return _complete; }

  /** @brief Gets and returns the permutation matrix */
  const std::vector<std::vector<T>>& P()	{ return _P; }

  /** @brief Gets and returns the permutation vector, contains the ids of the rows */
  const std::vector<T>& Pivots()		{ return _pivots; }

  /** @brief Gets and returns the lower triangular matrix */
  const std::vector<std::vector<T>>& L()	{ return _L; }

  /** @brief Gets and returns the upper triangular matrix */
  const std::vector<std::vector<T>>& U()	{ return _U; }

 private:

  T				_fact;		///< Parity of the permutation (1 or -1)
  bool				_complete;	///< True if the decomposition was successful
  std::vector<std::vector<T>>	_P;		///< Permutation matrix
  std::vector<T>		_pivots;	///< Premutation vector (IDs of the rows)
  std::vector<std::vector<T>>	_L;		///< L Matrix
  std::vector<std::vector<T>>	_U;		///< U Matrix
};

#include "DecompLU-inl.hh"
