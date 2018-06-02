#ifndef MATRIX_HH
#define MATRIX_HH

// C++ includes
#include <iostream>
#include <cmath>
#include <vector>

// Other includes
#include "DecompLU.hh"
#include "DecompQR.hh"
#include "Vector.hh"
#include "Pitch.hh"
#include "Assert.hh"

/** @brief Object that can host any type of data in an matrix arrangement. 
 *
 *    	   The class definition also contains a multitude of linear algrebra
 *         methods that are useful to work with matrices.
 */
template<typename T> class Matrix {
 public:
  /** @brief Default constructor, initilizes an empty matrix of dimension 0x0 */
  Matrix();

  /** @brief Square constructor, initializes an empty nxn square matrix */
  Matrix(const size_t n);

  /** @brief Size constructor, initializes an empty mxn matrix */
  Matrix(const size_t m, const size_t n);

  /** @brief Constant constructor, initializes an mxn matrix filled with value */
  Matrix(const size_t m, const size_t n, const T& value);

  /** @brief Column vector constructor, initializes an nx1 matrix filled with v */
  Matrix(const std::vector<T>& v);

  /** @brief Std constructor, sets the underlying std matrix */
  Matrix(const std::vector<std::vector<T>>& mat);

  /** @brief Copy constructor */
  Matrix(const Matrix& mat);

  /** @brief Destructor */
  ~Matrix();

  /** @brief Equality operator */
  Matrix& operator=(const Matrix& mat);

  /** @brief Boolean operator */
  explicit operator bool() const;

  /** @brief Boolean not operator */
  bool operator!() const;

  /** @brief Overload the unary minus for matrix algebra */
  Matrix operator-() const;

  /** @brief Overload the matrix sum */
  Matrix operator+(const Matrix& mat) const;

  /** @brief Overload the matrix sum */
  void operator+=(const Matrix& mat);

  /** @brief Overload the matrix substraction */
  Matrix operator-(const Matrix& mat) const;

  /** @brief Overload the matrix substraction */
  void operator-=(const Matrix& mat);

  /** @brief Overload the matrix product */
  Matrix operator*(const Matrix& mat) const;

  /** @brief Overload the matrix product equality */
  void operator*=(const Matrix& mat);

  /** @brief Overload the matrix product with a scalar */
  Matrix operator*(const T& factor) const;

  /** @brief Overload the matrix product with a scalar */
  void operator*=(const T& factor);

  /** @brief Overload the matrix product with a vector, returns a vector */
  Vector<T> operator*(const Vector<T>& vec) const;

  /** @brief Overload the matrix ratio with a scalar */
  Matrix operator/(const T& factor) const;

  /** @brief Overload the matrix ratio with a scalar */
  void operator/=(const T& factor);

  /** @brief Overload the subscript operator, return an std::vector row */
  const std::vector<T>& operator[](const size_t i) const;

  /** @brief Overload the subscript operator, return an std::vector row, allows mods */
  std::vector<T>& operator[](const size_t i);

  /** @brief Returns the ith-jth cofactor of the matrix
   *
   *  @param	i	Row index
   *  @param	j	Column index
   *
   *  @return		Cofactor
   */
  T Cofactor(const size_t i, const size_t j) const;

  /** @brief Returns the inverse of the matrix
   *
   *  @return		Inverse matrix
   */
  Matrix<T> CofactorMatrix() const;

  /** @brief Returns the deterimant of the matrix
   *
   *  @return		Value of the determinant
   */
  T Determinant() const;

  /** @brief Returns the i^th column of the matrix as a vector
   *
   *  @return		i^th column as a vector
   */
  std::vector<T> Column(const size_t i) const;

  /** @brief Returns the matrix of which the columns are the eigenvectors of the matrix
   *
   *  @param	lambdas	Eigenvalues
   *  @param	eps	Precision criterion to stop the iteration
   *
   *  @return		Matrix of eigenvectors
   */
  Matrix<T> EigenVectors(std::vector<T>& lambdas, const double eps=1e-6) const;

  /** @brief Fills the matrix with a certain value
   *
   *  @param	value	Value to fill the matrix with
   */
  void Fill(const T& value);

  /** @brief Sets the matrix to an identity matrix */
  void Identity();

  /** @brief Returns the inverse of the matrix
   *
   *  @return		Inverse matrix
   */
  Matrix<T> Inverse() const;

  /** @brief Inverts the matrix in place */
  void Invert();

  /** @brief Checks whether the matrix is symmetric or not
     @return		True is symmetric
   */
  bool IsSymmetric() const;

  /** @brief Returns the number of rows
   *
   *  @return		Number of rows
   */
  const size_t Nrows() const;

  /** @brief Returns the number of columns
   *
   *  @return		Number of columns
   */
  const size_t Ncols() const;

  /** @brief Print the matrix */
  void Print()  const;

  /** @brief Removes the i^th column of the matrix */
  void RemoveColumn(const size_t i);

  /** @brief Removes the i^th row of the matrix */
  void RemoveRow(const size_t i);

  /** @brief Resizes the matrix to chosen dimension
   *
   *  @param	m	Number of rows
   *  @param	n	Number of columns
   */
  void Resize(const size_t m, const size_t n);

  /** @brief Returns the i^th row of the matrix as a vector
   *
   *  @return		i^th row as a vector
   */
  const std::vector<T>& Row(const size_t i) const;

  /** @brief Returns the i^th row of the matrix as a vector, allows mods
   *
   *  @return		i^th row as a vector
   */
  std::vector<T>& Row(const size_t i);

  /** @brief Solves A*x=b, provided with b
   *
   *  @return		Vector solution x, optimized
   */
  std::vector<T> Solve(const std::vector<T>& b) const;

  /** @brief Returns the underlying std vector
   *
   *  @return		underlying std vector
   */
  const std::vector<std::vector<T>>& std() const;

  /** @brief Returns the matrix transpose
   *
   *  @return		The transpose of the matrix
   */
  Matrix<T> Transpose() const;

 private:

  std::vector<std::vector<T>>	_mat;	///< Underlying std::vector<std::vector>>
};

#include "Matrix-inl.hh"

#endif
