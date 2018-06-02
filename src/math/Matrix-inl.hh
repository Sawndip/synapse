template<typename T> Matrix<T>::Matrix() :
  _mat(0) {
}

template<typename T> Matrix<T>::Matrix(const size_t n) :
  _mat(0) {

  Resize(n, n);
}

template<typename T> Matrix<T>::Matrix(const size_t m, const size_t n) :
  _mat(0) {

  Resize(m, n);
}

template<typename T> Matrix<T>::Matrix(const size_t m, const size_t n, const T& value) :
  _mat(0) {

  Resize(m, n);
  Fill(value);
}

template<typename T> Matrix<T>::Matrix(const std::vector<T>& v) :
  _mat(0) {

  Resize(v.size(), 1);
  size_t i;
  for (i = 0; i < v.size(); i++)
      _mat[i][0] = v[i];
}

template<typename T> Matrix<T>::Matrix(const std::vector<std::vector<T>>& mat) :
  _mat(mat) {
}

template<typename T> Matrix<T>::Matrix(const Matrix& mat) {
  *this = mat;
}

template<typename T> Matrix<T>::~Matrix() {}

template<typename T> Matrix<T>& Matrix<T>::operator=(const Matrix<T>& mat) {
  if ( this == &mat )
      return *this;

  _mat = mat._mat;
  return *this;
}

template<typename T> Matrix<T>::operator bool() const {

  return _mat.size();
}

template<typename T> bool Matrix<T>::operator!() const {

  return !_mat.size();
}

template<typename T> Matrix<T> Matrix<T>::operator-() const {

  // Change the sign of every elements
  Matrix<T> opp(Nrows(), Ncols());
  size_t i, j;
  for (i = 0; i < Nrows(); i++)
    for (j = 0; j < Ncols(); j++)
	opp[i][j] = -_mat[i][j];

  return opp;
}

template<typename T> Matrix<T> Matrix<T>::operator+(const Matrix<T>& mat) const {

  // Check that the sizes of the matrices to be added are compatible
  Assert::SameSize("Matrix<T>::operator+", "Matrices", _mat, mat.std());

  // Add the elements one-to-one
  Matrix<T> sum(Nrows(), Ncols());
  size_t i, j;
  for (i = 0; i < Nrows(); i++)
    for (j = 0; j < Ncols(); j++)
	sum[i][j] = _mat[i][j]+mat[i][j];

  return sum;
}

template<typename T> void Matrix<T>::operator+=(const Matrix<T>& mat) {

  // Check that the sizes of the matrices to be added are compatible
  Assert::SameSize("Matrix<T>::operator+=", "Matrices", _mat, mat.std());

  // Add the elements one-to-one
  size_t i, j;
  for (i = 0; i < _mat.size(); i++)
    for (j = 0; j < _mat[0].size(); j++)
	_mat[i][j] += mat[i][j];
}

template<typename T> Matrix<T> Matrix<T>::operator-(const Matrix<T>& mat) const {

  // Check that the sizes of the matrices to be substracted are compatible
  Assert::SameSize("Matrix<T>::operator-", "Matrices", _mat, mat.std());

  // Substract the elements one-to-one
  Matrix<T> diff(Nrows(), Ncols());
  size_t i, j;
  for (i = 0; i < Nrows(); i++)
    for (j = 0; j < Ncols(); j++)
	diff[i][j] = _mat[i][j]-mat[i][j];

  return diff;
}

template<typename T> void Matrix<T>::operator-=(const Matrix<T>& mat) {

  // Check that the sizes of the matrices to be substracted are compatible
  Assert::SameSize("Matrix<T>::operator-=", "Matrices", _mat, mat.std());

  // Substract the elements one-to-one
  size_t i, j;
  for (i = 0; i < _mat.size(); i++)
    for (j = 0; j < _mat[0].size(); j++)
	_mat[i][j] -= mat[i][j];
}

template<typename T> Matrix<T> Matrix<T>::operator*(const Matrix<T>& mat) const {

  // Check that the sizes of the matrices to be multiplied are compatible
  Assert::Multiplicable("Matrix<T>::operator*", "Matrices", _mat, mat.std());

  // Multiply the two matrices
  Matrix<T> prod(Nrows(), mat.Ncols());
  size_t i, j, k;
  for (i = 0; i < Nrows(); i++)
    for (j = 0; j < mat.Ncols(); j++)
      for (k = 0; k < Ncols(); k++)
	  prod[i][j] += _mat[i][k]*mat[k][j];

  return prod;
}

template<typename T> void Matrix<T>::operator*=(const Matrix<T>& mat) {

  // Check that the sizes of the matrices to be multiplied are compatible
  Assert::Multiplicable("Matrix<T>::operator*=", "Matrices", _mat, mat.std());

  // Update the parent matrix
  std::vector<std::vector<T>> temp = _mat;
  Resize(Nrows(), mat.Ncols());
  size_t i, j, k;
  for (i = 0; i < temp.size(); i++)
    for (j = 0; j < mat.Ncols(); j++) {
      _mat[i][j] = 0; 
      for (k = 0; k < temp[i].size(); k++)
	  _mat[i][j] += temp[i][k]*mat[k][j];
    }
}

template<typename T> Matrix<T> Matrix<T>::operator*(const T& factor) const {

  // Multiply the matrix with a scalar to its right
  Matrix<T> prod(Nrows(), Ncols());
  size_t i, j;
  for (i = 0; i < Nrows(); i++)
    for (j = 0; j < Ncols(); j++)
	prod[i][j] = factor*_mat[i][j];

  return prod;
}

template<typename T> void Matrix<T>::operator*=(const T& factor) {

  // Multiply in place the matrix with a scalar to its right
  size_t i, j;
  for (i = 0; i < Nrows(); i++)
    for (j = 0; j < Ncols(); j++)
	_mat[i][j] *= factor;
}

template<typename T> Matrix<T> operator*(const T& factor, const Matrix<T> mat) {

  // Multiply the matrix with a scalar to its left
  Matrix<T> prod(mat.Nrows(), mat.Ncols());
  size_t i, j;
  for (i = 0; i < mat.Nrows(); i++)
    for (j = 0; j < mat.Ncols(); j++)
	prod[i][j] = factor*mat[i][j];

  return prod;
}

template<typename T> Vector<T> Matrix<T>::operator*(const Vector<T>& vec) const {

  // Check that the sizes of the matrix and vector to be multiplied are compatible
  Assert::SameSize("Matrix<T>::operator*", "Matrix and vector", _mat[0], vec.std());

  // Multiply the matrix with a vector to its right
  Vector<T> prod(Nrows());
  size_t i, j;
  for (i = 0; i < Nrows(); i++)
    for (j = 0; j < Ncols(); j++)
	prod[i] += _mat[i][j]*vec[j];

  return prod;
}

template<typename T> Matrix<T> Matrix<T>::operator/(const T& factor) const {

  // Check that the divisor is non-zero
  Assert::IsNonZero("Matrix<T>::operator/", "Divisor", factor);

  // Multiply the matrix with a scalar to its right
  Matrix<T> ratio(Nrows(), Ncols());
  size_t i, j;
  for (i = 0; i < Nrows(); i++)
    for (j = 0; j < Ncols(); j++)
	ratio[i][j] = _mat[i][j]/factor;

  return ratio;
}

template<typename T> void Matrix<T>::operator/=(const T& factor) {

  // Check that the divisor is non-zero
  Assert::IsNonZero("Matrix<T>::operator/=", "Divisor", factor);

  // Multiply the matrix with a scalar to its right
  size_t i, j;
  for (i = 0; i < Nrows(); i++)
    for (j = 0; j < Ncols(); j++)
	_mat[i][j] /= factor;
}

template<typename T> const std::vector<T>& Matrix<T>::operator[](const size_t i) const {

  Assert::IsWithinBounds("Matrix<T>::operator[]", Nrows(), i);
  return _mat[i];
}

template<typename T> std::vector<T>& Matrix<T>::operator[](const size_t i) {

  Assert::IsWithinBounds("Matrix<T>::operator[]", Nrows(), i);
  return _mat[i];
}

template<typename T> std::ostream& operator<<(std::ostream& os, const Matrix<T>& mat) {

  size_t i, j;
  os << "[[";
  for (i = 0; i < mat.Nrows(); i++) {
    if ( i )
        os << "[";
    for (j = 0; j < mat.Ncols(); j++) {
      os << mat[i][j];
      if ( j < mat.Ncols()-1 )
	os << ", ";
    }
    if ( i < mat.Nrows()-1 )
	os << "], ";
  }
  os << "]]";
  return os;
}

template<typename T> T Matrix<T>::Cofactor(const size_t i, const size_t j) const {

  // Check that the matrix is square and of at least dimension 1
  Assert::IsSquare("Matrix<T>::Cofactor", "Matrix", _mat);
  Assert::IsGreater("Matrix<T>::Cofactor", "Dimension", Nrows(), size_t(2));

  // Copy the matrix into a new vector of vector and simply remove the undesired elements
  std::vector<std::vector<T>> sub = _mat;
  sub.erase(sub.begin()+i);
  size_t k;
  for (k = 0; k < sub.size(); k++)
      sub[k].erase(sub[k].begin()+j);

  // The cofactor is the determinant of sub with a prefactor
  return pow(-1, i+j)*Matrix<T>(sub).Determinant();
}

template<typename T> Matrix<T> Matrix<T>::CofactorMatrix() const {

  // Define a matrix of dimension n*n and fill it with cofactors
  Matrix<T> cof(Nrows(), Ncols());
  try {
    size_t i, j;
    for (i = 0; i < Nrows(); i++)
      for (j = 0; j < Ncols(); j++)
	  cof[i][j] = Cofactor(i, j);
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "Matrix<T>::CofactorMatrix"));
  }

  return cof;
}

template<typename T> std::vector<T> Matrix<T>::Column(const size_t i) const {

  Assert::IsWithinBounds("Matrix<T>::Column", Ncols(), i);
  std::vector<T> column(Nrows());
  size_t j;
  for (j = 0; j < Nrows(); j++)
      column[j] = _mat[j][i];

  return column;
}

template<typename T> T Matrix<T>::Determinant() const {

  // Cannot compute the determinant of a non square matrix
  Assert::IsNotEmpty("Matrix<T>::Determinant", "Matrix", _mat);
  Assert::IsSquare("Matrix<T>::Determinant", "Matrix", _mat);

  // Get rid of the simple cases without using the heavy artillery
  if ( Nrows() == 1 ) {
    return _mat[0][0];
  } else if ( Nrows() == 2 ) {
    return _mat[0][0]*_mat[1][1] - _mat[0][1]*_mat[1][0];
  } else if ( Nrows() == 3 ) {
    return   _mat[0][0]*(_mat[1][1]*_mat[2][2]-_mat[1][2]*_mat[2][1])
	   - _mat[0][1]*(_mat[1][0]*_mat[2][2]-_mat[1][2]*_mat[2][0])
	   + _mat[0][2]*(_mat[1][0]*_mat[2][1]-_mat[1][1]*_mat[2][0]);
  }

  // For larger matrices, use LU decomposition, fastest
  DecompLU<T> decLU(_mat);
  return decLU.Determinant();
}

template<typename T> Matrix<T> Matrix<T>::EigenVectors(std::vector<T>& lambdas,
						       const double eps) const {

  // Check that the matrix is square
  Assert::IsNotEmpty("Matrix<T>::EigenVectors", "Matrix", _mat);
  Assert::IsSquare("Matrix<T>::EigenVectors", "Matrix", _mat);

  // Initialize the containers for the decomposition
  size_t dim = _mat.size();			// Dimension of the square matrix
  size_t i;

  std::vector<T> tempvals(dim);			// Array to store the temporary eigenvalues
  for (i = 0; i < dim; i++)
      tempvals[i] = 1.;
  std::vector<T> eigenvals(dim);

  std::vector<std::vector<T>> Mk = _mat; 	// k^th iteration of the similar matrix
  Matrix<T> Qk, Q(dim, dim, 0);			// k^th iteration of the Q matrix and its composite
  for (i = 0; i < dim; i++)
      Q[i][i] = 1;

  // QR decompose the Matrix M, compute Qk^T*M*Qk, repeat until convergence
  try {
    bool done = false;
    while ( !done ) {
      DecompQR<T> decQR(Mk);
      Qk = Matrix<T>(decQR.Q());
      Q *= Qk;
      Mk = (Qk.Transpose()*Matrix<T>(Mk)*Qk).std();

      // Update and check the evolution of the eigenvalues
      done = true;
      for (i = 0; i < Mk.size(); i++) {
        eigenvals[i] = Mk[i][i];
        if ( fabs((eigenvals[i]-tempvals[i])/tempvals[i]) > eps )
	    done = false;
      }
      tempvals = eigenvals;
    }
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not eigendecompose"+std::string(e.what()),
	  "Matrix<T>::EigenVectors"));
  }

  // Return the eigen values and vectors
  lambdas = eigenvals;
  return Q;
}

template<typename T> void Matrix<T>::Fill(const T& value) {

  size_t i, j;
  for (i = 0; i < Nrows(); i++)
    for (j = 0; j < Ncols(); j++)
	_mat[i][j] = value;
}

template<typename T> void Matrix<T>::Identity() {

  Assert::IsSquare("Matrix<T>::Identity", "Matrix", _mat);
  size_t i, j;
  for (i = 0; i < Nrows(); i++)
    for (j = 0; j < Ncols(); j++)
	_mat[i][j] = (i == j) ? 1 : 0;
}
 
template<typename T> Matrix<T> Matrix<T>::Inverse() const {

  // Initialize a matrix from the base matrix, invert it in place and return
  Matrix<T> inv(_mat);
  inv.Invert();
  return inv; 
}

template<typename T> void Matrix<T>::Invert() {

  // Check that the matrix to invert is square and non-zero dimensional
  Assert::IsNotEmpty("Matrix<T>::Invert", "Matrix", _mat);
  Assert::IsSquare("Matrix<T>::Invert", "Matrix", _mat);

  // Initialize the augmented matrix by adding an identity to it
  size_t i, j, k;
  size_t n = Nrows();
  for (i = 0; i < n; i++) {
    _mat[i].resize(2*n);
    _mat[i][n+i] = 1;
  }

  // Perform the Gauss-Jordan elimination
  size_t pivot;
  T max, temp;
  for (k = 0; k < n; k++) {
 
    // Find the pivot column, the one with the highest k^th element on the
    pivot = k;
    max = fabs(_mat[k][k]);
    for (j = k+1; j < n; j++) {
      if ( max < fabs(_mat[j][k]) ) {
        max = fabs(_mat[j][k]);
        pivot = j;
      }
    }

    // If the the pivot element is 0, return that the matrix is singular, increment otherwise
    if ( !max ) {
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
			"The matrix is singular", "Matrix<T>::Invert"));
      return;
    }

    // Swap the k^th and pivot^th lines if not identical
    if ( pivot != k )
    	std::swap(_mat[pivot], _mat[k]);

    // Divide the pivot line by its element on the diagonal
    temp = _mat[k][k];
    for (i = 0; i < 2*n; i++)
        _mat[k][i] /= temp;

    // Simplify all the lines below the current pivot row
    for (i = 0; i < n; i++) {
      if ( i != k ) {
	temp = _mat[i][k];
        for (j = k+1; j < 2*n; j++)
	    _mat[i][j] -= temp*_mat[k][j];
      }
      _mat[i][k] = 0;
    }
  }

  // Remove the first n columns of the augmented matrix
  for (i = 0; i < n; i++)
      _mat[i].erase(_mat[i].begin(), _mat[i].begin()+n);
}

template<typename T> bool Matrix<T>::IsSymmetric() const {

  Assert::IsSquare("Matrix<T>::IsSymmetric", "Matrix", _mat);
  size_t i, j;
  for (i = 0; i < Nrows(); i++)
    for (j = i+1; j < Ncols(); j++)
      if ( _mat[i][j] != _mat[j][i] )
	  return false;

  return true;
}

template<typename T> const size_t Matrix<T>::Nrows() const {
  return _mat.size();
}

template<typename T> const size_t Matrix<T>::Ncols() const {
  if ( !_mat.size() ) {
    return _mat.size();
  } else {
    return _mat[0].size();
  }
}

template<typename T> void Matrix<T>::Print() const {
  std::cerr << "Matrix (" << Nrows() << ", " << Ncols()  << ") is as follows" << std::endl;
  size_t i, j;
  for (i = 0; i < Nrows(); i++) {
    for (j = 0; j < Ncols(); j++)
	std::cerr << _mat[i][j] << "\t";
    std::cerr << std::endl;
  }
}

template<typename T> void Matrix<T>::RemoveColumn(const size_t i) {

  Assert::IsWithinBounds("Matrix<T>::RemoveColumn", Ncols(), i);
  size_t j;
  for (j = 0; j < Nrows(); j++)
      _mat[j].erase(_mat[j].begin()+i);
}

template<typename T> void Matrix<T>::RemoveRow(const size_t i) {

  Assert::IsWithinBounds("Matrix<T>::RemoveRow", Nrows(), i);
  _mat.erase(_mat.begin()+i);
}

template<typename T> void Matrix<T>::Resize(const size_t m, const size_t n) {

  _mat.resize(m);
  size_t i;
  for (i = 0; i < m; i++)
      _mat[i].resize(n);
}

template<typename T> const std::vector<T>& Matrix<T>::Row(const size_t i) const {

  Assert::IsWithinBounds("Matrix<T>::Row", Nrows(), i);
  return _mat[i];
}

template<typename T> std::vector<T>& Matrix<T>::Row(const size_t i) {

  Assert::IsWithinBounds("Matrix<T>::Row", Nrows(), i);
  return _mat[i];
}

template<typename T> std::vector<T> Matrix<T>::Solve(const std::vector<T>& b) const {

  // Check that the matrix to invert is square and non-zero dimensional
  Assert::IsNotEmpty("Matrix<T>::Solve", "Matrix", _mat);
  Assert::IsSquare("Matrix<T>::Solve", "Matrix", _mat);

  // Find the LU decomposition of the matrix. If the matrix is singular, do not bother
  size_t n = Nrows();
  DecompLU<T> decLU(_mat);
  if ( !decLU.Complete() ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
		"The matrix is singular", "Matrix<T>::Solve"));
    return std::vector<T>(n);
  }
  std::vector<std::vector<T>> L(decLU.L()), U(decLU.U());
  std::vector<T> pivots(decLU.Pivots());

  // Take care of the permutation, define a new b'
  std::vector<T> bprime(n);
  size_t i;
  for (i = 0; i < n; i++)
      bprime[i] = b[pivots[i]];

  // Solve Ld = b'
  std::vector<T> x(n);
  T s;
  int j, k;
  for (k = 0; k < (int)n; k++) {
    s = 0;
    for (j = 0; j < k; j++)
	s += L[k][j]*x[j];
    x[k] = (bprime[k] - s)/L[k][k];
  }

  // Back substitution to solve Ux = d
  for (k = n-1; k > -1; k--) {
    s = 0;
    for (j = k+1; j < (int)n; j++)
	s += U[k][j]*x[j];
    x[k] = x[k] - s;
  }

  return x;
}

template<typename T> const std::vector<std::vector<T>>& Matrix<T>::std() const {

  return _mat;
}

template<typename T> Matrix<T> Matrix<T>::Transpose() const {

  Matrix<T> trans(Ncols(), Nrows());
  size_t i, j;
  for (i = 0; i < Ncols(); i++)
    for (j = 0; j < Nrows(); j++)
	trans[i][j] = _mat[j][i];

  return trans;
}
