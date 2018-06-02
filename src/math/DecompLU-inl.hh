template<typename T> DecompLU<T>::DecompLU()
  : _fact(1), _complete(false), _P(0), _pivots(0), _L(0), _U(0) {}

template<typename T> DecompLU<T>::DecompLU(std::vector<std::vector<T>> mat)
  : _fact(1), _complete(false), _P(0), _pivots(0), _L(0), _U(0) {

  try {
    DecomposeCrout(mat);
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "DecompLU<T>::DecompLU"));
  }
}

template<typename T> DecompLU<T>::DecompLU(const DecompLU& decLU) {
  *this = decLU;
}

template<typename T> DecompLU<T>& DecompLU<T>::operator=(const DecompLU<T>& decLU) {
  if ( this == &decLU )
      return *this;

  _fact = decLU._fact;
  _complete = decLU._complete;
  _P = decLU._P;
  _pivots = decLU._pivots;
  _L = decLU._L;
  _U = decLU._U;
  return *this;
}

template<typename T> DecompLU<T>::~DecompLU() {}

template<typename T> void DecompLU<T>::DecomposeCrout(std::vector<std::vector<T>> mat) {

  // Check the input matrix
  Assert::IsNotEmpty("DecompLU<T>::DecomposeCrout", "Matrix", mat);
  Assert::IsSquare("DecompLU<T>::DecomposeCrout", "Matrix", mat);

  // Use the Crout matrix decomposition algorithm to get the rest of the elements
  size_t pivot;		// Current pivot row
  T max(0);		// Maximum possible diagonal value

  // Fill the matrices P, L and U with default values in case of singularity
  size_t n = mat.size();
  _P = NewMatrix(n, n);
  _L = NewMatrix(n, n);
  _U = NewMatrix(n, n);
  size_t i, j, k;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      _P[i][j] = ( i == j ) ? 1 : 0;
      _L[i][j] = 0;
      _U[i][j] = ( i == j ) ? 1 : 0;
    }
  }

  _pivots.resize(n);	// Pivot rows
  for (i = 0; i < n; i++)
      _pivots[i] = i;

  // For each matrix row k = 0, ..., n-1,
  for (k = 0; k < n; k++) {
 
    // Find the pivot row, the one with the highest k^th element
    pivot = k;
    max = fabs(mat[k][k]);
    for (j = k + 1; j < n; j++) {
      if ( max < fabs(mat[j][k]) ) {
        max = fabs(mat[j][k]);
        pivot = j;
      }
    }

    // If the pivot row differs from the current row, then interchange the two rows.
    if ( pivot != k ) {
      std::swap(mat[pivot], mat[k]);
      std::swap(_P[pivot], _P[k]);
      std::swap(_pivots[pivot], _pivots[k]);
      _fact *= -1;
    }

    // If the matrix is singular, return true anyway but flag incomplete
    if ( !mat[k][k] )
	return;

    // Otherwise find the upper triangular matrix elements for row k.
    for (j = k+1; j < n; j++)
         mat[k][j] /= mat[k][k];

    // Update remaining matrix
    for (i = k+1; i < n; i++)
       for (j = k+1; j < n; j++)
            mat[i][j] -= mat[i][k]*mat[k][j];

  }

  // Fill the matrices L and U
  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      // Fill the lower triangular matrix
      _L[j][i] = mat[j][i];

      // Fill the upper triangular matrix, the diag elements are 1
      _U[i][j] = ( i == j ) ? 1 : mat[i][j];
    }
  }

  _complete = true;
}

template<typename T> T DecompLU<T>::Determinant() {
  T det(1);
  size_t n = _L.size();
  for (size_t i = 0; i < n; i++)
      det *= _L[i][i];

  return _fact*det;
}

template<typename T> std::vector<std::vector<T>> DecompLU<T>::NewMatrix(size_t m, size_t n) {

  std::vector<std::vector<T>> mat(m);
  size_t i;
  for (i = 0; i < m; i++)
      mat[i].resize(n);

  return mat;
}
