template<typename T> DecompQR<T>::DecompQR()
  : _Q(0), _R(0) {}

template<typename T> DecompQR<T>::DecompQR(const std::vector<std::vector<T>>& mat)
  : _Q(0), _R(0) {

  try {
    DecomposeGramSchmidt(mat);
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "DecompQR<T>::DecompQR"));
  }
}

template<typename T> DecompQR<T>::DecompQR(const DecompQR& decQR) {
  *this = decQR;
}

template<typename T> DecompQR<T>& DecompQR<T>::operator=(const DecompQR<T>& decQR) {
  if ( this == &decQR )
      return *this;

  _Q = decQR._Q;
  _R = decQR._R;
  return *this;
}

template<typename T> DecompQR<T>::~DecompQR() {}

template<typename T> void DecompQR<T>::DecomposeGramSchmidt(const std::vector<std::vector<T>>& mat) {

  // Check the input matrix
  Assert::IsNotEmpty("DecompQR<T>::DecomposeGramSchmidt", "Matrix", mat);
  Assert::IsSquare("DecompQR<T>::DecomposeGramSchmidt", "Matrix", mat);

  // Fill the matrices P, L and U with default values in case of singularity
  size_t n = mat.size();
  _Q = NewMatrix(n, n);
  _R = NewMatrix(n, n);
  size_t i, j, k;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      _Q[i][j] = 0;
      _R[i][j] = 0;
    }
  }

  // For each matrix column k = 0, ..., n-1, find the unscaled Q matrix
  T fact, norm;
  std::vector<T> matcol(n), Qcol(n);
  for (k = 0; k < n; k++) {

    // Extract the k^th column of the input matrix and set the k^th column of Q
    for (j = 0; j < n; j++) {
      matcol[j] = mat[j][k];
      _Q[j][k] = matcol[j];
    }

    // For each i, find the i^th column of Q and project the current column on it.
    // Increment by substracting the projection from the current column
    for (i = 0; i < k; i++) {
      for (j = 0; j < n; j++)
	  Qcol[j] = _Q[j][i];
      fact = DotProduct(matcol, Qcol)/DotProduct(Qcol, Qcol);

      for (j = 0; j < n; j++)
	  _Q[j][k] -= fact*Qcol[j];
    }
  }

  // Normalize the columns
  for (k = 0; k < n; k++) {
    for (j = 0; j < n; j++)
	Qcol[j] = _Q[j][k];

    norm = sqrt(DotProduct(Qcol, Qcol));
    for (j = 0; j < n; j++)
	_Q[j][k] /= norm;
  }

  // For each matrix row k = 0, ..., n-1, compute R = Q^T*M
  for (k = 0; k < n; k++)
    for (j = k; j < n; j++)
      for (i = 0; i < n; i++)
  	  _R[k][j] += _Q[i][k]*mat[i][j];
}

template<typename T> std::vector<std::vector<T>> DecompQR<T>::NewMatrix(size_t m, size_t n) {

  std::vector<std::vector<T>> mat(m);
  size_t i;
  for (i = 0; i < m; i++)
      mat[i].resize(n);

  return mat;
}

template<typename T> T DecompQR<T>::DotProduct(const std::vector<T>& a, const std::vector<T>& b) {

  size_t i;
  T product(0);
  for (i = 0; i < a.size(); i++)
      product += a[i]*b[i];

  return product;
}
