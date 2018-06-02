template<typename T> Vector<T>::Vector() :
  _vec(0) {
}

template<typename T> Vector<T>::Vector(const size_t& n) :
  _vec(std::vector<T>(n)) {
}

template<typename T> Vector<T>::Vector(const std::vector<T>& vec) :
  _vec(vec) {
}

template<typename T> Vector<T>::Vector(const T* vec, const size_t& n) :
  _vec() {

  size_t i;
  for (i = 0; i < n; i++)
      _vec.push_back(vec[i]); 
}

template<typename T> Vector<T>::Vector(const size_t& n, const T& value) :
  _vec(std::vector<T>(n)) {

  std::fill(_vec.begin(), _vec.end(), value);
}

template<typename T> Vector<T>::Vector(const Vector& vec) {
  *this = vec;
}

template<typename T> Vector<T>::~Vector() {}

template<typename T> Vector<T>& Vector<T>::operator=(const Vector<T>& vec) {
  if ( this == &vec )
      return *this;

  _vec = vec._vec;
  return *this;
}

template<typename T> Vector<T> Vector<T>::operator-() const {
  // Change the sign of every element
  Vector<T> opp(size());
  size_t i;
  for (i = 0; i < size(); i++)
      opp[i] = -_vec[i];

  return opp;
}

template<typename T> Vector<T> Vector<T>::operator+(const Vector<T>& vec) const {

  // Check that the size of the vector to be added is compatible
  Assert::SameSize("Vector<T>::operator+", "The two vectors", _vec, vec.std());

  // Add the elements one to one
  Vector<T> sum(_vec.size());
  size_t i;
  for (i = 0; i < _vec.size(); i++)
	sum[i] = _vec[i]+vec[i];

  return sum;
}

template<typename T> void Vector<T>::operator+=(const Vector<T>& vec) {

  // Check that the size of the vector to be added is compatible
  Assert::SameSize("Vector<T>::operator+=", "The two vectors", _vec, vec.std());

  // Add the elements one to one
  size_t i;
  for (i = 0; i < _vec.size(); i++)
      _vec[i] += vec[i];
}

template<typename T> Vector<T> Vector<T>::operator-(const Vector<T>& vec) const {

  // Check that the size of the vector to be substracted is compatible
  Assert::SameSize("Vector<T>::operator-", "The two vectors", _vec, vec.std());

  // Add the elements one to one
  Vector<T> diff(_vec.size());
  size_t i;
  for (i = 0; i < _vec.size(); i++)
	diff[i] = _vec[i]-vec[i];

  return diff;
}

template<typename T> void Vector<T>::operator-=(const Vector<T>& vec) {

  // Check that the size of the vector to be substracted is compatible
  Assert::SameSize("Vector<T>::operator-=", "The two vectors", _vec, vec.std());

  // Add the elements one to one
  size_t i;
  for (i = 0; i < _vec.size(); i++)
      _vec[i] -= vec[i];
}

template<typename T> T Vector<T>::operator*(const Vector<T>& vec) const {

  // Check that the size of the vector to be multipled is compatible
  Assert::SameSize("Vector<T>::operator*", "The two vectors", _vec, vec.std());

  // Multiply the two vectors
  T prod(0);
  size_t i;
  for (i = 0; i < _vec.size(); i++)
      prod += _vec[i]*vec[i];

  return prod;
}

template<typename T> Vector<T> Vector<T>::operator*(const T& factor) const {

  // Multiply each of the vector element by the factor
  Vector<T> prod(_vec.size());
  size_t i;
  for (i = 0; i < _vec.size(); i++)
      prod[i] = factor*_vec[i];

  return prod;
}

template<typename T> void Vector<T>::operator*=(const T& factor) {

  // Multiply each of the vector element by the factor
  size_t i;
  for (i = 0; i < _vec.size(); i++)
      _vec[i] *= factor;
}

template<typename T> Vector<T> Vector<T>::operator/(const T& factor) const {

  // Check that the divisor is non-zero
  Assert::IsNonZero("Vector<T>::operator/", "Divisor", factor);

  // Divide each of the vector element by the factor
  Vector<T> div(_vec.size());
  size_t i;
  for (i = 0; i < _vec.size(); i++)
      div[i] = _vec[i]/factor;

  return div;
}

template<typename T> void Vector<T>::operator/=(const T& factor) {

  // Check that the divisor is non-zero
  Assert::IsNonZero("Vector<T>::operator/=", "Divisor", factor);

  // Multiply each of the vector element by the factor
  size_t i;
  for (i = 0; i < _vec.size(); i++)
      _vec[i] /= factor;
}

template<typename T> Vector<T> operator*(const T& factor, const Vector<T> vec) {

  // Multiply the two vecrices
  Vector<T> prod(vec.size());
  size_t i;
  for (i = 0; i < vec.size(); i++)
      prod[i] = factor*vec[i];

  return prod;
}

template<typename T> const T& Vector<T>::operator[](const size_t i) const {

  Assert::IsWithinBounds("Vector<T>::operator[]", size(), i);
  return _vec[i];
}

template<typename T> T& Vector<T>::operator[](const size_t i) {

  Assert::IsWithinBounds("Vector<T>::operator[]", size(), i);
  return _vec[i];
}

template<typename T> std::ostream& operator<<(std::ostream& os, const Vector<T>& vec) {

  size_t i;
  os << "[";
  for (i = 0; i < vec.size()-1; i++)
      os << vec[i] << ", ";
  os << vec.back() << "]";
  return os;
}

template<typename T> T Vector<T>::mag() const {

  return sqrt(mag2());
}

template<typename T> T Vector<T>::mag2() const {

  T mag2(0);
  size_t i;
  for (i = 0; i < _vec.size(); i++)
      mag2 += _vec[i]*_vec[i];

  return mag2;
}

template<typename T> T Vector<T>::dot(const Vector& vec) const {

  // Check that the size of the vector to be multipled is compatible
  Assert::SameSize("Vector<T>::dot", "The two vectors", _vec, vec.std());

  // Multiply the two vectors
  T prod(0);
  size_t i;
  for (i = 0; i < _vec.size(); i++)
      prod += _vec[i]*vec[i];

  return prod;
}

template<typename T> Vector<T> Vector<T>::cross(const Vector& vec) const {

  // Check that the size of the vector to be multipled is compatible
  if ( vec.size() != 3 || vec.size() != _vec.size() )
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	               		  "The vectors are not of dimension 3",
				  "Vector<T>::cross"));

  // Multiply the two vectors
  Vector<T> cross(3);
  cross[0] = _vec[1]*vec[2]-_vec[2]*vec[1];
  cross[1] = _vec[2]*vec[0]-_vec[0]*vec[2];
  cross[2] = _vec[0]*vec[1]-_vec[1]*vec[0];

  return cross;
}

template<typename T> const T& Vector<T>::back() const {
  return _vec.back();
}

template<typename T> const T& Vector<T>::front() const {
  return _vec.front();
}

template<typename T> T Vector<T>::norm() const {
  return this->mag();
}

template<typename T> T Vector<T>::norm2() const {
  return this->mag2();
}

template<typename T> void Vector<T>::resize(const size_t n) {
  _vec.resize(n);
}

template<typename T> void Vector<T>::push_back(const T& element) {
  _vec.push_back(element);
}

template<typename T> void Vector<T>::push_front(const T& element) {
  _vec.push_front(element);
}

template<typename T> const size_t Vector<T>::size() const {
  return _vec.size();
}

template<typename T> const std::vector<T>& Vector<T>::std() const {
  return _vec;
}
