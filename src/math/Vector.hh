#ifndef VECTOR_HH
#define VECTOR_HH

// C++ includes
#include <iostream>
#include <cmath>
#include <vector>

// Other includes
#include "Assert.hh"

/** @brief Object that can host any type of data in a vector arrangement. 
 *
 * 	   The class definition also contains a multitude of linear algebra methods
 * 	   that are useful to work with vectors.
 */
template<typename T> class Vector {
 public:
  /** @brief Default constructor, initilizes an empty vector of dimension 0 */
  Vector();

  /** @brief Size constructor, initializes an empty vector of dimension n */
  Vector(const size_t& n);

  /** @brief Std constructor, sets the underlying std::vector */
  Vector(const std::vector<T>& vec);

  /** @brief Array constructor, sets the underlying std::vector */
  Vector(const T* vec, const size_t& n);

  /** @brief Constant constructor, initializes a vector of dimension n filled with value */
  Vector(const size_t& n, const T& value);

  /** @brief Copy constructor */
  Vector(const Vector& vec);

  /** @brief Destructor */
  ~Vector();

  /** @brief Equality operator */
  Vector& operator=(const Vector& vec);

  /** @brief Overload the unary minus for vector algebra */
  Vector operator-() const;

  /** @brief Overload the vector sum */
  Vector operator+(const Vector& vec) const;

  /** @brief Overload the vector sum */
  void operator+=(const Vector& vec);

  /** @brief Overload the v substraction */
  Vector operator-(const Vector& vec) const;

  /** @brief Overload the v substraction */
  void operator-=(const Vector& vec);

  /** @brief Overload the vector product (dot product) */
  T operator*(const Vector& vec) const;

  /** @brief Overload the vector product with a scalar */
  Vector operator*(const T& factor) const;

  /** @brief Overload the vector product with a scalar */
  void operator*=(const T& factor);

  /** @brief Overload the vector ratio with a scalar */
  Vector operator/(const T& factor) const;

  /** @brief Overload the vector ratio with a scalar */
  void operator/=(const T& factor);

  /** @brief Overload the subscript operator, return a value */
  const T& operator[](const size_t i) const;

  /** @brief Overload the subscript operator, return a value, allows mods */
  T& operator[](const size_t i);

  /** @brief Returns the last element of the vector */
  const T& back() const;

  /** @brief Returns the first element of the vector */
  const T& front() const;

  /** @brief Returns the magnitude of the vector */
  T mag() const;

  /** @brief Returns the square of the magnitude of the vector */
  T mag2() const;

  /** @brief Returns the norm of the vector */
  T norm() const;

  /** @brief Returns the square of the norm of the vector */
  T norm2() const;

  /** @brief Returns the dot product with another vecror */
  T dot(const Vector& vec) const;

  /** @brief Returns the cross product with another vecror */
  Vector cross(const Vector& vec) const;

  /** @brief Adds element at the back */
  void push_back(const T& element);

  /** @brief Adds element at the back */
  void push_front(const T& element);

  /** @brief Resizes the vector to chosen dimension */
  void resize(const size_t n);

  /** @brief Returns dimension of the vector */
  const size_t size() const;

  /** @brief Returns the vector as an std::vector */
  const std::vector<T>& std() const;

 private:

  std::vector<T>	_vec;	///< Underlying std::vector
};

#include "Vector-inl.hh"

#endif
