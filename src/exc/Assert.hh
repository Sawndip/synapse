#ifndef ASSERT_HH
#define ASSERT_HH

// C++ includes
#include <iostream>
#include <vector>
#include <string>
#include <map>

// Other includes
#include "Exception.hh"

/** @brief List of common assertions that throw an Exceptions::Exception if not satistified.
 *
 * 	   Simplifies the code down the line. All the assertions take two common inputs:
 *	    - the function that called the assert (location);
 *	    - the name of the object to check (what).
 */
namespace Assert {

/** @brief Asserts that a value is non zero
 *
 *  @param	location	Name of the function that called the assert
 *  @param	what		Name of the object to check
 *  @param	value		Value to check
 */
template<typename T> void IsNonZero(const std::string& location,
	       			    const std::string& what,
 	       			    const T& value);

/** @brief Asserts that a value is greater than another
 *
 *  @param	location	Name of the function that called the assert
 *  @param	what		Name of the object to check
 *  @param	value		Value to check
 *  @param	limit		Minimum value
 *  @param	strictly	If true, must be strictly greater
 */
template<typename T> void IsGreater(const std::string& location,
	       			    const std::string& what,
 	       			    const T& value,
 	       			    const T& limit,
				    const bool strictly=false);

/** @brief Asserts that a value is lower than another
 *
 *  @param	location	Name of the function that called the assert
 *  @param	what		Name of the object to check
 *  @param	value		Value to check
 *  @param	limit		Maximum value
 *  @param	strictly	If true, must be strictly lower
 */
template<typename T> void IsLower(const std::string& location,
	       			  const std::string& what,
 	       			  const T& value,
 	       			  const T& limit,
				  const bool strictly=false);

/** @brief Asserts that the values are equal
 *
 *  @param	location	Name of the function that called the assert
 *  @param	what		Name of the object to check
 *  @param	a		First value
 *  @param	b		Second value
 */
template<typename T> void IsEqual(const std::string& location,
	       			  const std::string& what,
 	       			  const T& a,
 	       			  const T& b);

/** @brief Asserts that a vector is not empty
 *
 *  @param	location	Name of the function that called the assert
 *  @param	what		Name of the object to check
 *  @param	vector		Vector to check
 */
template<typename T> void IsNotEmpty(const std::string& location,
				     const std::string& what,
				     const std::vector<T>& vector);

/** @brief Asserts that the index that is accessed is within bounds
 *
 *  @param	location	Name of the function that called the assert
 *  @param	n		Size of the vector
 *  @param	i		Index
 */
void IsWithinBounds(const std::string& location,
		    const size_t n,
		    const size_t i);

/** @brief Asserts that a vector has at least n elements
 *
 *  @param	location	Name of the function that called the assert
 *  @param	what		Name of the object to check
 *  @param	vector		Vector to check
 *  @param	n		Minimum number of elements
 */
template<typename T> void HasAtLeast(const std::string& location,
				     const std::string& what,
				     const std::vector<T>& vector,
				     const size_t n);

/** @brief Asserts that a vector has at most n elements
 *
 *  @param	location	Name of the function that called the assert
 *  @param	what		Name of the object to check
 *  @param	vector		Vector to check
 *  @param	n		Maximum number of elements
 */
template<typename T> void HasAtMost(const std::string& location,
				    const std::string& what,
				    const std::vector<T>& vector,
				    const size_t n);
/** @brief Asserts that two vectors have the same size
 *
 *  @param	location	Name of the function that called the assert
 *  @param	what		Name of the object to check
 *  @param	va		First vector
 *  @param	vb		Second vector
 */
template<typename T> void SameSize(const std::string& location,
				   const std::string& what,
				   const std::vector<T>& va,
				   const std::vector<T>& vb);

/** @brief Asserts that two matrices have the same size
 *
 *  @param	location	Name of the function that called the assert
 *  @param	what		Name of the object to check
 *  @param	ma		First matrix
 *  @param	mb		Second matrix
 */
template<typename T> void SameSize(const std::string& location,
				   const std::string& what,
				   const std::vector<std::vector<T>>& ma,
				   const std::vector<std::vector<T>>& mb);

/** @brief Asserts that two matrices are multiplicable
 *
 *  @param	location	Name of the function that called the assert
 *  @param	what		Name of the object to check
 *  @param	ma		First matrix
 *  @param	mb		Second matrix
 */
template<typename T> void Multiplicable(const std::string& location,
				   	const std::string& what,
				   	const std::vector<std::vector<T>>& ma,
				   	const std::vector<std::vector<T>>& mb);

/** @brief Asserts that a matrix is square
 *
 *  @param	location	Name of the function that called the assert
 *  @param	what		Name of the object to check
 *  @param	matrix		Matrix to check
 */
template<typename T> void IsSquare(const std::string& location,
				   const std::string& what,
				   const std::vector<std::vector<T>>& matrix);

/** @brief Asserts that a matrix is symmetric
 *
 *  @param	location	Name of the function that called the assert
 *  @param	what		Name of the object to check
 *  @param	matrix		Matrix to check
 */
template<typename T> void IsSymmetric(const std::string& location,
				      const std::string& what,
				      const std::vector<std::vector<T>>& matrix);

/** @brief Asserts that a map contains a specific entry
 *
 *  @param	location	Name of the function that called the assert
 *  @param	what		Name of the object to check
 *  @param	map		Map to check
 *  @param	key		Key to look for
 */
template<typename T1, typename T2> void Contains(const std::string& location,
				      		 const std::string& what,
				      		 const std::map<T1, T2>& map,
				    		 const T1& key);

/** @brief Asserts that the value represents a probability (0 < p < 1)
 *
 *  @param	location	Name of the function that called the assert
 *  @param	what		Name of the object to check
 *  @param	p		Value to check
 *  @param	nonzero		Check that it is non zero
 *  @param	noncert		Check that it is non certain
 */
void IsProbability(const std::string& location,
		   const std::string& what,
		   const double& p,
		   const bool nonzero=true,
		   const bool noncert=true);

#include "Assert-inl.hh"

} // namespace Assert

#endif
