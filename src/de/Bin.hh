#ifndef BIN_HH
#define BIN_HH

// C++ includes
#include <stdlib.h>
#include <vector>
#include <cmath>

/** @brief Defines the most general n-dimensional bin.
 *
 * 	   It includes the histogram indices of the bin, the boundaries of an n-orthotope
 * 	   and its content.
 */
class Bin {
 public:
  /** @brief Default constructor, initilizes everything to 0 */
  Bin();

  /** @brief Constant constructor, sets the indices and boundaries of the bin and its content 
   *
   *  @param	indices		Indices of the bin within the histogram it is part of
   *  @param	lower		Array of lower bounds of the bin
   *  @param	upper		Array of upper bounds of the bin
   *  @param	content		Content of the bin
   */
  Bin(const std::vector<size_t>& indices,
      const std::vector<double>& lower,
      const std::vector<double>& upper,
      const double& content);

  /** @brief Copy constructor */
  Bin(const Bin& bin);

  /** @brief Equality operator */
  Bin& operator=(const Bin& bin);

  /** @brief Destructor */
  ~Bin();

  /** @brief Returns the dimension of the space in which the bin lies */
  const size_t& GetDimension() const			{ return _dim; }

  /** @brief Returns the bin indices */
  const std::vector<size_t>& GetIndices() const		{ return _indices; }

  /** @brief Returns the bin index in dimension i */
  const size_t& GetIndex(const size_t i) const		{ return _indices[i]; }

  /** @brief Returns the array of lower bounds in each dimension */
  const std::vector<double>& GetLowerBoundArray() const	{ return _lower; }

  /** @brief Returns the lower bound in dimension i */
  const double& GetLowerBound(const size_t i) const	{ return _lower[i]; }

  /** @brief Returns the array of upper bounds in each dimension */
  const std::vector<double>& GetUpperBoundArray() const	{ return _upper; }

  /** @brief Returns the upper bound in dimension i */
  const double& GetUpperBound(const size_t i) const	{ return _upper[i]; }

  /** @brief Returns the content of the bin */
  const double& GetContent() const			{ return _content; }

  /** @brief Returns the barycentre of the bin */
  std::vector<double> GetBarycentre() const;

  /** @brief Returns the centre of the bin on the ith axis */
  double GetCentre(const size_t i) const;

  /** @brief Returns the volume of the bin */
  double GetVolume() const;

  /** @brief Returns the width of the bin in a particular dimension */
  double GetWidth(const size_t i) const;

 private:
  size_t		_dim;		///< Dimension of the space
  std::vector<size_t>	_indices;	///< Indices of the bin in the histogram
  std::vector<double>	_lower;		///< Lower boundaries of the bin
  std::vector<double>	_upper;		///< Upper boundaries of the bin
  double		_content;	///< Content of the bin
};

#endif
