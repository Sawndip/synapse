#ifndef LINEARINTERPOLATOR_HH
#define LINEARINTERPOLATOR_HH

// C++ includes
#include <iostream>
#include <ctime>
#include <cmath>
#include <vector>

// Additional includes
#include "Grid.hh"

/** @brief Performs the linear interpolation between discrete samples
 *	   of a scalar field organised on an n-orthotope.
 */
class LinearInterpolator {
 public:
  /** @brief Default constructor, initilizes everything to 0 */
  LinearInterpolator();

  /** @brief Normal constructor, sets matrix of vertices and their mapping
   *  @param	grid	2^n vertices of an n-orthotope
   */
  LinearInterpolator(const Grid& grid);

  /** @brief Copy constructor */
  LinearInterpolator(const LinearInterpolator& lint);

  /** @brief Equality operator */
  LinearInterpolator& operator=(const LinearInterpolator& lint);

  /** @brief Function operator, executes evaluate */
  double operator() (const std::vector<double>& v) const	{ return Evaluate(v); }

  /** @brief Destructor */
  ~LinearInterpolator();

  /** @brief Fit the 2^n vertices with an (n-1)-spline, store the parameters
   *
   *  @return		True if successful
   */
  void Initialize(const Grid& grid);

  /** @brief Performs the linear interpolation at an n-point P
   *
   *  @param	v	Vector of coordinates of the n-point
   *
   *  @return		Interpolated value of the function in P
   */
  double Evaluate(const std::vector<double>& v) const;

 private:

  size_t		_dim;		///< Dimension of the space
  std::vector<double>	_params;	///< 2^n-matrix vector
};

#endif
