#ifndef SIMPLEXINTERPOLATOR_HH
#define SIMPLEXINTERPOLATOR_HH

// C++ includes
#include <iostream>
#include <ctime>
#include <cmath>
#include <vector>

// Additional includes
#include "Matrix.hh"
#include "Vertex.hh"

/** @brief Performs the simplical interplation between n+1 discrete samples in n dimensions.
 */
class SimplexInterpolator {
 public:
  /** @brief Default constructor, initilizes everything to 0 */
  SimplexInterpolator();

  /** @brief Normal constructor, sets the vector of vertices and their mapping
   *
   *  @param	vertices	n+1 simplex n-vertices
   */
  SimplexInterpolator(const std::vector<Vertex>& vertices);

  /** @brief Copy constructor */
  SimplexInterpolator(const SimplexInterpolator& sint);

  /** @brief Equality operator */
  SimplexInterpolator& operator=(const SimplexInterpolator& sint);

  /** @brief Function operator, executes evaluate */
  double operator() (const std::vector<double>& v) const	{ return Evaluate(v); }

  /** @brief Destructor */
  ~SimplexInterpolator();

  /** @brief Fit the n+1 vertices with a hyperplane, store the parameters
   *
   *  @param	vertices	n+1 simplex n-points
   *
   *  @return			True if successful
   */
  void Initialize(const std::vector<Vertex>& vertices);

  /** @brief Performs the simplex interpolation at an n-point v
   *
   *  @param	v	n-vector of coordinates of the n-point
   *
   *  @return		Interpolated value of the function in v
   */
  double Evaluate(const std::vector<double>& v) const;

 private:

  std::vector<double>	_params;	///< (n+1)-vector of fitting parameters
};

#endif
