#ifndef INTERPOLATOR_HH
#define INTERPOLATOR_HH

// C++ includes
#include <iostream>
#include <ctime>
#include <cmath>
#include <vector>

// Additional includes
#include "Matrix.hh"
#include "LinearInterpolator.hh"
#include "SimplexInterpolator.hh"
#include "Delaunay.hh"

/** @brief Interpolates between discrete samples of a scalar field.
 *
 * 	   Interpolates in x in an n-dimensional space provided with either of the following:
 * 	    - n-rectangular Grid of regularly space-points and their mapping (LinearInterpolator);
 *	    - set of randomly placed points and their mapping (SimplexInterpolator).
 */
class Interpolator {
 public:
  /** @brief Default constructor, initilizes everything to 0 */
  Interpolator();

  /** @brief Grid constructor, sets matrix of vertices and their mapping on a grid 
   *
   *  @param	grid		Grid of N_1*...*N_n vertices and their mapping, n = dimension
   *  @param 	ex		Extrapolate the outer cells to infinity if requested
   *  @param 	algo		Interpolation algorithm
   */
  Interpolator(const Grid& grid,
	       const bool ex=true,
	       const std::string algo="linear");

  /** @brief Cloud constructor, sets vector of vertices and their mapping for a point cloud
   *
   *  @param	vertices	N_1*...*N_n vector of vertices, n = dimension
   *  @param 	ex		Extrapolate outside the outmost points if requested
   *  @param 	algo		Interpolation algorithm
   */
  Interpolator(const std::vector<Vertex>& vertices,
	       const bool ex=true,
	       const std::string algo="simplex");

  /** @brief Copy constructor */
  Interpolator(const Interpolator& inter);

  /** @brief Equality operator */
  Interpolator& operator=(const Interpolator& inter);

  /** @brief Overloaded function operator, calls Evaluate
   *
   *  @param	v	Vector of coordinates
   */
  double operator()(const std::vector<double>& v) const	{ return Evaluate(v); }

  /** @brief Overloaded function operator, calls Evaluate
   *
   *  @param	v	One dimensional value
   */
  double operator()(const double& v) const		{ return Evaluate(v); }

  /** @brief Overloaded function operator, calls Evaluate
   *
   *  @param	v	Pointer to an array of coordinates
   */
  double operator()(const double* v) const		{ return Evaluate(v); }

  /** @brief Destructor */
  ~Interpolator();

  /** @brief Performs the interpolation at an n-point v
   *
   *  @param	v	Vector of coordinates of the n-point
   */
  double Evaluate(const std::vector<double>& v) const;

  /** @brief Performs the interpolation at an n-point v
   *
   *  @param	v	One dimensional value of the point
   */
  double Evaluate(const double& v) const;

  /** @brief Performs the interpolation at an n-point v
   *
   *  @param	v	Pointer to an array of coordinates of the n-point
   */
  double Evaluate(const double* v) const;

  /** @brief Returns the dimension of the space in which the interpolator lives */
  const size_t& GetDimension() const			{ return _dim; }

  /** @brief Returns the grid on which the interpolator is built */
  const Grid& GetGrid() const				{ return _grid; }

  /** @brief Returns the tesselation on which the interpolator is built */
  const Delaunay& GetTesselation() const		{ return _del; }

  /** @brief Returns the meshing, used to build the interpolation, as polygons
   *
   *  @return		Meshing
   */
  std::vector<TPolyLine*> Meshing() const;

 private:

  size_t		_dim;		///< Dimension of the space the interpolator lives in
  std::string		_algo;		///< Algorithm to use for the interpolation
  Grid			_grid;		///< Grid of N vertices and their mapping
  Delaunay		_del;		///< Delaunary triangulation of the N vertices
  bool			_ex;		///< Extrapolation switch
};

#endif
