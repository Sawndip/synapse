#ifndef TDE_HH
#define TDE_HH

// C++ includes
#include <iostream>
#include <ctime>
#include <cmath>
#include <vector>
#include <map>

// Additional includes
#include "Voronoi.hh"
#include "Interpolator.hh"
#include "Assert.hh"

/** @brief Computes the Tesselation Density Estimator (TDE) of a set of points.
 *
 *	   The TDE uses the inverse of the volume of each surrounding Voronoi cell as an
 *	   estimator for density. The estimator supports simplex interpolation,
 *	   CVT and PCVT (as defined in the Voronoi class).
 */
class TDE {
 public:
  /** @brief Default constructor, initilizes everything to 0 */
  TDE();

  /** @brief Cloud constructor, sets the vector of points
   *
   *  @param	points		N points in a _dim dimensional space
   *  @param	interp		Interpolates between the vertices if requested
   *  @param	bounded		Add 2^_dim box vertices to contain the tessellation if requested
   *  @param	alpha		Penalised Centroidal TDE Tessellation parameter
   *  @param	eps		Precision required for the pure CVT (alpha = 1.)
   */
  TDE(std::vector<std::vector<double>> points,
      const bool interp=true,
      const bool bounded=true,
      const double alpha=0.,
      const double eps=1e-6);

  /** @brief Copy constructor */
  TDE(const TDE& tde);

  /** @brief Equality operator */
  TDE& operator=(const TDE& tde);

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
  ~TDE();

  /** @brief Performs the density estimation at an n-point v
   *
   *  @param	v	Vector of coordinates of the n-point
   */
  double Evaluate(const std::vector<double>& v) const;

  /** @brief Performs the density estimation at an n-point v
   *
   *  @param	v	One dimensional value of the point
   */
  double Evaluate(const double& v) const;

  /** @brief Performs the density estimation at an n-point v
   *
   *  @param	v	Pointer to an array of coordinates of the n-point
   */
  double Evaluate(const double* v) const;

  /** @brief Initializes the interpolator if requested */
  void Initialize();

  /** @brief Returns the array of points on which the triangulation is constructed */
  const std::vector<std::vector<double>>& GetPointArray() const	{ return _vor.GetPointArray(); }

  /** @brief Returns a specific vertex */
  const std::vector<double>& GetVertex(const size_t i) const	{ return _vor.GetPoint(i); }

  /** @brief Returns the array of cells corresponding to the array of points */
  const std::vector<Hull>& GetCellArray() const		 	{ return _vor.GetCellArray(); }

  /** @brief Returns a specific cell */
  const Hull& GetCell(const size_t i) const			{ return _vor.GetCell(i); }

  /** @brief Returns the dimension of the tesselated space */
  const size_t& GetDimension() const				{ return _vor.GetDimension(); }

  /** @brief Returns true if box points have been added */
  const bool& IsBounded() const					{ return _vor.IsBounded(); }

  /** @brief Returns a rectangle to be drawn in ROOT (only 2D) */
  TPolyLine* BoundingBox() const				{ return _vor.BoundingBox(); }

  /** @brief Returns an array of polygons to be drawn in ROOT (only 2D) */
  std::vector<TPolyLine*> Polygons() const		 	{ return _vor.Polygons(); }

  /** @brief Returns a histogram where each bin is a cell and its content is its density */
  TH2Poly* DensityProfile() const			 	{ return _vor.DensityProfile(); }

  /** @brief Draws the underlying Voronoi tesselation on a canvas
   *  @param	name		Will print the canvas as "VT_name.pdf"
   *  @param	color		Implements a color scale inverse proportional to the cell volume
   *  @param	points		Draws the points on which the tesselation is constructed
   *  @param	axes		Draws axes
   */
  void Draw(std::string name="",
	    bool color=true,
	    bool points=false,
	    bool axes=false) const 			{ _vor.Draw(name, color, points, axes); }

  /** @brief Returns the meshing, used to build the interpolation, as polygons */
  std::vector<TPolyLine*> Meshing() const;

 private:

  Voronoi 		_vor;		///< Voronoi tesselation of the space
  Interpolator		_int;		///< Simplexial interpolator
  std::vector<double>	_weights;	///< Weights of each of the point, if repetition
  bool			_interp;	///< Apply simplexial interpolation if true
};

#endif
