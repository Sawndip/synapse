#ifndef PBATDE_HH
#define PBATDE_HH

// C++ includes
#include <iostream>
#include <ctime>
#include <cmath>
#include <vector>
#include <map>

// ROOT includes
#include "TStyle.h"
#include "TCanvas.h"
#include "TPolyLine.h"
#include "TH2Poly.h"
#include "TGraph.h"
#include "TCanvas.h"

// Additional includes
#include "Statistics.hh"
#include "TDE.hh"

/** @brief Computes the Penalised Bootstrap Aggregate Tesselation Density Estimator (PBATDE)
 *	   of a set of points.
 *
 *	   It is a density estimation method that divides the input sample of point
 *	   in a certain amount of subsamples (aggregates), applies TDE to each one of them,
 *	   then bootstraps the estimations into one, the mean estimate. The best number of
 *	   subsamples is achieved when the AIC of BIC is minimized (information criterion).
 */
class PBATDE {
 public:
  /** @brief Default constructor, initilizes everything to 0 */
  PBATDE();

  /** @brief Cloud constructor, sets the vector of points
   *
   *  @param	points		N points in a _dim dimensional space
   *  @param	interp		Interpolates between the vertices if requested
   *  @param	J		Number of points in each subsample
   */
  PBATDE(const std::vector<std::vector<double>>& points,
	 const bool interp=false,
	 const size_t J=100);

  /** @brief Copy constructor */
  PBATDE(const PBATDE& pbatde);

  /** @brief Equality operator */
  PBATDE& operator=(const PBATDE& pbatde);

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
  ~PBATDE();

  /** @brief Initialize the PBATDE in the N points */
  void Initialize();

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

  /** @brief Returns the array of points on which the triangulation is constructed */
  const std::vector<std::vector<double>>& GetPointArray() const	{ return _points; }

  /** @brief Returns a specific point */
  const std::vector<double>& GetPoint(const size_t i) const	{ return _points[i]; }

  /** @brief Returns the dimension of the tesselated space */
  const size_t& GetDimension()	const				{ return _dim; }

 private:

  size_t				_dim;		///< Dimension of the space the estimator lives in
  size_t				_M;		///< Number of subsamples to produce
  size_t				_J;		///< Number of points in each subsample
  std::vector<std::vector<double>>	_points;	///< Array of input points
  std::vector<TDE> 			_tdes;		///< Array of bootstrapped TDEs
  bool					_interp;	///< Apply simplexial interpolation if true
};

#endif
