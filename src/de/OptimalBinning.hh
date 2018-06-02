#ifndef OPTIMALBINNING_HH
#define OPTIMALBINNING_HH

// C++ includes
#include <iostream>
#include <ctime>
#include <cmath>
#include <vector>

// ROOT includes
#include "TROOT.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TTree.h"

// Additional includes
#include "Statistics.hh"
#include "Geometry.hh"
#include "Histogram.hh"
#include "Interpolator.hh"
#include "Assert.hh"

/** @brief Computes the Optimal Binning density estimator of a set of points
 *
 *  	   This method finds the optimal binning of the space by an array of different methods. 
 * 	   The amount of data points in each bin is directly related to the density in the bin.
 */
class OptimalBinning {
 public:
  /** @brief Default constructor, initilizes everything to 0 */
  OptimalBinning();

  /** @brief Cloud constructor, sets vector of points and their mapping for a point cloud
   *
   *  @param	points		N points in a _dim dimensional space
   *  @param	algo		Algorithm to be used for binning optimization
   *  @param	interp		Perform an interpolation on the histogram if requested
   *  @param 	ex		Extrapolate outside the outmost points if requested
   */
  OptimalBinning(const std::vector<std::vector<double>>& points,
		 const std::string algo="scott",
		 const bool interp=true,
	       	 const bool ex=true);

  /** @brief Copy constructor */
  OptimalBinning(const OptimalBinning& ob);

  /** @brief Equality operator */
  OptimalBinning& operator=(const OptimalBinning& ob);

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
  ~OptimalBinning();

  /** @brief Initializes the histogram with the optimal binning
   *
   *  @param	points	\f$N_1\times...\times N_n\f$ vector of points of n dimensions
   *  @param	algo	Algorithm to be used for binning optimization
   */
  void Initialize(const std::vector<std::vector<double>>& points,
		  const std::string& algo);

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

  /** @brief Returns the array of lower bounds in each dimension */
  const std::vector<double>& GetLowerBoundArray() const	{ return _hist.GetLowerBoundArray(); }

  /** @brief Returns the lower bound in dimension i */
  const double& GetLowerBound(const size_t i=0) const	{ return _hist.GetLowerBound(i); }

  /** @brief Returns the array of upper bounds in each dimension */
  const std::vector<double>& GetUpperBoundArray() const	{ return _hist.GetUpperBoundArray(); }

  /** @brief Returns the upper bound in dimension i */
  const double& GetUpperBound(const size_t i=0) const	{ return _hist.GetUpperBound(i); }

  /** @brief Returns the amount of bins in one projection */
  size_t GetNbins(const size_t i=0) const	{ return _hist.GetNbins(i); }

  /** @brief Returns the meshing, used to build the interpolation, as polygons */
  std::vector<TPolyLine*> Meshing() const;

 private:

  size_t		_dim;		///< Dimension of the space the estimator lives in
  Histogram		_hist;		///< Optimal histogram for density estimation
  Interpolator		_int;		///< Grid interpolator
  bool			_interp;	///< Apply grid interpolation if true
  bool			_ex;		///< Extrapolate if true
};

#endif
