#ifndef DENSITYESTIMATOR_HH
#define DENSITYESTIMATOR_HH

// C++ includes
#include <iostream>
#include <ctime>
#include <cmath>
#include <vector>

// ROOT includes
#include "TRandom3.h"

// Additional includes
#include "Vector.hh"
#include "Matrix.hh"
#include "Geometry.hh"
#include "OptimalBinning.hh"
#include "KNearestNeighbours.hh"
#include "LocalReachability.hh"
#include "TDE.hh"
#include "PBATDE.hh"
#include "ProbabilityContour.hh"

/** @brief Computes the non parametric density estimator of a set of points.
 *
 * 	   Builds a non parametric density estimator based on a training set of points in
 *	   n dimensions, sampled from the true parent distribution to estimate.
 *
 *	   This class encompasses a broad array of non parametric methods to produce an estimate.
 */
class DensityEstimator {
 public:
  /** @brief Default constructor, initilizes everything to 0 */
  DensityEstimator();

  /** @brief Cloud constructor, sets vector of points on which to perform the density estimation
   *
   *  @param	points		N points in a _dim dimensional space
   *  @param 	algo		Defines the density estimation algorithm to be used
   *  @param 	ex		Extrapolate outside the outmost points if requested
   *  @param 	norm		Normalization factor for the density (if truncated distribution)
   */
  DensityEstimator(const std::vector<std::vector<double>>& points,
	       	   const std::string algo="bin",
	       	   const bool ex=true,
	       	   const double norm=1.);

  /** @brief Copy constructor */
  DensityEstimator(const DensityEstimator& de);

  /** @brief Equality operator */
  DensityEstimator& operator=(const DensityEstimator& de);

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
  ~DensityEstimator();

  /** @brief Initializes the requested density estimator */
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

  /** @brief Returns the volume of an alpha contour of the probability density function
   *
   *  @param	alpha		Integrated probability to achieve
   *  @param	method		Contour volume calculation method
   *  @param	draw		Creates a grid to be able to draw the contour if requested
   *  @param	lower		Lower bounds of the bounding box in each dimension
   *  @param	upper		Upper bounds of the bounding box in each dimension
   *  @param	nsamples	Number of samples to produce for the MC estimate
   *
   *  @return			Volume of the alpha contour
   */
  double ContourVolume(double alpha,
		       const std::string method="mc",
		       const bool draw=true,
		       std::vector<double> lower=std::vector<double>(),
		       std::vector<double> upper=std::vector<double>(),
		       const size_t nsamples=1e6);

  /** @brief Returns the uncertainty on the uncertainty on the reconstructed contour volume
   *
   *  @return			Contour volume uncertainty
   */
  double ContourVolumeError();

  /** @brief Returns segments that represent the 1D contour on a function
   *
   *  @param	idx	ID of the axis in which to plot the contour
   *
   *  @return		Set of lines representing the contours
   */
  TH1F* Contour1D(size_t idx=0) const			{ return _prob.Contour1D(idx); }

  /** @brief Returns the requested 2D contour on a function
   *
   *  @param	idx	ID of the first axis in which to plot the contour
   *  @param	idy	ID of the second axis in which to plot the contour
   *
   *  @return		Histogram at a certain contour level
   */
  TH2F* Contour2D(size_t idx=0, size_t idy=1) const	{ return _prob.Contour2D(idx, idy); }

  /** @brief Returns the requested 3D contour on a function
   *
   *  @param	idx	ID of the first axis in which to plot the contour
   *  @param	idy	ID of the second axis in which to plot the contour
   *  @param	idz	ID of the third axis in which to plot the contour
   *
   *  @return		Isosurface of a TF3 at a certain contour level
   */
  TH3F* Contour3D(size_t idx=0, size_t idy=1, size_t idz=2) const
							{ return _prob.Contour3D(idx, idy, idz); }

  /** @brief Returns the dimension of the space in which the estimator is constructed */
  const size_t& GetDimension() const			{ return _dim; }

  /** @brief Returns the points on which the DensityEstimator is constructed */
  const std::vector<std::vector<double>>& GetPoints() const
							{ return _points; }

  /** @brief Returns the CVT points of the space */
  const std::vector<std::vector<double>>& GetCVTpoints() const
							{ return _tde.GetPointArray(); }

  /** @brief Returns a grid built on the DensityEstimator, which is defined everywhere
   *
   *  @param	number	Number of steps between min and max values
   *  @param	lower	Lower bounds of the grid in each dimension
   *  @param	upper	Upper bounds of the grid in each dimension
   *
   *  @return		Grid evaluated with the density estimator
   */
  Grid GetGrid(const std::vector<size_t>& number,
	       const std::vector<double> lower=std::vector<double>(),
	       const std::vector<double> upper=std::vector<double>()) const;

  /** @brief Returns the meshing used to build the interpolation as polygons
   *
   *  @return		Meshing
   */
  std::vector<TPolyLine*> Meshing() const;

  /** @brief Returns the tesselation of the space as polygons
   *
   *  @return		Tesselation
   */
  std::vector<TPolyLine*> Tesselation() const;

  /** @brief Sets the name of the density estimator */
  void SetName(const std::string& name)			{ _name = name+"_"+_algo; }

  /** @brief Returns a 1D graph of the estimator in the requested axis
   *
   *  @param	xmin	Minimum of the range of the coordinate to plot
   *  @param	xmax	Maximum of the range of the coordinate to plot
   *  @param	idx	ID of the first coordinate to represent
   *  @param	x	Point that the projection intersects
   */
  TGraph* Graph(double xmin, double xmax,
		size_t idx=0, std::vector<double> x=std::vector<double>()) const;

  /** @brief Returns a 2D histogram of the estimator in the requested axis
   *
   *  @param	xmin	Minimum of the range of the first coordinate to plot
   *  @param	xmax	Maximum of the range of the first coordinate to plot
   *  @param	ymin	Minimum of the range of the second coordinate to plot
   *  @param	ymax	Maximum of the range of the second coordinate to plot
   *  @param	idx	ID of the first coordinate to represent
   *  @param	idy	ID of the second coordinate to represent
   *  @param	x	Point that the projection intersects
   */
  TH2F* Graph2D(double xmin, double xmax, double ymin, double ymax,
		size_t idx=0, size_t idy=1, std::vector<double> x=std::vector<double>()) const;

  /** @brief Draws the estimator on whichever TCanvas is currently being used (ROOT) 
   *
   *  @param	xmin	Minimum of the range of the first coordinate to plot
   *  @param	xmax	Maximum of the range of the first coordinate to plot
   *  @param	ymin	Minimum of the range of the second coordinate to plot
   *  @param	ymax	Maximum of the range of the second coordinate to plot
   *  @param	opt	Drawing options
   *  @param	idx	ID of the first coordinate to represent
   *  @param	idy	ID of the second coordinate to represent
   *  @param	x	Point that the projection intersects
   */
  void Draw(double xmin, double xmax, double ymin=0, double ymax=0, const std::string opt="",
	    int idx=0, int idy=1, std::vector<double> x=std::vector<double>()) const;

  /** @brief Creates a canvas, draws the estimator on it and save it to a file
   *
   *  @param	xmin	Minimum of the range of the first coordinate to plot
   *  @param	xmax	Maximum of the range of the first coordinate to plot
   *  @param	ymin	Minimum of the range of the second coordinate to plot
   *  @param	ymax	Maximum of the range of the second coordinate to plot
   *  @param	opt	Drawing options
   *  @param	idx	ID of the first coordinate to represent
   *  @param	idy	ID of the second coordinate to represent
   *  @param	x	Point that the projection intersects
   */
  void Print(double xmin, double xmax, double ymin=0, double ymax=0, const std::string opt="",
	     int idx=0, int idy=1, std::vector<double> x=std::vector<double>()) const;

 private:

  size_t				_dim;		///< Dimension of the space
  std::string				_name;		///< Name of the density estimator
  std::string				_algo;		///< Density estimation algorithm
  std::vector<std::vector<double>>	_points;	///< Array of input points
  OptimalBinning			_ob;		///< Optimal binning density estimator
  KNearestNeighbours    		_knn;		///< kNN density estimator
  LocalReachability		    	_lrd;		///< LRD density estimator
  TDE					_tde;		///< Tesselation density estimator
  PBATDE				_pbatde;	///< PBA tesselation density estimator
  ProbabilityContour			_prob;		///< Stores the last requested contour
  bool					_ex;		///< Extrapolation switch
  double				_norm;		///< Norm. factor (truncated parent)
};

#endif
