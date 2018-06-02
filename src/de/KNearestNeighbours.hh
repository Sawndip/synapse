#ifndef KNEARESTNEIGHBOURS_HH
#define KNEARESTNEIGHBOURS_HH

// C++ includes
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

// ROOT includes
#include "TMatrixD.h"

// NanoFLANN includes
#include <nanoflann.hpp>
#include "KDTreeVectorOfVectorsAdaptor.h"

// Additional includes
#include "Assert.hh"
#include "Statistics.hh"

typedef KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<double>>>	KDTree;

/** @brief Computes the k-Nearest Neighbour (kNN) density estimator of a set of points
 *
 *	   Computes the volume and density of the hypersphere containing the k nearest neighbours.
 */
class KNearestNeighbours {
 public:
  /** @brief Default constructor, initilizes everything to 0 */
  KNearestNeighbours();

  /** @brief Normal constructor, fills the cloud with the provided points
   *  @param	points		N points to fill the cloud with
   *  @param	k		Number of neighbours to probe
   *  @param	rotate		If true, use the metric of the covariance matrix
   */
  KNearestNeighbours(const std::vector<std::vector<double>>& points,
		     const size_t k,
		     bool rotate=false);

  /** @brief Copy constructor */
  KNearestNeighbours(const KNearestNeighbours& knn);

  /** @brief Equality operator */
  KNearestNeighbours& operator=(const KNearestNeighbours& knn);

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
  ~KNearestNeighbours();

  /** @brief Checks that the input is sensible */
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

  /** @brief Returns the dimension of the space */
  const size_t& GetDimension() const			{ return _dim; }

  /** @brief Returns the distance to the k^th NN of a point */
  double Distance(const std::vector<double>& point) const;

  /** @brief Returns the volume of the hypersphere that contains the k NN of a point */
  double Volume(const std::vector<double>& point) const;

  /** @brief Returns the density of points in the hypersphere that contains the k NN of a point */
  double Density(const std::vector<double>& point) const;

 private:

  size_t				_dim;		///< Dimension of the space
  std::vector<std::vector<double>>	_points;	///< Cloud of points
  size_t				_k;		///< Number of neighbours to consider
  KDTree*				_kdtree;	///< k-d Tree of the points to find NNs
  bool					_rotate;	///< Whether or not a rotated metric is used
  TMatrixD				_metric;	///< Metric of the points
  double				_scale;		///< Scaling factor of the metric
};

#endif
