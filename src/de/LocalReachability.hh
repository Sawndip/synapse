#ifndef LOCALREACHABILITY_HH
#define LOCALREACHABILITY_HH

// C++ includes
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

// NanoFLANN includes
#include <nanoflann.hpp>
#include "KDTreeVectorOfVectorsAdaptor.h"

// Additional includes
#include "Assert.hh"

typedef KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<double>>>	KDTree;

/** @brief Computes the Local Reachability density estimator of a set of points
 *
 * 	   Computes reachability of the of the sample point by the cloud and
 *	   its corresponding density
 */
class LocalReachability {
 public:
  /** @brief Default constructor, initilizes everything to 0 */
  LocalReachability();

  /** @brief Normal constructor, fills the cloud with the provided points
   *
   *  @param	points	Vector of points
   *  @param	k	Number of neighbours to probe
   */
  LocalReachability(const std::vector<std::vector<double>>& points, const size_t k);

  /** @brief Copy constructor */
  LocalReachability(const LocalReachability& lr);

  /** @brief Equality operator */
  LocalReachability& operator=(const LocalReachability& lr);

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
  ~LocalReachability();

  /** @brief Checks that the input is sensible */
  void Initialize();

  /** @brief Performs the density estimation at an n-point v
   *
   *  @param	v	Vector of coordinates of the n-point
   */
  double Evaluate(const std::vector<double>& v) const	{ return RDensity(v); };

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
  const size_t& GetDimension() const				{ return _dim; }

  /** @brief Returns the the distance between a point A and a point B */
  double Distance(const std::vector<double>& A, const std::vector<double>& B) const;

  /** @brief Returns the distance to the k^th NN of a point A */
  double KDistance(const std::vector<double>& A) const;

  /** @brief Returns the the reachability distance of a point A by a point B */
  double RDistance(const std::vector<double>& A, const std::vector<double>& B) const;

  /** @brief Returns the volume of the hypersphere that contains of reachability radius */
  double RVolume(const double& radius) const;

  /** @brief Returns the the local reachability density of a point A */
  double RDensity(const std::vector<double>& A) const;

 private:

  size_t				_dim;		///< Dimension of the space
  std::vector<std::vector<double>>	_points;	///< Cloud of points
  size_t				_k;		///< Number of neighbours to consider
  KDTree*				_kdtree;	///< k-d Tree of the points to find NNs
};

#endif
