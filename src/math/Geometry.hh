#ifndef GEOMETRY_HH
#define GEOMETRY_HH

// C++ includes
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>

// ROOT includes
#include "TF1.h"

// Other includes
#include "Vector.hh"
#include "Matrix.hh"
#include "Statistics.hh"

/** @file Geometry.hh
 *
 *  @brief Set of functions that compose the Geometry package.
 *
 *	   Contains a broad array of useful mathematical methods in the context
 *         of computational geometry.
 */

/** @brief Mathematics environment (includes Statistics and Geometry packages)
 */
namespace Math {

/** @brief Computes the equation parameters of a plane given three points
 *
 *  @param	v0, v1, v2	Coordinates of the 3 points
 *
 *  @return			Array of parameters n_x, n_y, n_z, r_0
 */
std::vector<double> PlaneParameters(const Vector<double>& v0,
				    const Vector<double>& v1,
				    const Vector<double>& v2);

/** @brief Computes the equation parameters of a plane given three points
 *
 *  @param	points		Vector of coordinates of the 3 points
 *
 *  @return			Array of parameters n_x, n_y, n_z, r_0
 */
std::vector<double> PlaneParameters(const std::vector<Vector<double>>& points);

/** @brief Computes the distance between a plane of given equation parameters and a point
 *
 *  @param	pars		Array of parameters n_x, n_y, n_z, r_0
 *  @param	point		Coordinates of the point
 *
 *  @return			Distance
 */
double PlaneDistance(const std::vector<double>& pars,
                     const std::vector<double>& point);

/** @brief Returns true if the a lies to the left of b in the 2d plane
 *
 *  @param	a		First vector
 *  @param	b		Second
 *
 *  @return			True if to a to the left of b
 */
bool LeftHanded(const Vector<double>& a,
		const Vector<double>& b);

/** @brief Returns true if the a lies to the right of b in the 2d plane
 *
 *  @param	a		First vector
 *  @param	b		Second
 *
 *  @return			True if to a to the right of b
 */
bool RightHanded(const Vector<double>& a,
		 const Vector<double>& b);

/** @brief Sorts the input points clockwise in the the 2D plane
 *
 *  @param	points		Array of points to sort
 */
void SortCW(std::vector<Vector<double>>& points);

/** @brief Sorts the input points anticlockwise in the the 2D plane
 *
 *  @param	points		Array of points to sort
 */
void SortACW(std::vector<Vector<double>>& points);

/** @brief Builds a bounding box around a set of points with a certain buff
 *
 *  @param	points		Array of points to bound
 *  @param	lower		Array of lower bounds to set
 *  @param	upper		Array of upper bounds to set
 *  @param	buff		Fraction of the range in the buffer space
 *
 *  @return			Volume of the bounding box
 */
double BoundingBox(const std::vector<std::vector<double>>& points,
		   Vector<double>& lower,
		   Vector<double>& upper,
		   double buff=0.);

/** @brief Returns the n+1 vertices of an n-simplex
 *
 *  @param	n	Dimension of the space
 *
 *  @return		Array of n+1 points with n coordinates
 */
std::vector<std::vector<double>> SimplexVertices(const size_t n);

/** @brief Returns true if a point is inside a provided simplex
 *
 *  @param	point	Test point
 *  @param	simplex Test simplex
 *
 *  @return		True if the point is inside the simplex
 */
bool IsInsideSimplex(const std::vector<double>& point,
		     const std::vector<std::vector<double>>& simplex);

/** @brief Returns the lp-norm of a vector to the power of p
 *
 *
 *  @param	v	Vector
 *  @param	p	Parameter of the Lebesgue norm
 *
 *  @return		||v||_p^p
 */
double LpNormp(const std::vector<double>& v, const double& p);

/** @brief Returns the lp-norm of a vector
 *
 *  @param	v	Vector
 *  @param	p	Parameter of the norm
 *
 *  @return		||v||_p
 */
double LpNorm(const std::vector<double>& v, const double& p);

/** @brief Returns the volume of the unit d-ball
 *
 *  @param	d	Dimension
 *  @param	p	Lebesgue Metric

 *  @return		Volume of B_d^p
 */
double UnitBallVolume(const size_t& d, const double p=2.);

/** @brief Returns the hull volume bias factor for N points uniformly distributed inside a ball
 *
 *  @param	d	Dimension of the space
 *  @param	p	Parameter of the Lebesgue norm
 *
 *  @return		Bias factor
 */
double HullVolumeUniformFactor(const double& d, const double& p);

/** @brief Returns the relative bias (i.e. E(^V-V)/V) on the volume of the convexe hull of 
 *         N points uniformly distributed inside a p-unit ball in d dimensions
 *
 *  @param	d	Dimension of the space
 *  @param	p	Parameter of the Lebesgue norm
 *  @param	N	Number of samples
 *
 *  @return		Bias factor (1+beta)
 */
double HullVolumeUniformRelativeBias(const double& d, const double& p, const size_t N);

/** @brief Returns the relative bias (i.e. E(^V-V)/V) on the volume of the convexe hull of 
 *         of a fraction alpha of N points selected from a multivariate Gaussian in dimension d
 *
 *  @param	d	Dimension of the space
 *  @param	alpha	Fraction of the total number of points
 *  @param	N	Total number of samples
 *
 *  @return		Bias factor (1+beta)
 */
double HullVolumeTGausRelativeBias(const double& d, const double& alpha, const size_t N);

/** @brief Returns the relative RMS (i.e. RMS(^V)/^V) on the volume of the convexe hull of
 *         N points uniformly distributed inside a p-unit ball in d dimensions
 *
 *  @param	d	Dimension of the space
 *  @param	p	Parameter of the Lebesgue norm
 *  @param	N	Number of samples
 *
 *  @return		Relative RMS
 */
double HullVolumeUniformRelativeRMS(const double& d, const double& p, const size_t N);

/** @brief Returns the relative RMS (i.e. RMS(^V)/^V) on the volume of the convexe hull of
 *         of a fraction alpha of N points selected from a multivariate Gaussian in dimension d
 *
 *  @param	d	Dimension of the space
 *  @param	alpha	Fraction of the total number of points
 *  @param	N	Total number of samples
 *
 *  @return		Relative RMS
 */
double HullVolumeTGausRelativeRMS(const double& d, const double& alpha, const size_t N);
} // namespace Geometry

#endif
