#ifndef HULL_HH
#define HULL_HH

// C++ includes
#include <iostream>
#include <ctime>
#include <cmath>
#include <vector>

// ROOT includes
#include "TPolyLine.h"
#include "TCanvas.h"

// Additional includes
#include "Matrix.hh"
#include "Geometry.hh"

/** @brief Data stucture of a convex hull.
 *
 *	   The object orgQhull::Qhull cannot be copied and hence makes it hard to memory handle.
 *	   This class stores the basic parameters computed by Qhull to be restituted on command.
 *
 *         Qhull performs the heavy lifting, this is a data stucture class.
 */
class Hull {
 public:
  /** @brief Default constructor, initilizes everything to 0 */
  Hull();

  /** @brief Normal constructor, sets the vector of points
   *
   *  @param	points		Vector of points
   */
  Hull(const std::vector<std::vector<double>>& points);

  /** @brief Normal constructor, sets the vector of points
   *
   *  @param	points		Vector of points in line
   *  @param	dim		Dimension of the space
   */
  Hull(const std::vector<double>& points,
       const size_t& dim);

  /** @brief Copy constructor */
  Hull(const Hull& hull);

  /** @brief Equality operator */
  Hull& operator=(const Hull& hull);

  /** @brief Destructor */
  ~Hull();

  /** @brief Find the Hull triangulation of the N points */
  void Initialize(const std::vector<std::vector<double>>& points);

  /** @brief Sets the array of points on the hull */
  void SetPointArray(const std::vector<std::vector<double>>& points)	{ this->Initialize(points); }

  /** @brief Sets the array of points on the hull */
  void SetPointArray(const std::vector<double>& points, const size_t& dim);

  /** @brief Sets the array of points on the hull */
  void AddPoint(const std::vector<double>& point)			{ _points.push_back(point); }

  /** @brief Returns the array of points on which the triangulation is constructed */
  const std::vector<std::vector<double>>& GetPointArray() const		{ return _points; }

  /** @brief Returns a specific point */
  const std::vector<double>& GetPoint(const size_t i) const		{ return _points[i]; }

  /** @brief Sets the volume of the hull */
  void SetVolume(const double& vol)					{ _vol = vol; }

  /** @brief Returns the volume of the hull */
  const double& GetVolume() const					{ return _vol; }

  /** @brief Returns the points ordered clockwise with respect to the CMF */
  Matrix<double> GetCWPoints() const;

  /** @brief Returns a polygon of the hull to be drawn in ROOT (only 2D) */
  TPolyLine* Polygon() const;


 private:

  size_t 			   	_dim;		///< Dimension of the space
  std::vector<std::vector<double>> 	_points;	///< N points on the hull
  double			   	_vol;		///< Volume of the convex hull
};

#endif
