#ifndef VORONOI_HH
#define VORONOI_HH

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

// QHULL includes
#include "PointCoordinates.h"
#include "RboxPoints.h"
#include "QhullError.h"
#include "Qhull.h"
#include "QhullQh.h"
#include "QhullFacet.h"
#include "QhullFacetList.h"
#include "QhullFacetSet.h"
#include "QhullLinkedList.h"
#include "QhullVertex.h"
#include "QhullSet.h"
#include "QhullVertexSet.h"

// NanoFLANN includes
#include <nanoflann.hpp>
#include "KDTreeVectorOfVectorsAdaptor.h"

// Additional includes
#include "Vector.hh"
#include "Matrix.hh"
#include "Hull.hh"

typedef KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<double>>>	KDTree;

/** @brief Computes the Voronoi tessellation of a set of points.
 *
 *	   The Voronoi tessellation is obtained by computing the (n+1)-dimensional convex
 *	   hull of the points by adding the L2 vector norm as an additional coordinate.
 *	   This forms a paraboloid which projected bottom end corresponds to the Delaunay
 *	   triangulation. The Voronoi tessellation is the dual graph of the Delaunay triangulation.
 *
 *         Qhull performs the heavy lifting, this is a front end class.
 */
class Voronoi {
 public:
  /** @brief Default constructor, initilizes everything to 0 */
  Voronoi();

  /** @brief Cloud constructor, sets the vector of points
   *
   *  @param	points		N points in a _dim dimensional space
   *  @param	bounded		Add 2^_dim box points to contain the tessellation
   *  @param	alpha		Penalised Centroidal Voronoi Tessellation parameter
   *  @param	eps		Precision required for the pure CVT (alpha = 1.)
   */
  Voronoi(const std::vector<std::vector<double>>& points,
	  const bool bounded=true,
	  const double alpha=0.,
	  const double eps=1e-6);

  /** @brief Copy constructor */
  Voronoi(const Voronoi& vor);

  /** @brief Equality operator */
  Voronoi& operator=(const Voronoi& vor);

  /** @brief Destructor */
  ~Voronoi();

  /** @brief Find the Voronoi tessellation of the N points */
  void Initialize(const double alpha, const double eps);

  /** @brief Find the Voronoi tessellation of the N points in 1D (special case) */
  void Initialize1D(const double alpha, const double eps);

  /** @brief Returns the array of points on which the triangulation is constructed */
  const size_t GetNpoints() const				{ return _points.size(); }

  /** @brief Returns the array of points on which the triangulation is constructed */
  const std::vector<std::vector<double>>& GetPointArray() const	{ return _points; }

  /** @brief Returns a specific vertex */
  const std::vector<double>& GetPoint(const size_t i) const	{ return _points[i]; }

  /** @brief Returns the array of cells corresponding to the array of vertice */
  const std::vector<Hull>& GetCellArray() const			{ return _cells; }

  /** @brief Returns a specific cell */
  const Hull& GetCell(const size_t i) const			{ return _cells[i]; }

  /** @brief Returns the dimension of the tesselated space */
  const size_t& GetDimension()	const				{ return _dim; }

  /** @brief Returns true if box points have been added */
  const bool& IsBounded() const					{ return _bounded; }

  /** @brief Finds the closest cell and returns its hull */
  size_t GetBestCellID(const std::vector<double>& v) const;

  /** @brief Finds the closest cell and returns its hull */
  const Hull& GetBestCell(const std::vector<double>& v) const;

  /** @brief Returns the volume of a simplex in the current space */
  double SimplexVolume(const std::vector<double>& point) const;

  /** @brief Returns a rectangle to be drawn in ROOT (only 2D) */
  TPolyLine* BoundingBox() const;

  /** @brief Returns an array of polygons to be drawn in ROOT (only 2D) */
  std::vector<TPolyLine*> Polygons() const;

  /** @brief Returns a histogram where each bin is a cell and its content is its density */
  TH2Poly* DensityProfile() const;

  /** @brief Draws the Voronoi tesselation on a canvas
   *
   *  @param	name		Will print the canvas as "VT_name.pdf"
   *  @param	color		Implements a color scale inverse proportional to the cell volume
   *  @param	points		Draws the points on which the tesselation is constructed
   *  @param	axes		Draws axes
   */
  void Draw(std::string name="",
	    bool color=true,
	    bool points=false,
	    bool axes=false) const;

 private:

  size_t 				_dim;		///< Dimension of the space
  std::vector<std::vector<double>>	_points;	///< Vector of N points 
  std::vector<Hull>			_cells;		///< Voronoi cells (subhulls)
  bool					_bounded;	///< Switch to bound the points
  std::vector<double>			_lower;		///< Lower bounds in each dimension
  std::vector<double>			_upper;		///< Upper bounds in each dimension
  KDTree*				_kdtree;	///< k-d Tree of the points to find NNs
};

#endif
