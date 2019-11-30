#ifndef ALPHACOMPLEX_HH
#define ALPHACOMPLEX_HH

// C++ includes
#include <iostream>
#include <vector>

// ROOT includes
#include "TPolyLine.h"
#include "TGraph.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoMatrix.h"
#include "TGeoBBox.h"
#include "TGeoArb8.h"
#include "TView.h"
#include "TCanvas.h"
#include "TPolyMarker3D.h"

// QHULL includes
#include "Qhull.h"
#include "QhullVertexSet.h"
#include "QhullFacetList.h"

// Additional includes
#include "Vector.hh"
#include "Matrix.hh"
#include "Geometry.hh"

/** @brief Computes the &alpha;-complex of a set of points
 *
 *	   An &alpha;-complex is the most general definition of a concave hull. The parameter alpha
 *	   defines the curvature of the circumcircle that limits the maximum distance between
 *	   two points on the hull.
 */
class AlphaComplex {
 public:
  /** @brief Default constructor, initilizes everything to 0 */
  AlphaComplex();

  /** @brief Normal constructor, sets the vector of points
   *
   *  @param	points		N points
   *  @param	alpha		Curvature alpha
   */
  AlphaComplex(const std::vector<std::vector<double>>& points,
	       const double alpha=0.);

  /** @brief Copy constructor */
  AlphaComplex(const AlphaComplex& alc);

  /** @brief Equality operator */
  AlphaComplex& operator=(const AlphaComplex& alc);

  /** @brief Destructor */
  ~AlphaComplex();

  /** @brief Find the &alpha;-complex triangulation of the N points */
  void Initialize();

  /** @brief Returns the dimension of the space in which the tesselation is constructed */
  const size_t& GetDimension() const				{ return _dim; }

  /** @brief Returns the array of points on which the &alpha;-complex is constructed */
  const std::vector<std::vector<double>>& GetPointArray() const	{ return _points; }

  /** @brief Returns a specific point */
  const std::vector<double>& GetPoint(const size_t i) const	{ return _points[i]; }

  /** @brief Returns the volume of the &alpha;-complex */
  const double& GetVolume() const				{ return _vol; }

  /** @brief Returns the list of volumes of the cells in the &alpha;-complex */
  const std::vector<double>& GetCellVolumeArray() const		{ return _areas; }

  /** @brief Sets the curvature of the maximum circucircle radius; &alpha;=1/R */
  void SetAlpha(const double alpha);

  /** @brief Returns an array of polygons to be drawn in ROOT (only 2D)
   *
   *  @param	fill	Fills the cells and remove edges if requested
   *
   *  @return		Array of polygons
   */
  std::vector<TPolyLine*> Polygons(const bool fill=false) const;

  /** @brief Returns a parent volume that contains all of the simplices (only 3D)
   *
   *  @return		Parent volume that contains all the polyhedra
   */
  TGeoVolume* Polyhedra() const;

  /** @brief Draw the &alpha;-complex with the requested properties
   *
   *  @param	fill	Fills the cells and remove edges if requested
   *  @param	points	Draws the points if requested
   */
  void Draw(const bool fill=false,
	    const bool points=false) const;

  /** @brief Paints the &alpha;-complex onto a TCanvas and save it
   *
   *  @param	name	Name of the canvas
   *  @param	fill	Fills the cells and remove edges if requested
   *  @param	points	Draws the points if requested
   *  @param	exts	Array of extensions of the files to produce
   */
  void Paint(const std::string name="",
	     const bool fill=false,
	     const bool points=false,
	     const std::vector<std::string> exts=std::vector<std::string>({"pdf"})) const;

 private:

  /** @brief Returns the circumradius of a simplical facet
   *
   *  @param	facet	Facet to compute the circumradius of
   */
  double Circumradius(orgQhull::QhullFacet& facet) const;

  /** @brief Returns a TGeoVolume discribing the simplex in 3D
   *
   *  @param	facet	Facet to represent
   *  @param	H	Transformation matrix to place the facet
   */
  TGeoVolume* Polyhedron(const std::vector<std::vector<double>>& facet,
			 TGeoHMatrix*& H) const;

  size_t 					_dim;		///< Dimension of the space
  double					_alpha;		///< Curvature alpha
  std::vector<std::vector<double>>		_points;	///< Vector of N vertices
  std::vector<std::vector<std::vector<double>>>	_facets;	///< Array of simplices
  std::vector<double>				_circumradii;	///< Array of facet circumradii
  std::vector<double>				_areas;		///< Array of facet areas
  double					_vol;		///< Volume of the complex
};

#endif
