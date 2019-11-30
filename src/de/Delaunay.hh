#ifndef DELAUNAY_HH
#define DELAUNAY_HH

// C++ includes
#include <iostream>
#include <vector>
#include <map>

// ROOT includes
#include "TPolyLine.h"

// QHULL includes
#include "Qhull.h"
#include "QhullVertexSet.h"
#include "QhullFacetList.h"

// Additional includes
#include "Vertex.hh"
#include "Geometry.hh"

/** @brief Computes the Delaunay triangulation of a set of points.
 *
 *	   The Delaunay triangulation is obtained by computing the (n+1)-dimensional
 *	   convex hull of the points by adding the L2 vector norm as an additional coordinate.
 *	   This forms a paraboloid which projected bottom end corresponds to the triangulation.
 *     Qhull performs the heavy lifting, this is a front end class.
 */
class Delaunay {
 public:
  /** @brief Default constructor, initilizes everything to 0 */
  Delaunay();

  /** @brief Normal constructor, sets the vector of vertices
   *
   *  @param	vertices	Vector of vertices
   */
  Delaunay(const std::vector<Vertex>& vertices);

  /** @brief Normal constructor, sets the vector of points
   *
   *  @param	vertices	Vector of points
   */
  Delaunay(const std::vector<std::vector<double>>& vertices);

  /** @brief Copy constructor */
  Delaunay(const Delaunay& del);

  /** @brief Equality operator */
  Delaunay& operator=(const Delaunay& del);

  /** @brief Destructor */
  ~Delaunay();

  /** @brief Find the Delaunay triangulation of the N vertices */
  bool Initialize();

  /** @brief Returns the dimension of the space in which the tesselation is constructed */
  const size_t& GetDimension() const			{ return _dim; }

  /** @brief Returns the array of vertices on which the triangulation is constructed */
  const std::vector<Vertex>& GetVertexArray() const	{ return _vertices; }

  /** @brief Returns a specific vertex */
  const Vertex& GetVertex(const size_t i) const		{ return _vertices[i]; }

  /** @brief Finds the closest facet and returns its vertices */
  const std::vector<Vertex>& GetBestFacet(const std::vector<double>& v, bool* isin=NULL) const;

  /** @brief Returns an array of polygons to be drawn in ROOT (only 2D) */
  std::vector<TPolyLine*> Polygons() const;

 private:

  size_t 				_dim;		///< Dimension of the space
  std::vector<Vertex>			_vertices;	///< Vector of N vertices
  std::map<size_t, std::vector<Vertex>> _facets;	///< Map of facets (simplices)
  orgQhull::Qhull*			_qhull;		///< Pointer to qhull object, owns the memory
};

#endif
