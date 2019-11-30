#ifndef MST_HH
#define MST_HH

// C++ includes
#include <iostream>
#include <vector>
#include <map>

// ROOT includes
#include "TLine.h"
#include "TGraph.h"
#include "TPolyLine3D.h"
#include "TPolyMarker3D.h"
#include "TCanvas.h"
#include "TView.h"

// Additional includes
#include "Geometry.hh"

typedef std::pair<size_t, size_t> EdgeIndex;

/** @brief Computes the Minimum Spanning Tree (MST) of a set of points
 *
 *	   A Minimum Spanning Tree (MST) or minimum weight spanning tree is a subset
 *     of the edges of a connected, edge-weighted (un)directed graph that
 *     connects all the vertices together, without any cycles and with the
 *     minimum possible total edge weight.
 */
class MST {
 public:
  /** @brief Default constructor, initilizes everything to 0 */
  MST();

  /** @brief Normal constructor, sets the vector of points
   *
   *  @param	points		N points
   */
  MST(const std::vector<std::vector<double>>& points);

  /** @brief Copy constructor */
  MST(const MST& mst);

  /** @brief Equality operator */
  MST& operator=(const MST& mst);

  /** @brief Destructor */
  ~MST();

  /** @brief Find the MST of the N points */
  void Initialize();

  /** @brief Returns the dimension of the space in which the tree is constructed */
  const size_t& GetDimension() const                            { return _dim; }

  /** @brief Returns the array of points on which the tree is constructed */
  const std::vector<std::vector<double>>& GetPointArray() const { return _points; }

  /** @brief Returns a specific point */
  const std::vector<double>& GetPoint(const size_t i) const     { return _points[i]; }

  /** @brief Get the length of the full tree (sum of the lengths of individual edges) */
  const double& GetLength() const                               { return _length; }

  /** @brief Returns an array of lines to be drawn in ROOT (only 2D)
   *
   *  @param	ndc	Use the normalised device coordinate (NDC) system
   *
   *  @return		Array of lines
   */
  std::vector<TLine*> Lines(bool ndc=false) const;

  /** @brief Returns an array of 3d lines to be drawn in ROOT (only 3D)
   *
   *  @return		Parent volume that contains all the polyhedra
   */
  std::vector<TPolyLine3D*> Lines3D() const;

  /** @brief Draw the MST with the requested properties
   *
   *  @param	points	Draws the points if requested
   */
  void Draw(const bool points=false) const;

  /** @brief Paints the MST onto a TCanvas and save it
   *
   *  @param	name	Name of the canvas
   *  @param	points	Draws the points if requested
   *  @param	exts	Array of extensions of the files to produce
   */
  void Paint(const std::string name="",
             const bool points=false,
             const std::vector<std::string> exts=std::vector<std::string>({"pdf"})) const;

 private:

  /** @brief Returns distance between vertices i and j
   *
   *  @param	i	Index of the first vertex
   *  @param	j	Index of the second vertex
   *
   *  @return		Distance between the vertices
   */
  double Distance(size_t i, size_t j) const;

  size_t                              _dim;		  ///< Dimension of the space
  std::vector<std::vector<double>>    _points;  ///< Vector of N vertices
  std::vector<EdgeIndex>			        _edges;   ///< Array of edge indices
  std::vector<double>				          _lengths; ///< Array of edge lengths
  double					                    _length;  ///< Total length of the tree
};

#endif
