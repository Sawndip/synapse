#ifndef PROBABILITYCONTOUR_HH
#define PROBABILITYCONTOUR_HH

// C++ includes
#include <iostream>
#include <ctime>
#include <cmath>
#include <vector>
#include <algorithm>

// ROOT includes
#include "TString.h"
#include "TLine.h"
#include "TPolyLine.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TMath.h"

// Other includes
#include "Grid.hh"
#include "AlphaComplex.hh"

/** @brief Reconstructs the optimal &alpha;-contour of a Probability Density Function 
 *
 *	   Provided with a scalar field normalised to 1 and defined everywhere, the
 *	   algorithm finds the manifold that bounds the &alpha;-contour and reconstructs its
 *	   volume. The &alpha;-contour is the smallest hypersurface that contains a fraction
 *	   &alpha; of the total probability.
 */
class ProbabilityContour {
 public:
  /** @brief Default constructor, initilizes everything to 0 */
  ProbabilityContour();

  /** @brief Hull constructor, sets the vertices and their mapping, computes hull volumes
   *
   *  @param	vertices	Density estimation input vertices
   *  @param	grid		Grid on which the distrubution is evaluated (to draw only)
   *  @param	algo	Probability contour estimation algorithm
   */
  ProbabilityContour(const std::vector<Vertex>& vertices,
		     const Grid grid=Grid(),
		     const std::string algo = "hull");

  /** @brief Grid constructor, sets a grid of points at which the function is known
   *
   *  @param	grid	Grid on which the distrubution is evaluated
   *  @param	algo	Probability contour estimation algorithm
   */
  ProbabilityContour(const Grid& grid,
		     const std::string algo = "rectangle");

  /** @brief MC constructor, sets the vertices, mcvertices and box volume 
   *
   *  @param	vertices	Density estimation input vertices
   *  @param	mcvertices	MC vertices uniformly distributed inside a box
   *  @param	boxvol		Scaled volume of the box in which the MC vertices were produced
   *  @param	grid		Grid on which the distrubution is evaluated (to draw only)
   */
  ProbabilityContour(const std::vector<Vertex>& vertices,
		     const std::vector<Vertex>& mcvertices,
		     const double& boxvol,
		     const Grid grid=Grid());

  /** @brief Copy constructor */
  ProbabilityContour(const ProbabilityContour& cont);

  /** @brief Equality operator */
  ProbabilityContour& operator=(const ProbabilityContour& cont);

  /** @brief Destructor */
  ~ProbabilityContour();

  /** @brief Returns segments that represent the 1D contour on a function
   *
   *  @param	idx	ID of the axis in which to plot the contour
   *
   *  @return		Set of lines representing the contours
   */
  TH1F* Contour1D(size_t idx=0) const;

  /** @brief Returns the requested 2D contour on a function
   *
   *  @param	idx	ID of the first axis in which to plot the contour
   *  @param	idy	ID of the second axis in which to plot the contour
   *
   *  @return		Histogram at a certain contour level
   */
  TH2F* Contour2D(size_t idx=0, size_t idy=1) const;

  /** @brief Returns the requested 2D contour on a function
   *
   *  @param	idx	ID of the first axis in which to plot the contour
   *  @param	idy	ID of the second axis in which to plot the contour
   *  @param	idz	ID of the third axis in which to plot the contour
   *
   *  @return		Histogram at a certain contour level
   */
  TH3F* Contour3D(size_t idx=0, size_t idy=1, size_t idz=2) const;

  /** @brief Computes the optimal contour using the specified algorithm
   *
   *  @param	alpha	Integrated probability to achieve
   *  @param	eps	Resolution with which to achieve the contour
   *
   *  @return		True if successful
   */
  void Evaluate(double alpha, double eps = 1e-3);

  /** @brief Returns the array of indices of the vertices used */
  bool IsInitialized() const					{ return _vertices.size(); }

  /** @brief Returns the fraction of the function for which the contour was computed */
  const double& GetAlpha() const				{ return _alpha; }

  /** @brief Returns the volume of the box containing the MC vertices */
  const double& GetBoxVolume() const				{ return _boxvol; }

  /** @brief Returns the volume of the current contour */
  const double& GetContourVolume() const			{ return _volume; }

  /** @brief Returns the volume of the current contour */
  double GetContourVolumeError() const;

  /** @brief Returns the array of indices of grid cells used to build the contour */
  const std::vector<std::vector<size_t>>& GetIndexArray() const	{ return _indices; }

  /** @brief Returns the array of MC vertices used to determine the contour */
  const std::vector<Vertex>& GetMCVertexArray() const		{ return _mcvertices; }

  /** @brief Returns the array of input vertices */
  const std::vector<Vertex>& GetVertexArray() const		{ return _vertices; }

  /** @brief Returns the volume of the contour */
  const double& GetVolume() const				{ return _volume; }

  /** @brief Sets the array of MC vertices used to determine the contour */
  void SetMCVertexArray(const std::vector<Vertex>& mcvertices)	{ _mcvertices = mcvertices; }

  /** @brief Sets the array of input vertices */
  void SetVertexArray(const std::vector<Vertex>& vertices)	{ _vertices = vertices; }

  /** @brief Sets the array of input vertices */
  void SetBoxVolume(const double& boxvol)			{ _boxvol = boxvol; }

 private:

  /** @brief Computes the optimal contour using the "alphac" method
   *
   *  @param	alpha	Integrated probability to achieve
   *  @param	eps	Resolution with which to achieve the contour
   */
  void AlphaComplexMethod(double alpha, double eps);

  /** @brief Computes the optimal contour using the "hull" method
   *
   *  @param	alpha	Integrated probability to achieve
   *  @param	eps	Resolution with which to achieve the contour
   */
  void HullMethod(double alpha, double eps);

  /** @brief Computes the optimal contour using the "mc" method
   *
   *  @param	alpha	Integrated probability to achieve
   *  @param	eps	Resolution with which to achieve the contour
   */
  void MCMethod(double alpha, double eps);

  /** @brief Computes the optimal contour using the "rectangle" method
   *
   *  @param	alpha	Integrated probability to achieve
   *  @param	eps	Resolution with which to achieve the contour
   */
  void RectangleMethod(double alpha, double);

  /** @brief Computes the optimal contour using the "trapezoid" method
   *
   *  @param	alpha	Integrated probability to achieve
   *  @param	eps	Resolution with which to achieve the contour
   */
  void TrapezoidMethod(double alpha, double eps = 1e-3);


  std::string				_algo;		///< Contour computation algorithm
  Grid					_grid;		///< Grid of points are their value
  std::vector<std::vector<size_t>>	_indices; 	///< Vector of indices of the vertices used
  double				_alpha;		///< Fraction of the function
  double				_volume;  	///< Volume of the contour
  double				_level;		///< Level of the contour
  std::vector<Vertex>			_vertices;	///< Density estimation input vertices
  std::vector<Vertex>			_mcvertices;	///< Monte Carlo vertices on a box
  double				_boxvol;	///< Volume of the box in which the MC is run
};

#endif
