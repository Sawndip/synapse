// C++ includes
#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
#include <cfloat>
#include <cmath>

// ROOT includes
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TText.h"

// Other includes
#include "Assert.hh"
#include "Pitch.hh"
#include "GeometryHandler.hh"

/** @brief Environment of beam related objects and methods */
namespace Beam {

/** @brief Holds a set of beam line apertures.
 *
 *         The object stores and manages a set of apertures encountered in the beam line. It
 *	   contains definitions of methods associated with the apertures.
 */
class Aperture {
 public:
  /** @brief Default constructor, initilizes everything to 0 */
  Aperture();

  /** @brief Normal constructor, sets apertures 
   *
   *  @param	starts	Starting points of the apertures [mm]
   *  @param	ends	End points of the apertures [mm]
   *  @param	radii	Radii of the apertures [mm]
   *  @param	names	Names of the apertures
   */
  Aperture(std::vector<double> starts,
	   std::vector<double> ends,
	   std::vector<double> radii,
	   std::vector<std::string> names = std::vector<std::string>());

  /** @brief Copy constructor */
  Aperture(const Aperture& Aperture);

  /** @brief Equality operator */
  Aperture& operator=(const Aperture& Aperture);

  /** @brief Destructor */
  ~Aperture();

  /** @brief Initializes the apertures
   *
   *  @param	starts	Starting points of the apertures [mm]
   *  @param	ends	End points of the apertures [mm]
   *  @param	radii	Radii of the apertures [mm]
   *  @param	names	Names of the apertures
   */
  void Initialize(std::vector<double> starts = std::vector<double>(),
	   	  std::vector<double> ends = std::vector<double>(),
	   	  std::vector<double> radii = std::vector<double>(),
	   	  std::vector<std::string> names = std::vector<std::string>());

  /** @brief Sets the default MICE apertures 
   *
   *  @param	filename	Path to the MICE geometry file
   */
  void SetMICEDefault(const std::string& filename);

  /** @brief Adds a single aperture to the object
   *
   *  @param	start	Starting point of the aperture [mm]
   *  @param	end	End point of the aperture [mm]
   *  @param	radius	Radius of the aperture [mm]
   *  @param	name	Name of the aperture
   */
  void Add(double start, double end, double radius, std::string name = "");

  /** @brief Checks the aperture radius at a given position along the beam line
   *
   *  @param	z	Position along the beam line [mm]
   *
   *  @return		Radius of the aperture at z [mm]
   */
  const double& Radius(double z) const;

  /** @brief Checks if the particle is within the apertures
   *
   *  @param	x	Horizontal coordinate [mm]
   *  @param	y	Vertical coordinate [mm]
   *  @param	z	Position along the beam line [mm]
   *
   *  @return		True if within the aperture
   */
  bool IsIn(double x, double y, double z) const;

  /** @brief Draws the apertures on a canvas
   *
   *  @param	name	Name of the beam line
   *  @param	dz	Stepping size between points [mm]
   */
  void Draw(std::string name="", double dz=1) const;

 private:
  std::vector<double>		_apertures;	///< Summary vector used to fetch the smallest
  std::vector<double>		_starts;	///< Starting points vector of the apertures
  std::vector<double>		_ends;		///< End points vector of the apertures
  std::vector<double>		_radii;		///< Radii of the apertures
  std::vector<std::string>	_names;		///< Names of the apertures
};
} // namespace Beam
