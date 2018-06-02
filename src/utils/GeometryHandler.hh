#ifndef GEOMETRYHANDLER_HH
#define GEOMETRYHANDLER_HH

// Cpp includes
#include <iostream>
#include <string>
#include <fstream>
#include <iterator>
#include <map>

// ROOT includes
#include "TVector3.h"

// Other includes
#include "Exception.hh"

// Mapping between the names and geometry tags
static std::map<std::string, std::string> geo_mapping =
	{{"tof0","TOF0"}, {"tof1","TOF1"}, {"tof2","TOF2"},
	 {"ckova","Ckov1"}, {"ckovb","Ckov2"}, {"tku","Tracker0"}, {"tkd", "Tracker1"},
	 {"kl","KL"}, {"emr","EMR"}, {"abs","Disk_LiH LH2-Empty LH2-Full NeAbs"}};

// Tracker stations offset from centre point
static std::vector<double> station_offsets = {-550., -350., -100., 200., 550.};

/** @brief Stores the positions of a MICE detector in global coordinates.
  *
  * 	   This object extracts the globals of the requested module from the MICE .dat geometry
  *	   file. It retains its positions and rotations in the MICE global coordinate system.
 */
class DetectorGlobals {
 public:

  /** @brief Default constructor, initilizes everything to 0 */
  DetectorGlobals();

  /** @brief Import the requested detector 
   *
   *  @param	filename	Path to the MICE geometry file
   *  @param	name		Name of the detector
   */
  DetectorGlobals(const std::string& filename,
		  const std::string& name);

  /** @brief Copy constructor */
  DetectorGlobals(const DetectorGlobals& detector_globals);

  /** @brief Equality operator */
  DetectorGlobals& operator=(const DetectorGlobals& detector_globals);

  /** @brief Destructor */
  ~DetectorGlobals();

  /** @brief Overload the subscript operator, return a value */
  const TVector3& operator[](const size_t i) const	{ return _subpos[i]; }

  /** @brief Overload the subscript operator, return a value, allows mods */
  TVector3& operator[](const size_t i)			{ return _subpos[i]; }

  /** @brief Loads the global constants from the parent geometry file 
   *
   *  @param	filename	Path to the MICE geometry file
   *  @param	name		Name of the module
   */
  void Load(const std::string& filename,
	    const std::string& name);

  /** @brief Sets the global position of the detector */
  void SetPosition(const TVector3& pos);

  /** @brief Returns the global position of the detector */
  const TVector3& Position() const		{ return _pos; }

  /** @brief Sets the global rotation angles of the detector */
  void SetAngles(const TVector3& angles)	{ _angles = angles; }

  /** @brief Returns the global rotation angles of the detector */
  const TVector3& Angles() const 		{ return _angles; }

  /** @brief Returns the global x coordinate of the detector */
  double x() const				{ return _pos.x(); }

  /** @brief Returns the global y coordinate of the detector */
  double y() const				{ return _pos.y(); }

  /** @brief Returns the global z coordinate of the detector */
  double z() const				{ return _pos.z(); }

  /** @brief Returns the global pitch of the detector */
  double alpha() const				{ return _angles.x(); }

  /** @brief Returns the global yaw of the detector */
  double beta()	const				{ return _angles.y(); }

  /** @brief Returns the global roll of the detector */
  double gamma() const				{ return _angles.z(); }

 private:

  TVector3			_pos;				///< Position of module
  TVector3 			_angles;			///< Rotation angles of module
  std::vector<TVector3>		_subpos;			///< Positions of submodules
};

/** @brief Stores the positions of all the requested MICE detectors
  *
  * 	   This object extracts sets up a DetectorGlobals for each of the requested modules.
 */
class GeometryHandler {
 public:

  /** @brief Default constructor, initilizes everything to 0 */
  GeometryHandler();

  /** @brief Import all of the MICE detectors and elements
   *
   *  @param	filename	Path to the MICE geometry file
   */
  GeometryHandler(const std::string& filename);

  /** @brief Import only the list of requested detectors
   *
   *  @param	filename	Path to the MICE geometry file
   *  @param	names		List of MICE modules to look for
   */
  GeometryHandler(const std::string& filename,
		  const std::vector<std::string>& names);

  /** @brief Copy constructor */
  GeometryHandler(const GeometryHandler& geohandler);

  /** @brief Equality operator */
  GeometryHandler& operator=(const GeometryHandler& geohandler);

  /** @brief Destructor */
  ~GeometryHandler();

  /** @brief Overload the subscript operator, return a value */
  const DetectorGlobals& operator[](const std::string& name) const;

  /** @brief Overload the subscript operator, return a value, allows mods */
  DetectorGlobals& operator[](const std::string& name);

 private:

  std::vector<std::string>			_names;		///< Short names of modules
  std::map<std::string, DetectorGlobals>	_detectors;	///< List of detector globals
};

#endif
