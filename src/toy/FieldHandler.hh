// C++ includes
#include <iostream>
#include <ctime>
#include <cmath>
#include <map>

// ROOT includes
#include "TF1.h"
#include "TVector3.h"

// Additional includes
#include "Assert.hh"
#include "Pitch.hh"

/** @brief Data structure that stores the magnet type and its parameters
 */
struct FieldElement {
  std::string			magtype;	///< Type of magnet
  std::map<std::string, double> chars;		///< Field parameters
};

/** @brief Builds, stores and reads the fields present in the beam line
 *
 * 	   Builds theoretical magnetic fields based on the current in the coils,
 *	   returns the x, y and z components of the field at a requested point (x, y, z)
 */
class FieldHandler {
 public:
  /** @brief Default constructor, initilizes everything to 0 */
  FieldHandler();

  /** @brief Copy constructor */
  FieldHandler(const FieldHandler& fh);

  /** @brief Equality operator */
  FieldHandler& operator=(const FieldHandler& fh);

  /** @brief Destructor */
  ~FieldHandler();

  /** @brief Adds a solenoid to the dictionary of fields
   *
   *  @param	name	Name of the field/magnet
   *  @param	a	Radius of the solenoid
   *  @param	L	Length of the solenoid
   *  @param	I	Current in the solenoid
   *
   *  @return		True if successful
   */
  bool AddSolenoid(const std::string& name,
		   const double& a,
		   const double& L,
		   const double& I);

  /** @brief Gets the total magnetic field at a location x, y, z in cartesian coordinates
   *
   *  @param	x	Horizontal offset in which to evaluate the field
   *  @param	y	Vertical offset in which to evaluate the field
   *  @param	z	Position along the beam line in which to evaluate the field
   *
   *  @return		\f$\vec{B}\f$ in cartesian coordinates
   */
  TVector3 GetMagneticField(const double& x,
			    const double& y,
			    const double& z) const;

 private:

  /** @brief Cartesian 3-vector of the solenoidal magnetic field in x, y, z
   *
   *  @param	fel	Field element to get the field from
   *  @param	x	Horizontal offset in which to evaluate the field
   *  @param	y	Vertical offset in which to evaluate the field
   *  @param	z	Position along the beam line in which to evaluate the field
   *
   *  @return		\f$\vec{B}\f$ in cartesian coordinates
   */
  TVector3 GetSolenoidField(const FieldElement& fel,
			    const double& x,
			    const double& y,
			    const double& z) const;

  /** @brief Perpendicular component of a solenoidal field
   *
   *  @param	a	Radius of the solenoid
   *  @param	L	Length of the solenoid
   *  @param	I	Current in the solenoid
   *  @param	rho	Radius of the position in which to evaluate the field
   *  @param	z	Longitudinal position in which to evaluate the field
   *
   *  @return		 \f$B_\perp(\rho, 0, z)\f$ (cylindrical coords)
   */
  double SolenoidBperp(const double& a,
		       const double& L,
		       const double& I,
		       const double& rho,
		       const double& z) const;

  /** @brief Longitudinal component of a solenoidal field
   *
   *  @param	a	Radius of the solenoid
   *  @param	L	Length of the solenoid
   *  @param	I	Current in the solenoid
   *  @param	rho	Radius of the position in which to evaluate the field
   *  @param	z	Longitudinal position in which to evaluate the field
   *
   *  @return		\f$B_z(\rho, 0, z)\f$ (cylindrical coords)
   */
  double SolenoidBlong(const double& a,
		       const double& L,
		       const double& I,
		       const double& rho,
		       const double& z) const;

  /** @brief Elliptical integral of the first kind with parameter m */
  double K(const double& m) const;

  /** @brief Elliptical integral of the second kind with parameter m */
  double E(const double& m) const;

  /** @brief Elliptical integral of the third kind with parameter n and m */
  double Pi(const double& n, const double& m) const;

  double				_mu0;		///< Vacuum permeability
  std::map<std::string, FieldElement>	_fields;	///< List of fields in the beam line
};
