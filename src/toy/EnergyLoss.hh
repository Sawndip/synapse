// C++ includes
#include <iostream>
#include <ctime>
#include <cmath>
#include <vector>
#include <map>

// Additional includes
#include "Assert.hh"
#include "Pitch.hh"


/** @brief Data structure that stores the material properties.
 *
 * 	   Needs to be provided with the material characteristics.
 */
struct Material {
  std::string	name;	///< Name of the element/molecule/crystal
  double	Z;	///< Effective atomic number
  double	A;	///< Atomic mass [g/mol]
  double	rho;	///< Density [g/cm^3]
  double	I;	///< Mean excitation potential [eV]
  double	k1;	///< First factor of the energy loss [MeV/mm]
  double	k2;	///< Second factor of the energy loss
};


/** @brief Computes the ionisation energy loss of a particle through a specfied Material
 *
 * 	   Needs to be provided with the material characteristics and its thickness.
 */
class EnergyLoss {
 public:
  /** @brief Default constructor, initilizes everything to 0 */
  EnergyLoss();

  /** @brief Proper constructor, initilizes everything to the provided material
   *
   *  @param	mat		Properties of the absorber material
   *  @param	thickness	Thickness of the absorber [mm]
   */
  EnergyLoss(const Material& mat,
	     const double& thickness);

  /** @brief Copy constructor */
  EnergyLoss(const EnergyLoss& eloss);

  /** @brief Equality operator */
  EnergyLoss& operator=(const EnergyLoss& eloss);

  /** @brief Destructor */
  ~EnergyLoss();

  /** @brief Initialize the material for the energy loss */
  void Initialize();

  /** @brief Returns the mean momentum loss for a given particle
   *
   *  @param	p	Input momentum
   *  @param	m	Input particle mass
   *
   *  @return		Momentum loss [MeV/c]
   */
  double GetMomentumLoss(const double& p,
			 const double& m) const;

  /** @brief Returns the most probable momentum loss for a given particle
   *
   *  @param	p	Input momentum
   *  @param	m	Input particle mass
   *
   *  @return		Momentum loss [MeV/c]
   */
  double GetMostProbableMomentumLoss(const double& p,
			 	     const double& m) const;

  /** @brief Performs the energy loss on a beam (dictionary that maps variables onto vectors)
   *
   *  @param	beam	Input beam: x, (y), px, (py), pz [mm, MeV/c]
   *  @param	p	Total input beam momentum [MeV/c]
   *  @param	m	Particle mass [MeV/c^2]
   */
  void IonizeBunch(std::map<std::string, std::vector<double>>& beam,
		   const double& p,
		   const double& m) const;

  /** @brief Returns the Bethe-Bloch stopping formula
   *
   *  @param	bg	Relativistic &beta;&gamma; (p/m)
   *
   *  @return		Value of the formula in &beta;&gamma; [MeV/mm]
   */
  double BetheBloch(const double& bg) const;

  /** @brief Gets the thickness in amount of radiation length */
  const Material& GetMaterial()	const		{ return _mat; }

  /** @brief Sets the thickness in amount of radiation length */
  void SetStepSize(const double& stepsize)	{ _stepsize = stepsize; }

  /** @brief Gets the thickness in amount of radiation length */
  const double& GetStepSize()	const		{ return _stepsize; }

  /** @brief Sets the thickness in amount of radiation length */
  void SetThickness(const double& thickness)	{ _thickness = thickness; }

  /** @brief Gets the thickness in amount of radiation length */
  const double& GetThickness()	const		{ return _thickness; }

 private:

  Material _mat;	///< Description of the absorber material
  double _stepsize;	///< Size of a step through the absorber
  double _thickness;	///< Thickness of the absorber
};
