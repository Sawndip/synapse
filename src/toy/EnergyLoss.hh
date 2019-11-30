#ifndef ENERGYLOSS_HH
#define ENERGYLOSS_HH

// C++ includes
#include <iostream>
#include <ctime>
#include <cmath>
#include <vector>
#include <map>

// ROOT includes
#include "TF1.h"
#include "TCanvas.h"

// Additional includes
#include "Assert.hh"
#include "Pitch.hh"
#include "MaterialDefinitions.hh"

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
   *  @param	mat		Properties of the material
   *  @param	thickness	Thickness of the material [mm]
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

  /** @brief Returns the mean momentum loss for a given particle
   *
   *  @param	p	Input momentum
   *  @param	m	Input particle mass
   *
   *  @return		Momentum loss [MeV/c]
   */
  double GetApproxMomentumLoss(const double& p,
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

  /** @brief Returns the total stopping power of the material
   *
   *  @param	bg	Relativistic &beta;&gamma; (p/m)
   *
   *  @return		Value of the stopping power in &beta;&gamma; [MeV/mm]
   */
  double StoppingPower(const double& bg) const;

  /** @brief Returns the Bethe-Bloch stopping formula
   *
   *  @param	bg	Relativistic &beta;&gamma; (p/m)
   *  @param	el	Element through which the particle is going
   *
   *  @return		Value of the formula in &beta;&gamma; [MeV/mm]
   */
  double BetheBloch(const double& bg,
		    const Element& el) const;

  /** @brief Returns the delta function, a fit to the density effect
   *
   *  @param	bg	Relativistic &beta;&gamma; (p/m)
   *  @param	el	Element through which the particle is going
   *
   *  @return		Value of the function in &beta;&gamma
   */
  double DeltaFunction(const double& bg,
		       const Element& el) const;

  /** @brief Sets the material definition */
  void SetMaterial(const Material& mat);

  /** @brief Gets the material definition */
  const Material& GetMaterial()	const		{ return _mat; }

  /** @brief Sets the thickness in amount of radiation length */
  void SetThickness(const double& thickness)	{ _thickness = thickness; }

  /** @brief Gets the thickness in amount of radiation length */
  const double& GetThickness() const		{ return _thickness; }

  /** @brief Draws and saves the Bethe-Bloch function as a function of &beta;&gamma; */
  void DrawBetheBloch() const;

 private:

  Material 		_mat;		///< Description of the materials
  double 		_thickness;	///< Thickness of the absorber
  std::vector<double> 	_params;	///< ROOT material parameters
  double		_Kstar;		///< Approximate energy loss prefactor
};

#endif
