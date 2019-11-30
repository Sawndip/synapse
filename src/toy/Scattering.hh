#ifndef SCATTERING_HH
#define SCATTERING_HH

// C++ includes
#include <iostream>
#include <ctime>
#include <cmath>
#include <map>

// ROOT includes
#include "TRandom3.h"

// Additional includes
#include "Assert.hh"

/** @brief Randomly scatters a particle across a specified material
 *
 *	   It uses the Gaussian approximation to scattering, with a width entirely defined
 * 	   by the material thickness in number \f$X_0\f$.
 *
 * 	   It is function of the input particle momentum and velocity.
 */
class Scattering {
 public:
  /** @brief Default constructor, initilizes everything to 0 */
  Scattering();

  /** @brief Normal constructor, sets the thickness and the scattering factor 
   *
   *  @param	thickness	Thickness of the absorber in \f$X_0\f$
   */
  Scattering(const double& thickness);

  /** @brief Copy constructor */
  Scattering(const Scattering& scat);

  /** @brief Equality operator */
  Scattering& operator=(const Scattering& scat);

  /** @brief Destructor */
  ~Scattering();

  /** @brief Performs the scattering on a particle 
   *
   *  @param	p	Total input momentum
   *  @param	m	Particle mass
   *
   *  @return		Scattering angle
   */
  double Scatter(const double& p,
		 const double& m);

  /** @brief Performs the scattering on a beam of particles 
   *
   *  @param	beam	Input beam: x, (y), px, (py), pz [mm, MeV/c]
   *  @param	p	Total input beam momentum [MeV/c]
   *  @param	m	Particle mass [MeV/c^2]
   */
  void ScatterBunch(std::map<std::string, std::vector<double>>& beam,
		    const double& p,
		    const double& m);

  /** @brief Sets the seed */
  void SetSeed(const size_t& seed)		{ _rdmzer = TRandom3(seed); }

  /** @brief Sets the thickness in amount of radiation length */
  void SetThickness(const double& thickness)	{ _thickness = thickness; }

  /** @brief Gets the thickness in amount of radiation length */
  const double& GetThickness()	const		{ return _thickness; }

  /** @brief Sets the scattering factor [MeV/c] */
  void SetFactor(const double& factor)		{ _factor = factor; }

  /** @brief Gets the scattering factor [MeV/c] */
  const double& GetFactor() const		{ return _factor; }

 private:

  /** @brief Compute the scattering factor from the thickness [MeV/c] */
  void ComputeFactor(const double& thickness);

  TRandom3 	_rdmzer;	///< Pseudorandom number generator
  double 	_thickness;	///< Thickness of the material in amount of \f$X_0\f$
  double 	_factor;	///< Prefactor of the scattering formul
};

#endif
