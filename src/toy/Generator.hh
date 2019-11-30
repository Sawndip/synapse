#ifndef GENERATOR_HH
#define GENERATOR_HH

// C++ includes
#include <iostream>
#include <vector>
#include <map>
#include <ctime>

// ROOT includes
#include "TRandom3.h"

// Additional includes
#include "Assert.hh"
#include "Pitch.hh"
#include "DGaus.hh"
#include "DSpiral.hh"
#include "Bunch.hh"

/** @brief Computes the ionisation energy loss of a particle through a specfied Material
 *
 * 	   Needs to be provided with the material characteristics and its thickness.
 */
class Generator {
 public:
  /** @brief Default constructor, initilizes everything to 0 */
  Generator();

  /** @brief Proper constructor, initilizes everything to the provided material
   *
   *  @param	dim		Dimension of the transverse phase space
   *  @param	mass		Mass of the particle to generate [MeV/c^2]
   */
  Generator(const size_t dim,
	    const double& mass);

  /** @brief Copy constructor */
  Generator(const Generator& gen);

  /** @brief Equality operator */
  Generator& operator=(const Generator& gen);

  /** @brief Destructor */
  ~Generator();

  /** @brief Initializes the generator **/
  void Initialize();

  /** @brief Returns the central totam momentum of the beam generator **/
  const double& GetMomentum() const				{ return _p; }

  /** @brief Sets the central totam momentum of the beam generator **/
  void SetMomentum(const double& p)				{ _p = p; }

  /** @brief Returns the covariance matrix of the beam generator **/
  const Matrix<double>& GetCovarianceMatrix() const		{ return _S; }

  /** @brief Sets the covariance matrix of the beam generator **/
  void SetCovarianceMatrix(const Matrix<double>& S)		{ _S = S; }

  /** @brief Build the covariance matrix from the Penn parametrization of it
   *
   *  @param	p	Total beam momentum [MeV/c]
   *  @param	eps	Transverse normalised emittance [mm]
   *  @param	beta	Transverse beta function [mm]
   *  @param	alpha	Transverse alpha function (correlation)
   **/
  void SetMatrixParametrisation(const double& p,
				const double& eps,
		      	  	const double& beta,
		      	  	const double alpha=0.);

  /** @brief Returns the underlying data of a bunch sampled from a Gaussian of covariance S
   *
   *  @param	n	Number of particles to generate
   *
   *  @return		Dictionary that maps phase space coordinates to vectors
   */
  Beam::BunchMap GaussianBunchMap(const size_t n) const;

  /** @brief Returns a beam bunch sampled from a Gaussian of covariance S
   *
   *  @param	n	Number of particles to generate
   *
   *  @return		Beam bunch
   */
  Beam::Bunch GaussianBunch(const size_t n) const;

  /** @brief Returns the underlying data of a bunch sampled from a Spiral of covariance S
   *
   *  @param	n	Number of particles to generate
   *  @param	t	Turning strength
   *
   *  @return		Dictionary that maps phase space coordinates to vectors
   */
  Beam::BunchMap SpiralBunchMap(const size_t n,
				const double& t) const;

  /** @brief Returns a beam bunch sampled from a Spiral of covariance S
   *
   *  @param	n	Number of particles to generate
   *  @param	t	Turning strength
   *
   *  @return		Dictionary that maps phase space coordinates to vectors
   */
  Beam::Bunch SpiralBunch(const size_t n,
			  const double& t) const;

 private:
  size_t		_dim;		///< Dimension of the transverse phase space
  double		_mass;		///< Particle mass
  double		_p;		///< Central momentum of the beam
  Matrix<double>	_S;		///< Covariance matrix of the generator
};

#endif
