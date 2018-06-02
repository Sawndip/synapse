#ifndef BUNCH_HH
#define BUNCH_HH

// C++ includes
#include <stdlib.h>
#include <map>
#include <vector>

// QHULL includes
#include "Qhull.h"

// ROOT includes
#include "TEllipse.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"

// Other includes
#include "Statistics.hh"
#include "DGaus.hh"
#include "DChiSquared.hh"
#include "DensityEstimator.hh"
#include "ProbabilityContour.hh"
#include "BunchData.hh"
#include "Variable.hh"
#include "Drawer.hh"

/** @brief Environment of beam related objects and methods */
namespace Beam {

/** @brief Holds a bunch of particles at a specific position along the beam line.
 *
 *	   Quantities associated with the bunch such as Twiss parameters, amplitudes,
 * 	   emittance, volume, etc. may be computed and stored here. This avoids unnecessarily
 * 	   repeating operations such as computing the particle amplitudes.
 *
 *	   It is designed to hold a single species of particles. 
 */
class Bunch {
 public:
  /** @brief Default constructor, initilizes everything to 0 */
  Bunch();

  /** @brief Map constructor, sets the particle sample from maps (no measurement errors)
   *
   *  @param	samples		Dictionary that maps each variable onto a vector of values
   *  @param	name		Name of the bunch
   *  @param	pos		z position of the bunch [mm]
   *  @param	mass		Mass of the particles that compose the bunch \f$[\textrm{MeV}/c^2]\f$
   */
  Bunch(const BunchMap& samples,
	const double pos = 0,
        const std::string& name = "",
        const double mass = 105.66);

  /** @brief Map constructor, sets the particle sample from maps
   *
   *  @param	samples		Dictionary that maps each variable onto a vector of values
   *  @param	errors		Dictionary that maps each variable onto a vector of uncertainties
   *  @param	pos		z position of the bunch [mm]
   *  @param	name		Name of the bunch
   *  @param	mass		Mass of the particles that compose the bunch \f$[\textrm{MeV}/c^2]\f$
   */
  Bunch(const BunchMap& samples,
        const BunchMap& errors,
	const double pos = 0,
        const std::string& name = "",
        const double mass = 105.66);

  /** @brief Copy constructor */
  Bunch(const Bunch& bunch);

  /** @brief Equality operator */
  Bunch& operator=(const Bunch& bunch);

  /** @brief Destructor */
  ~Bunch();

  /** @brief Returns the name of the bunch */
  const std::string& Name() const			{ return _name; }

  /** @brief Sets the name of the bunch */
  void SetName(const std::string& name)			{ _name = name; }

  /** @brief Returns the dimension of space in which the bunch lives */
  const size_t& Dimension() const			{ return _dim; }

  /** @brief Returns the z position at which the bunch was sampled [mm] */
  const double& Position() const			{ return _pos; }

  /** @brief Returns the mass of the particles that make up the beam [MeV/c^2] */
  const double& Mass() const				{ return _mass; }

  /** @brief Sets the position of the bunch */
  void SetPosition(const double& pos)			{ _pos = pos; }

  /** @brief Returns the size of the whole sample */
  const size_t Size() const				{ return _data.Size(); }

  /** @brief Returns the size of the core sample */
  const size_t CoreSize() const				{ return _core.Size(); }

  /** @brief Returns the array of samples of variable var */
  const std::vector<double>& Samples(const std::string var) const;

  /** @brief Returns the array of errors on the samples of variable var */
  const std::vector<double>& Errors(const std::string var) const;

  /** @brief Returns the radius of the i^th sample */
  double Radius(const size_t i) const;

  /** @brief Returns the array of radii of the samples */
  std::vector<double> Radii() const;

  /** @brief Returns the mean of the requested quantity */
  Variable Mean(const std::string& var) const;

  /** @brief Returns the covariance matrix of the bunch */
  const Matrix<double>& CovarianceMatrix() const	{ return _data.S(); }

  /** @brief Returns the Twiss parameter alpha of the bunch */
  const Variable& Alpha() const				{ return _alpha; }

  /** @brief Returns the Twiss parameter alpha of the bunch in one 2D projection */
  Variable Alpha(const std::string axis) const;

  /** @brief Returns the Twiss parameter beta of the bunch */
  const Variable& Beta() const				{ return _beta; }

  /** @brief Returns the Twiss parameter beta of the bunch in one 2D projection */
  Variable Beta(const std::string axis) const;

  /** @brief Returns the Twiss parameter gamma of the bunch */
  const Variable Gamma() const				{ return _gamma; }

  /** @brief Returns the Twiss parameter gamma of the bunch in one 2D projection */
  Variable Gamma(const std::string axis) const;

  /** @brief Returns the mechanical angular momentum of the bunch */
  const Variable& MechanicalL() const			{ return _mecl; }

  /** @brief Returns the RMS normalised emittance of the bunch */
  const Variable& Emittance() const			{ return _eps; }

  /** @brief Returns the RMS geometric emittance of the bunch */
  const Variable& NormEmittance() const			{ return _neps; }

  /** @brief Returns the subemittance of the core */
  const Variable& SubEmittance() const			{ return _seps; }

  /** @brief Returns the volume of the core */
  const Variable& FracEmittance() const			{ return _feps; }

  /** @brief Returns the maximum amplitude in the core */
  const Variable& AmplitudeQuantile() const		{ return _qamp; }

  /** @brief Returns the estimated normalised emittance from the measured fractional quantity */
  Variable NormEmittanceEstimate(const SumStat& stat) const;

  /** @brief Returns the requested summary statistic */
  Variable SummaryStatistic(const SumStat& stat) const;

  /** @brief Returns the array of amplitudes of the samples */
  const std::vector<double>& Amplitudes() const		{ return _amps; }

  /** @brief Returns the array of Voronoi volumes of the samples */
  const std::vector<double>& VoronoiVolumes() const	{ return _vols; }

  /** @brief Defines the fraction the core represents (in case of transmission losses) */
  void SetFraction(const double frac)			{ _fraction = frac; }

  /** @brief Defines the fraction of the beam contained in the core
   *
   *  @param	frac	Fraction of the initial bunch to be kept
   */
  void SetCoreFraction(const double frac);

  /** @brief Assigns a certain amount of core particles
   *
   *  @param	size	Number of particles to be kept
   */
  void SetCoreSize(const size_t size);

  /** @brief Computes the volume of the fraction of the beam contained in the core
   * 
   *  @param	algo	Volume reconstruction algorithm
   *  @param	norm	Normalisation constant if the bunch has scraped
   */
  void SetCoreVolume(const std::string algo="hull",
		     const double norm=1.);

  /** @brief Compute the corrected transverse amplitudes of every particle in the sample by 
	     only considering particles of lower amplitudes.
   */
  void SetCorrectedAmplitudes();

  /** @brief Compute the MCD transverse amplitude of every particle in the sample 
	     by selecting a subsample optimized based on its covariance determinant
   */
  void SetMCDAmplitudes();

  /** @brief Compute the generalised transverse amplitude of every particle in the sample 
	     by basing it on density estimation and contour volume rather than covariance.
   *
   *  @param 	norm	Normalisation constant (if the distribution is truncated)
   */
  void SetGeneralisedAmplitudes(const double norm=1.);

  /** @brief Compute the Voronoi volumes around each particle
   *
   *	     This is not done by default as it takes a long time (need to tesselate the 4D space).
   *
   *	     Does not do anything if it has already been called.
   */
  void SetVoronoiVolumes();

  /** @brief Returns the n-RMS ellipse of a sample of two variables
   *
   *  	     The default probability corresponds to the \f$1\sigma\f$ contour.
   *
   *  @param 	vara	First variable in which to compute the RMS ellipse
   *  @param 	varb	Second variable in which to compute the RMS ellipse
   *  @param	p	P-value of the distribution (probability content)
   *
   *  @return 		n-RMS TEllipse pointer
   */
  TEllipse* Ellipse(const std::string& vara,
		    const std::string& varb,
		    const double p=.39347) const;

  /** @brief Returns the n-RMS robust ellipse of a sample of two variables
   *
   *  @param 	vara	First variable in which to compute the RMS ellipse
   *  @param 	varb	Second variable in which to compute the RMS ellipse
   *  @param	p	P-value of the distribution (probability content)
   *
   *  @return 		Robust n-RMS TEllipse pointer
   */
  TEllipse* RobustEllipse(const std::string& vara,
			  const std::string& varb,
			  const double p=.39347) const;

  /** @brief Returns the 1D distrubution of the requested variable
   *
   *  @param 	var	Name of the variable of which to draw the distribution
   *  @param 	min	Minimum of the range of var
   *  @param 	max	Maximum of the range of var
   *
   *  @return 		TH1F pointer
   */
  TH1F* Histogram(const std::string& var,
		  double min=0,
		  double max=0) const;

  /** @brief Returns the 2D distrubution of the requested variable
   *
   *  @param 	vara	First variable of which to draw the distribution
   *  @param 	varb	Second variable of which to draw the distribution
   *  @param 	mina	Minimum of the range of vara
   *  @param 	maxa	Maximum of the range of vara
   *  @param 	minb	Minimum of the range of varb
   *  @param 	maxb	Maximum of the range of varb
   *
   *  @return 		ROOT TH2F pointer
   */
  TH2F* Histogram(const std::string& vara, 
		  const std::string& varb,
		  double mina=0, 
		  double maxa=0,
	 	  double minb=0,
		  double maxb=0) const;

  /** @brief Returns the 2D scatter plot of particle amplitudes
   *
   *  @param 	vara	First variable of which to draw the amplitude scatter plot
   *  @param 	varb	Second variable of which to draw the amplitude scatter plot
   *
   *  @return 		ScatterGraph pointer
   */
  ScatterGraph* AmplitudeScatter(const std::string& vara,
			     	 const std::string& varb) const;

 private:

  /** @brief Normal initializer, sets the particle sample and bunch variables
   *
   *  @param	samples		Dictionary that maps each variable onto a vector of values
   *  @param	errors		Dictionary that maps each variable onto a vector of errors
   */
  void Initialize(const BunchMap& samples,
		  const BunchMap& errors);

  /** @brief Checks if the variable is provided */
  bool Contains(const std::string& var) const		{ return _data.Contains(var); }

  /** @brief Checks if the variable is provided, throws if not */
  void AssertContains(const std::string& var) const	{ _data.AssertContains(var); }

  /** @brief Returns the value and unertainty arrays corresponding to the requested variable
   *
   *  @param	var	Variable of which to return the arrays
   *  @param	samples	Array of values
   *  @param	errors	Array of uncertainties
   */
  void Arrays(const std::string& var,
	      std::vector<double>& samples,
	      std::vector<double>& errors) const;

  /** @brief Returns the total momentum of a given sample */
  double TotalMomentum(const size_t i) const;

  /** @brief Returns the total momentum of all the samples */
  std::vector<double> TotalMomenta() const;

  /** @brief Returns the error on the total momentum of a given sample */
  double TotalMomentumError(const size_t i) const;

  /** @brief Returns the error on the total momentum of all the samples */
  std::vector<double> TotalMomentumErrors() const;

  /** @brief Computes the Twiss parameter alpha of the bunch */
  void SetTwiss();

  /** @brief Computes the measurement error on the requested Twiss parameter
   *
   *  @param 	name		Name of the Twiss parameter
   *  @param 	means		Map of the means for each variables
   *  @param 	detr		Errors on the determinant of the optical covariance matrix
   *  @param 	C		Cofactor matrix of the optical covariance matrix
   */
  double TwissMError(const std::string& name,
		     std::map<std::string, double>& means,
		     const double& detr,
		     const Matrix<double>& C) const;

  /** @brief Computes the statistical error on the provided Twiss parameter
   *
   *  @param 	twiss		Value of the twiss parameter
   */
  double TwissSError(const double& twiss) const;

  /** @brief Computes the RMS geometric emittance of the bunch */
  void SetEmittance();

  /** @brief Computes the RMS normalised emittance of the bunch */
  void SetNormEmittance();

  /** @brief Computes the subemittance of the core */
  void SetSubEmittance();

  /** @brief Sets the amplitude quantile and its uncertainty */
  void SetAmplitudeQuantile(const double& qamp);

  /** @brief Computes the RMS emittance provided with the covariance matrix */
  double NormEmittanceValue(const Matrix<double>& covmat) const;

  /** @brief Compute the transverse amplitudes of every particle in the sample */
  void SetAmplitudes();

  /** @brief Compute the transverse amplitude provided with S^-1 and the PS vector
   *
   *  @param 	invcovmat	Inverse of the covariance matrix
   *  @param 	means		Map of the means of each variables
   *  @param 	vec		Phase space vector of the particle
   *  @param 	neps		Normalised emittance of the bunch
   */
  double AmplitudeValue(const Matrix<double>& invcovmat,
			const std::vector<double>& means,
			const std::vector<double>& vec,
			const double& neps) const;

  std::string			_name;		///< Name of the bunch
  size_t			_dim;		///< Dimension of the transverse phase space
  double			_pos;		///< Position at which the bunch was sampled
  double			_mass;		///< Mass of the particles
  double			_fraction;	///< Fraction contained in the beam core
  BunchData			_data;		///< Underlying phase space vectors
  BunchData			_core;		///< Core fraction of beam
  Variable			_alpha;		///< Twiss parameter alpha
  Variable			_beta;		///< Twiss parameter beta
  Variable			_gamma;		///< Twiss parameter gamma 
  Variable			_mecl;		///< Normalised mechanical angular momentum
  Variable			_eps;		///< Geometric RMS emittance
  Variable			_neps;		///< Normalised RMS emittance
  Variable			_seps;		///< Subemittance
  Variable			_feps;		///< Fractional emittance
  Variable			_qamp;		///< Lagest amplitude in the core
  std::vector<double>  		_amps;		///< Amplitudes of the particles
  std::vector<double>		_vols;		///< Voronoi cell volumes of the particles
  DensityEstimator		_de;		///< Density estimation of the bunch
};
} // namespace Beam

#endif
