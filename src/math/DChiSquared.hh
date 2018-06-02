#ifndef DCHISQUARED_HH
#define DCHISQUARED_HH

// C++ includes
#include <iostream>
#include <vector>
#include <cmath>

// ROOT includes
#include "TMath.h"

// Other includes
#include "DFunction.hh"

/** @brief Defines the methods associated with the Chi-Squared distribution
 */
class DChiSquared : public DFunction {
 public:

  /** @brief Basic constructor, initilizes a Chi-Squared distribution of n DOF
   *
   *  @param 	ndf	Number of degrees of freedom
   */
  DChiSquared(const size_t ndf);

  /** @brief Copy constructor */
  DChiSquared(const DChiSquared& df);

  /** @brief Equality operator */
  DChiSquared& operator=(const DChiSquared& df);

  /** @brief Destructor */
  ~DChiSquared();

  /** @brief Cumulative distribution function at v
   *
   *  @param 	v	Value in which to evaluate the CDF
   */
  double CDF(const double& v) const;

  /** @brief If the function has an L1, L2 or Linf symmetry, returns its CDF at radius R
   *
   *  @param	R	Radius to evaluate at
   */
  double CDFRadial(const double& R) const;

  /** @brief Segment that encompass a certain p-value in 1 dimension
   *
   *  @param 	alpha	P-value
   */
  std::vector<TLine*> Contour1D(const double alpha) const;

  /** @brief Dummy definition for the 2D contour, chi^2 does not generalise to 2D
   *
   *  @param 	alpha	P-value
   */
  std::vector<TObject*> Contour2D(const double alpha) const;

  /** @brief Dummy definition for the 3D contour, chi^2 does not generalise to 3D
   *
   *  @param 	alpha	P-value
   */
  TF3* Contour3D(const double alpha) const;

  /** @brief Segment that encompass a certain p-value in radius
   *
   *  @param 	alpha	P-value
   */
  std::vector<TLine*> ContourRadial(const double alpha) const;

  /** @brief Volume of the smallest contour of a given p-value
   *
   *  @param 	alpha	P-value
   */
  double ContourVolume(const double alpha) const;

  /** @brief Returns the value of the distribution in v
   *
   *  @param	v	Vector in which to evaluate
   */
  double Evaluate(const std::vector<double>& v) const;

  /** @brief Returns the level of the alpha-contour of the distribution
   *
   *  @param 	alpha	P-value
   */
  double Level(const double alpha) const;

  /** @brief Returns the probability content within the defined interval */
  double Norm() const;

  /** @brief Returns the same as Evaluate as it is one dimensional */
  double Radial(const double R) const;

  /** @brief Returns the radii of the alpha-contour
   *
   *  @param 	alpha	P-value
   *  @param	r	Inner radius of the interval
   *  @param	R	Outer radius of the interval
   *  @param	l	Level of the contour
   */
  void Radii(const double alpha, double& r, double& R, double& l) const;

  /** @brief Returns a random value sampled from the Chi-squared distribution */
  double Random();

  /** @brief Returns a random vector sampled from the Chi-squared distribution */
  std::vector<double> RandomVector();

 private:

  /** @brief Initializes the distribution */
  bool Initialize();

  size_t	_ndf;		///< Number of degrees of freedom
  TF1		_frand;		///< Function used to generate random samples
};

#endif
