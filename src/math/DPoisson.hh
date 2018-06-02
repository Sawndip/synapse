#ifndef DPoisson_HH
#define DPoisson_HH

// C++ includes
#include <iostream>
#include <vector>
#include <cmath>

// ROOT includes
#include "TBox.h"
#include "TMath.h"

// Other includes
#include "DFunction.hh"

/** @brief Defines the methods associated with the Poisson distribution.
 */
class DPoisson : public DFunction {
 public:

  /** @brief Basic constructor, initilizes a Poisson distribution */
  DPoisson();

  /** @brief 1D constructor, initilizes a Poisson distribution of mean mu
   *
   *  @param 	mu	Mean of the distribution
   */
  DPoisson(const double mu);

  /** @brief Copy constructor */
  DPoisson(const DPoisson& df);

  /** @brief Equality operator */
  DPoisson& operator=(const DPoisson& df);

  /** @brief Destructor */
  ~DPoisson();

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

  /** @brief Dummy definition for the 2D contour
   *
   *  @param 	alpha	P-value
   */
  std::vector<TObject*> Contour2D(const double alpha) const;

  /** @brief Dummy definition for the 3D contour
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

  /** @brief Returns the centred and normalised version of Evaluate
   *
   *  @param	R	Radius to evaluate at
   */
  double Radial(const double R) const;

  /** @brief Returns the radii of the alpha-contour
   *
   *  @param 	alpha	P-value
   *  @param	r	Inner radius of the interval
   *  @param	R	Outer radius of the interval
   *  @param	l	Level of the contour
   */
  void Radii(const double alpha, double& r, double& R, double& l) const;

  /** @brief Returns a random value sampled from the Poisson distribution */
  double Random();

  /** @brief Returns a random vector sampled from the Poisson distribution */
  std::vector<double> RandomVector();

 private:

  /** @brief Initializes the distribution */
  bool Initialize();

  double	_mean;		///< Mean
};

#endif
