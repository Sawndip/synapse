#ifndef DCAUCHY_HH
#define DCAUCHY_HH

// C++ includes
#include <iostream>
#include <vector>
#include <cmath>

// ROOT includes
#include "TBox.h"

// Other includes
#include "DFunction.hh"

/** @brief Defines the methods associated with the Cauchy distribution. 
 *
 *         Only defined in 1D. It cannot be extended to n dimensions as it is not normalisable.
 */
class DCauchy : public DFunction {
 public:

  /** @brief Basic constructor, initilizes a Cauchy distribution */
  DCauchy();

  /** @brief Standard constructor, initilizes a Cauchy distribution of median mu and width sigma
   *
   *  @param 	mu	Median of the distribution
   *  @param 	sigma	Width of the distribution
   */
  DCauchy(const double mu,
	  const double sigma);

  /** @brief Copy constructor */
  DCauchy(const DCauchy& df);

  /** @brief Equality operator */
  DCauchy& operator=(const DCauchy& df);

  /** @brief Destructor */
  ~DCauchy();

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

  /** @brief Dummy definition for the 2D contour, Cauchy is not normalisable in 2D
   *
   *  @param 	alpha	P-value
   */
  std::vector<TObject*> Contour2D(const double alpha) const;

  /** @brief Dummy definition for the 3D contour, Cauchy is not normalisable in 3D
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

  /** @brief Returns a random value sampled from the Cauchy distribution */
  double Random();

  /** @brief Returns a random vector sampled from the Cauchy distribution */
  std::vector<double> RandomVector();

 private:

  /** @brief Initializes the distribution */
  bool Initialize();

  double	_mean;		///< Median
  double	_sigma;		///< Width
};

#endif
