#ifndef DEXPONENTIAL_HH
#define DEXPONENTIAL_HH

// C++ includes
#include <iostream>
#include <vector>
#include <cmath>

// ROOT includes
#include "TPolyLine.h"

// Other includes
#include "Statistics.hh"
#include "DFunction.hh"

/** @brief Defines the methods associated with the exponential distribution in any dimension
 */
class DExponential : public DFunction {
 public:
  /** @brief Basic constructor, initilizes an exponential distribution of dimension n
   *
   *  @param 	n	Dimension of the space
   */
  DExponential(const size_t n=1);

  /** @brief 1D constructor, initilizes an exponential distribution of mean mu and scale lambda
   *
   *  @param 	mean	Mean of the distribution
   *  @param 	lambda	Scale of the distribution
   */
  DExponential(const double mean,
	       const double lambda);

  /** @brief nD constructor, initilizes an exponential distibution of given means and scales
   *
   *  @param 	mean	Mean of the distribution
   *  @param 	lambda	Scales of the distribution
   */
  DExponential(const std::vector<double>& mean,
	       const std::vector<double>& lambda);

  /** @brief Copy constructor */
  DExponential(const DExponential& df);

  /** @brief Equality operator */
  DExponential& operator=(const DExponential& df);

  /** @brief Destructor */
  ~DExponential();

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

  /** @brief Rhombus that encompasses a certain p-value in 2 dimensions
   *
   *  @param 	alpha	P-value
   */
  std::vector<TObject*> Contour2D(const double alpha) const;

  /** @brief n-Rhombus that encompasses a certain p-value in 3 dimensions
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

  /** @brief Return the value of the Exponential for an R following the L1 symmetry
   *
   *  @param	R	Radius to evaluate at (L1 symmetry)
   */
  double Radial(const double R) const;

  /** @brief Returns the radius of the alpha-contour, following an L1 symmetry
   *
   *  @param 	alpha	P-value
   */
  double Radius(const double alpha) const;

  /** @brief Returns a random value sampled from the exponential distribution */
  double Random();

  /** @brief Returns a random vector sampled from the exponential distribution */
  std::vector<double> RandomVector();

 private:

  /** @brief Initializes the distribution */
  bool Initialize();

  /** @brief Retuns the 1D normalisation for an arbitrary interval
   *
   *  @param	i	Axis to consider
   */
  double Norm1D(const size_t i) const;

  /** @brief Returns the parameters necessary to produce a ROOT 3D contour
   *
   *  @param	C	Level of the contour
   */
  std::vector<double> OffsetParameters(const double C) const;

  std::vector<double>	_mean;		///< Vector mean
  std::vector<double>	_lambda;	///< Vector of scale factors
  double		_prod;		///< Product of the scale factors
};

#endif
