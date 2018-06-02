#ifndef DSPIRAL_HH
#define DSPIRAL_HH

// C++ includes
#include <iostream>
#include <vector>
#include <cmath>

// ROOT includes
#include "TMath.h"
#include "TEllipse.h"

// Other includes
#include "DFunction.hh"
#include "Matrix.hh"

/** @brief Defines the methods associated with the Spiral distribution in any dimension.
 *
 *	   This distribution is based on an n dimensional Gaussian that has a radial
 * 	   dependant rotation applied to it. The rotation preserves the normalisation of the
 *         Gaussian. It adopts the shape of a spiral in any 2 dimensional projections.
 *	   It is a Fermat spiral as the rotation angle grows as \f$\theta= r^2/a^2\f$.
 *
 *	   The radius of this density function is defined as the product,
 *	   \f$\vec{x}^{T}(T\Sigma^{-1}T^{T})\vec{x}\f$, with \f$x\f$ the phase space vector,
 *  	   \f$\Sigma\f$ the covariance matrix and \f$T\f$ the twisting matrix,
 *	   function of the \f$L^{2}\f$ squared radius, \f$R^2=\vec{x}^T\vec{x}\f$.
 */
class DSpiral : public DFunction {
 public:

  /** @brief Normal constructor, initilizes a Spiral distribution of dimension n
   *
   *  @param 	n	Dimension of the space
   *  @param	turn 	Turning strength
   */
  DSpiral(const size_t n=1,
	  const double turn=1);

  /** @brief 1D constructor, initilizes a Spiral of mean mu and RMS sigma
   *
   *	     This is identical to a 1D gaussian, a rotation in 1D is always the identity.
   *
   *  @param 	mu	Mean of the distribution
   *  @param 	sigma	RMS of the distribution
   *  @param 	turn	Turning strength
   */
  DSpiral(const double mu,
	  const double sigma,
	  const double turn=1);

  /** @brief Uncorrelated constructor, initilizes an n-Spiral of mean mu and RMS sigma
   *
   *  @param 	mu	Mean of the distribution
   *  @param 	sigma	RMS of the distribution in each dimension
   *  @param 	turn	Turning strength
   */
  DSpiral(const std::vector<double>& mu,
	  const std::vector<double>& sigma,
	  const double turn=1);

  /** @brief General constructor, intializes mean and the covariance matrix of the distribution
   *
   *  @param 	mu	Mean of the distribution
   *  @param 	cov	Covariance matrix of the distribution
   *  @param 	turn	Turning strength
   */
  DSpiral(const std::vector<double>& mu,
	  const Matrix<double>& cov,
	  const double turn=1);

  /** @brief Copy constructor */
  DSpiral(const DSpiral& df);

  /** @brief Equality operator */
  DSpiral& operator=(const DSpiral& df);

  /** @brief Destructor */
  ~DSpiral();

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

  /** @brief Ellipse that encompasses a certain p-value in 2 dimensions
   *
   *  @param 	alpha	P-value
   */
  std::vector<TObject*> Contour2D(const double alpha) const;

  /** @brief Ellipsoid that encompasses a certain p-value in 3 dimensions
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

  /** @brief Return the value of the spiral for an R following the L2 symmetry
   *
   *  @param	R	Radius to evaluate at (L2 symmetry)
   */
  double Radial(const double R) const;

  /** @brief Returns the radius of the alpha-contour, following an L2 symmetry
   *
   *  @param 	alpha	P-value
   */
  double Radius(const double alpha) const;

  /** @brief Returns a random value sampled from the spiral distribution */
  double Random();

  /** @brief Returns a random vector sampled from the spiral distribution */
  std::vector<double> RandomVector();

 private:

  /** @brief Initializes the distribution */
  bool Initialize();

  /** @brief Retuns the 1D normalisation for an arbitrary interval
   *
   *  @param	i	Axis to consider
   */
  double Norm1D(const size_t i) const;

  /** @brief Returns the parameters to produce a ROOT type spiral function */
  std::vector<double> Parameters() const;

  /** @brief Returns the parameters necessary to produce a ROOT 3D contour 
   *
   *  @param	C	Level of the contour
   */
  std::vector<double> OffsetParameters(const double C) const;

  /** @brief Returns the twisting matrix at a given squared radius R2
   *
   * @param	R2	Squared radius
   */
  Matrix<double> TwistingMatrix(const double R2) const;

  std::vector<double>	_mean;		///< Vector mean
  Matrix<double> 	_cov;		///< Covariance matrix
  Matrix<double> 	_invcov;	///< Inverse covariance matrix
  Matrix<double>	_J;		///< Jacobian of variable change
  double		_det;		///< Determinant of the covariance matrix
  double		_turn;		///< Controls the turning strength
};

#endif
