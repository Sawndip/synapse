#ifndef DMAXWELL_HH
#define DMAXWELL_HH

// C++ includes
#include <iostream>
#include <vector>
#include <cmath>

// ROOT includes
#include "TMath.h"
#include "TEllipse.h"

#include "TCanvas.h"
#include "TH1F.h"

// Other includes
#include "DFunction.hh"
#include "Matrix.hh"

/** @brief Defines the methods associated with the the Maxwell-Bolztmann distribution
	   in any dimension
 */
class DMaxwell : public DFunction {
 public:

  /** @brief Normal constructor, initilizes a normal distribution of dimension n
   *
   *  @param 	n	Dimension of the space
   */
  DMaxwell(const size_t n=1);

  /** @brief 1D constructor, initilizes a Gaussian mean mu and RMS sigma
   *
   *  @param 	mu	Mean of the distribution
   *  @param 	sigma	RMS of the distribution
   */
  DMaxwell(const double mu,
	   const double sigma);

  /** @brief Uncorrelated constructor, initilizes an n-Gaussian of mean mu and RMS sigma
   *
   *  @param 	mu	Mean of the distribution
   *  @param 	sigma	RMS of the distribution in each dimension
   */
  DMaxwell(const std::vector<double>& mu,
	   const std::vector<double>& sigma);

  /** @brief General constructor, intializes mean and the covariance matrix of the distribution
   *
   *  @param 	mu	Mean of the distribution
   *  @param 	cov	Covariance matrix of the distribution
   */
  DMaxwell(const std::vector<double>& mu,
	   const Matrix<double>& cov);

  /** @brief Copy constructor */
  DMaxwell(const DMaxwell& df);

  /** @brief Equality operator */
  DMaxwell& operator=(const DMaxwell& df);

  /** @brief Destructor */
  ~DMaxwell();

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

  /** @brief Segments that encompass a certain p-value in 1 dimension
   *
   *  @param 	alpha	P-value
   */
  std::vector<TLine*> Contour1D(const double alpha) const;

  /** @brief Concentric ellipses that encompass a certain p-value in 2 dimensions
   *
   *  @param 	alpha	P-value
   */
  std::vector<TObject*> Contour2D(const double alpha) const;

  /** @brief Concentric ellipsoids that encompasses a certain p-value in 3 dimensions
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

  /** @brief Return the value of the Maxwell distribution for an R following the L2 symmetry
   *
   *  @param	R	Radius to evaluate at (L2 symmetry)
   */
  double Radial(const double R) const;

  /** @brief Returns the radii of the alpha-contour, following an L2 symmetry
   *
   *  @param 	alpha	P-value
   *  @param	r	Inner radius of the shell
   *  @param	R	Outer radius of the shell
   *  @param	l	Level of the contour
   */
  void Radii(const double alpha, double& r, double& R, double& l) const;

  /** @brief Returns a random value sampled from the Maxwell distribution */
  double Random();

  /** @brief Returns a random vector sampled from the Maxwell distribution */
  std::vector<double> RandomVector();

 private:

  /** @brief Initializes the distribution */
  bool Initialize();

  /** @brief Retuns the 1D normalisation for an arbitrary interval
   *
   *  @param	i	Axis to consider
   */
  double Norm1D(const size_t i) const;

  /** @brief Returns the parameters to produce a ROOT type Gaussian function */
  std::vector<double> Parameters() const;

  /** @brief Returns the parameters necessary to produce a ROOT 3D contour 
   *
   *  @param	C	Level of the contour
   */
  std::vector<double> OffsetParameters(const double C) const;

  std::vector<double>	_mean;		///< Vector mean
  Matrix<double> 	_cov;		///< Covariance matrix
  Matrix<double> 	_invcov;	///< Inverse covariance matrix
  Matrix<double>	_J;		///< Jacobian of variable change
  double		_det;		///< Determinant of the covariance matrix
  TF1			_frand;		///< Function used to generate random 1D samples
};

#endif
