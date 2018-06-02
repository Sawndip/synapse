#ifndef DMULTIGAUS_HH
#define DMULTIGAUS_HH

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
#include "Statistics.hh"
#include "Geometry.hh"

/** @brief Defines the methods associated with a set of Gaussian distributions in any dimension
 */
class DMultiGaus : public DFunction {
 public:

  /** @brief Normal constructor, initilizes k normal distribution of dimension n
   *
   *  @param 	n	Dimension of the space
   *  @param 	k	Number of Gaussian peaks
   *  @param 	dist	Distance between successive peaks
   */
  DMultiGaus(const size_t n=1,
	     const size_t k=1,
	     const double dist=5);

  /** @brief 1D constructor, initilizes k Gaussians with provided means and sigmas
   *
   *  @param 	mu	Mean of the distribution
   *  @param 	sigma	RMS of the distribution
   */
  DMultiGaus(const std::vector<double>& mu,
	     const std::vector<double>& sigma);

  /** @brief nD constructor, initilizes k n-Gaussianswith provided means and sigmas
   *
   *  @param 	mu	Mean of the distribution
   *  @param 	sigma	RMS of the distribution in each dimension
   */
  DMultiGaus(const std::vector<std::vector<double>>& mu,
	     const std::vector<std::vector<double>>& sigma);

  /** @brief General constructor, intializes k n-Gaussians with provided means and cov. matrices
   *
   *  @param 	mu	Mean of the distribution
   *  @param 	cov	Covariance matrix of the distribution
   */
  DMultiGaus(const std::vector<std::vector<double>>& mu,
	     const std::vector<Matrix<double>>& cov);

  /** @brief Copy constructor */
  DMultiGaus(const DMultiGaus& df);

  /** @brief Equality operator */
  DMultiGaus& operator=(const DMultiGaus& df);

  /** @brief Destructor */
  ~DMultiGaus();

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

  /** @brief Adds a Gaussian peak to the current array of peaks
   *
   *  @param 	mu		Mean of the distribution
   *  @param 	sigma		RMS of the distribution
   */
  bool AddPeak(const std::vector<double>& mu, const std::vector<double>& sigma);

  /** @brief Adds a Gaussian peak to the current array of peaks
   *
   *  @param 	mu		Mean of the distribution
   *  @param 	cov		Covariance matrix of the distribution
   */
  bool AddPeak(const std::vector<double>& mu, const Matrix<double>& cov);

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

  /** @brief Dummy definition of radial contour, no radial symmetry
   *
   *  @param 	alpha	P-value
   */
  std::vector<TLine*> ContourRadial(const double alpha) const;

  /** @brief Ellipsoid that encompasses a certain p-value in 3 dimensions
   *
   *  @param 	alpha	P-value
   */
  TF3* Contour3D(const double alpha) const;

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

  /** @brief Radial function. The uniform distribution does not have a radial symmetry */
  double Radial(const double R) const;

  /** @brief Return the value of the Gaussian for an R following the L2 symmetry
   *
   *  @param	i	ID of the peak of which to compute the level
   *  @param	R	Radius to evaluate at (L2 symmetry)
   */
  double Radial(const size_t i,
		const double R) const;

  /** @brief Returns the radius of the alpha-contour, following an L2 symmetry
   *
   *  @param 	alpha	P-value
   */
  double Radius(const double alpha) const;

  /** @brief Returns the radii of each peak for the alpha-contour, following an L2 symmetry
   *
   *  @param 	alpha	P-value
   */
  std::vector<double> Radii(const double alpha) const;

  /** @brief Returns a random value sampled from the Gaussian distribution */
  double Random();

  /** @brief Returns a random vector sampled from the Gaussian distribution */
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

  size_t				_k;		///< Number of peaks
  std::vector<std::vector<double>>	_means;		///< Array of vector means
  std::vector<Matrix<double>> 		_covs;		///< Array of covariance matrices
  std::vector<Matrix<double>> 		_invcovs;	///< Array of inverse covariance matrices
  std::vector<Matrix<double>>		_Js;		///< Array of Jacobian of variable changes
  std::vector<double>			_dets;		///< Array of determinants
};

#endif
