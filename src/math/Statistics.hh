#ifndef STATISTICS_HH
#define STATISTICS_HH

// C++ includes
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <ctime>
#include <algorithm>

// ROOT includes
#include "TMath.h"
#include "TF1.h"
#include "TRobustEstimator.h"

// Other includes
#include "Assert.hh"
#include "Matrix.hh"

/** @file Statistics.hh 
 *
 *  @brief Set of functions that compose the Statistics package.
 */

/** @brief Mathematics environment (includes Statistics and Geometry packages)
 */
namespace Math {

/** @brief Computes the mean of a variable sample
 *
 *  @param	vx	Array of x
 *
 *  @return		Mean
 */ 
double Mean(const std::vector<double>& vx);

/** @brief Computes the error on the mean of a variable sample
 *
 *  @param	vx	Array of x
 *  @param	vex	Array of x errors
 *
 *  @return		Uncertainty on the mean
 */ 
double MeanError(const std::vector<double>& vx,
		 const std::vector<double>& vex=std::vector<double>());

/** @brief Computes the measurement error on the mean of a variable sample
 *
 *  @param	vex	Array of x errors
 *
 *  @return		Measurement uncertainty on the mean
 */ 
double MeanMError(const std::vector<double>& vex);

/** @brief Computes the statistical error on the mean of a variable sample
 *
 *  @param	vx	Array of x
 *
 *  @return		Statistical uncertainty on the mean
 */ 
double MeanSError(const std::vector<double>& vx);

/** @brief Computes the MCD robust mean of a variable sample
 *
 *  @param	vx	Array of x
 *
 *  @return		Mean
 */ 
double RobustMean(const std::vector<double>& vx);

/** @brief Removes an element from a precomputed mean
 *
 *  @param	n	Current number of elements in the mean
 *  @param	mean	Value of the mean of n elements
 *  @param	x	Value to remove from the sample
 */ 
void DecrementMean(const size_t& n,
		   double& mean,
		   const double& x);

/** @brief Adds an element to a precomputed mean
 *
 *  @param	n	Current number of elements in the mean
 *  @param	mean	Value of the mean of n elements
 *  @param	x	Value to add to the sample
 */ 
void IncrementMean(const size_t& n,
		   double& mean,
		   const double& x);

/** @brief Computes the trimmed mean of a variable sample
 *
 *  @param	vx	Array of x
 *  @param	frac	Fraction of lowest and highest values to discard
 *
 *  @return		Trimmed mean
 */ 
double TMean(const std::vector<double>& vx,
	     const double& frac);

/** @brief Computes the statistical error on the trimmed mean of a variable sample
 *
 *  @param	vx	Array of x
 *  @param	frac	Fraction of lowest and highest values to discard
 *
 *  @return		Statistical uncertainty on the trimmed mean
 */ 
double TMeanSError(const std::vector<double>& vx,
		   const double& frac);

/** @brief Computes the winsorized mean of a variable sample
 *
 *  @param	vx	Array of x
 *  @param	frac	Fraction of lowest and highest values to winsorize
 *
 *  @return		Winsorized mean
 */ 
double WMean(const std::vector<double>& vx,
	     const double& frac);

/** @brief Computes the variance of a variable sample
 *
 *  @param	vx	Array of x
 *
 *  @return		Variance
 */ 
double Variance(const std::vector<double>& vx);

/** @brief Computes the error on the variance of a variable sample
 *
 *  @param	vx	Array of x
 *  @param	vex	Array of x errors
 *
 *  @return		Uncertainty on the variance
 */ 
double VarianceError(const std::vector<double>& vx,
		     const std::vector<double>& vex=std::vector<double>());

/** @brief Computes the measurement error on the variance of a variable sample
 *
 *  @param	vx	Array of x
 *  @param	vex	Array of x errors
 *
 *  @return		Measurement uncertainty on the variance
 */ 
double VarianceMError(const std::vector<double>& vx,
		      const std::vector<double>& vex=std::vector<double>());

/** @brief Computes the statistical error on the variance of a variable sample
 *
 *  @param	vx	Array of x
 *
 *  @return		Statistical uncertainty on the variance
 */ 
double VarianceSError(const std::vector<double>& vx);

/** @brief Computes the standard deviation of a variable sample
 *
 *  @param	vx	Array of x
 *
 *  @return		Standard deviation
 */ 
double RMS(const std::vector<double>& vx);

/** @brief Computes the error on the standard deviation of a variable sample
 *
 *  @param	vx	Array of x
 *  @param	vex	Array of x errors
 *
 *  @return		Uncertainty on the standard deviation
 */ 
double RMSError(const std::vector<double>& vx,
		const std::vector<double>& vex=std::vector<double>());

/** @brief Computes the measurement error on the standard deviation of a variable sample
 *
 *  @param	vx	Array of x
 *  @param	vex	Array of x errors
 *
 *  @return		Measurement uncertainty on the standard deviation
 */ 
double RMSMError(const std::vector<double>& vx,
		 const std::vector<double>& vex);

/** @brief Computes the statistical error on the standard deviation of a variable sample
 *
 *  @param	vx	Array of x
 *
 *  @return		Statistical uncertainty on the standard deviation
 */ 
double RMSSError(const std::vector<double>& vx);

/** @brief Computes the winsorized standard deviation of a variable sample
 *
 *  @param	vx	Array of x
 *  @param	frac	Fraction of lowest and highest values to winsorize
 *
 *  @return		Standard deviation
 */ 
double WRMS(const std::vector<double>& vx,
	    const double& frac);

/** @brief Computes the covariance of two variable samples
 *
 *  @param	vx	Array of x
 *  @param	vy	Array of y
 *
 *  @return		Covariance
 */ 
double Covariance(const std::vector<double>& vx,
		  const std::vector<double>& vy);

/** @brief Computes the error on the covariance of two variable samples
 *
 *  @param	vx	Array of x
 *  @param	vy	Array of y
 *  @param	vex	Array of x errors
 *  @param	vey	Array of y errors
 *
 *  @return		Uncertainty on the covariance
 */ 
double CovarianceError(const std::vector<double>& vx,
		       const std::vector<double>& vy,
		       const std::vector<double>& vex=std::vector<double>(),
		       const std::vector<double>& vey=std::vector<double>());

/** @brief Computes the measurement error on the covariance of two variable samples
 *
 *  @param	vx	Array of x
 *  @param	vy	Array of y
 *  @param	vex	Array of x errors
 *  @param	vey	Array of y errors
 *
 *  @return		Measurement uncertainty on the covariance
 */ 
double CovarianceMError(const std::vector<double>& vx,
			const std::vector<double>& vy,
			const std::vector<double>& vex,
			const std::vector<double>& vey);

/** @brief Computes the statistical error on the covariance of two variable samples
 *
 *  @param	vx	Array of x
 *  @param	vy	Array of y
 *
 *  @return		Statistical uncertainty on the covariance
 */ 
double CovarianceSError(const std::vector<double>& vx,
			const std::vector<double>& vy);

/** @brief Adds an pair of elements from a precomputed covariance
 *
 *  @param	n	Current number of pair of elements in the covariance
 *  @param	cov	Value of the covariance of n elements
 *  @param	xmean	Current mean of the first variable
 *  @param	ymean	Current mean of the second variable
 *  @param	x	Value of the first variable to remove from the sample
 *  @param	y	Value of the second variable to remove from the sample
 */ 
void IncrementCovariance(const size_t& n,
		     	 double& cov,
		     	 const double& xmean,
		     	 const double& ymean,
		     	 const double& x,
		     	 const double& y);

/** @brief Removes an pair of elements from a precomputed covariance
 *
 *  @param	n	Current number of pair of elements in the covariance
 *  @param	cov	Value of the covariance of n elements
 *  @param	xmean	Current mean of the first variable
 *  @param	ymean	Current mean of the second variable
 *  @param	x	Value of the first variable to remove from the sample
 *  @param	y	Value of the second variable to remove from the sample
 */ 
void DecrementCovariance(const size_t& n,
		     	 double& cov,
		     	 const double& xmean,
		     	 const double& ymean,
		     	 const double& x,
		     	 const double& y);

/** @brief Returns the covariance matrix of a multivariate sample
 *
 *  @param	mat	Matrix that contains the relevant sample
 *
 *  @return		Covariance matrix
 */ 
Matrix<double> CovarianceMatrix(const Matrix<double>& mat);

/** @brief Returns the MCD robust covariance matrix of a multivariate sample
 *
 *  @param	mat	Matrix that contains the relevant sample
 *  @param 	hh	Number of elements in the subsamble of which to get the MCD
 *
 *  @return		Robust covariance matrix
 */ 
Matrix<double> RobustCovarianceMatrix(const Matrix<double>& mat,
				      const double hh=0);

/** @brief Removes an element array from a precomputed covariance matrix
 *
 *  @param	n	Current number of element arrays in the covariance matrix
 *  @param	covmat	Value of the covariance matrix
 *  @param	means	Vector of means
 *  @param	v	Element array to remove from the sample
 */ 
void DecrementCovarianceMatrix(const size_t& n,
		     	       Matrix<double>& covmat,
		     	       const std::vector<double>& means,
		     	       const std::vector<double>& v);

/** @brief Adds an element array to a precomputed covariance matrix
 *
 *  @param	n	Current number of element arrays in the covariance matrix
 *  @param	covmat	Value of the covariance matrix
 *  @param	means	Vector of means
 *  @param	v	Element array to add to the sample
 */ 
void IncrementCovarianceMatrix(const size_t& n,
		     	       Matrix<double>& covmat,
		     	       const std::vector<double>& means,
		     	       const std::vector<double>& v);

/** @brief Returns the exact error on the determinant of the covariance matrix
 *
 *  @param	mat		Matrix that contains the samples
 *  @param	mat_error	Matrix that contains the errors on those samples
 *
 *  @return			Uncertainty on the determinant of the covariance matrix
 */ 
double DeterminantError(const Matrix<double>& mat,
			const Matrix<double>& mat_error=Matrix<double>());

/** @brief Returns the measurement error on the determinant of the covariance matrix
 *
 *  @param	mat		Matrix that contains the samples
 *  @param	mat_error	Matrix that contains the errors on those samples
 *
 *  @return			Measurement uncertainty on the determinant of the covariance matrix
 */ 
double DeterminantMError(const Matrix<double>& mat,
			 const Matrix<double>& mat_error);

/** @brief Returns the statistical error on the determinant of the covariance matrix
 *
 *  @param	mat		Matrix that contains the samples
 *
 *  @return			Statistical uncertainty on the determinant of the covariance matrix
 */ 
double DeterminantSError(const Matrix<double>& mat);

/** @brief Computes the median of a variable sample
 *
 *  @param	vx	Array of x
 *
 *  @return		Median
 */ 
double Median(const std::vector<double>& vx);

/** @brief Computes the statistical error on the median of a variable sample
 *
 *  @param	n	Total number of samples points
 *  @param	density	Density at the point in question (essential)
 *
 *  @return		Median statistical error
 */ 
double MedianSError(const size_t n,
		    const double& density);

/** @brief Computes the Lebesgues width for a given metric (RMS = L2)
 *
 *  @param	vx	Array of x
 *  @param	p	Parameter p of the Lp space
 *
 *  @return		Sample mode
 */
double LWidth(const std::vector<double>& vx,
	      const double& p);

/** @brief Computes the quantile of a variable sample at probability alpha
 *
 *  @param	vx	Array of x
 *  @param	alpha	Probability value, at which the quantile is computed
 *
 *  @return		Quantile at probability alpha
 */
double Quantile(const std::vector<double>& vx,
		const double& alpha);

/** @brief Computes the statistical error on the quantile of a variable sample at probability alpha
 *
 *  @param	n	Total number of samples points
 *  @param	alpha	Probability value, at which the quantile is computed
 *  @param	density	Density at the point in question (essential)
 *
 *  @return		Quantile statistical error at probability alpha
 */
double QuantileSError(const size_t n,
		      const double& alpha,
		      const double& density);

/** @brief Computes the most probable value of the variable sample
 *
 *  @param	vx	Array of x
 *
 *  @return		Sample mode
 */
double Mode(const std::vector<double>& vx);

/** @brief Returns the minimum value in the array
 *
 *  @param	vx	Array of x
 *
 *  @return		Sample minimum
 */
double Min(const std::vector<double>& vx);

/** @brief Returns the maximum value in the array
 *
 *  @param	vx	Array of x
 *
 *  @return		Sample maximum
 */
double Max(const std::vector<double>& vx);

/** @brief Returns the Mahalanobis distance of a point to a distribution
 *
 *  @param 	invcov	Inverse covariance matrix of the data
 *  @param	vec	Vector of measurements
 *
 *  @return		Mahalanobis distance
 */
double Mahalanobis(const Matrix<double>& invcov,
		   const Matrix<double>& vec);

/** @brief Returns the squared Mahalanobis distance of a point to a distribution
 *
 *  @param 	invcov	Inverse covariance matrix of the data
 *  @param	vec	Vector of measurements
 *
 *  @return		Mahalanobis squared distance
 */
double MahalanobisSquared(const Matrix<double>& invcov,
			  const Matrix<double>& vec);

/** @brief Returns the generalised distance of a point to a distribution
 *
 *  @param	rhoR	Density of the distribution at the level of the point
 *  @param	rho0	Maximum density of the distribution
 *
 *  @return		Generalised distance
 */
double GeneralisedDistance(const double& rhoR,
			   const double& rho0);

/** @brief Returns the squared generalised distance of a point to a distribution
 *
 *  @param	rhoR	Density of the distribution at the level of the point
 *  @param	rho0	Maximum density of the distribution
 *
 *  @return		Generalised squared distance
 */
double GeneralisedDistanceSquared(const double& rhoR,
			 	  const double& rho0);

/** @brief Returns the factorial of the provided integer
 *
 *  @param	n	Integer
 *
 *  @return		Factorial
 */
size_t Factorial(const size_t n);

/** @brief Trim the n% lowest and n% largest values
 *
 *  @param	vx	Input sample to trim
 *  @param	frac	Fraction of lowest and highest values to discard
 */
template<typename T> void Trim(std::vector<T>& vx,
			       const double& frac);

/** @brief Trim the n% lowest and n% largest values
 *
 *  @param	vx	Input sample to trim
 *  @param	frac	Fraction of lowest and highest values to discard
 *
 *  @return		Trimmed sample
 */
template<typename T> std::vector<T> Trimmed(std::vector<T> vx,
			       		    const double& frac);

/** @brief Winsorize to keep only a fraction frac of the original sample
 *
 *  @param	vx	Input sample to winsorize
 *  @param	frac	Fraction of lowest and highest values to winsorize
 */
template<typename T> void Winsorize(std::vector<T>& vx,
			       	    const double& frac);

/** @brief Winsorize to keep only a fraction frac of the original sample
 *
 *  @param	vx	Input sample to winsorize
 *  @param	frac	Fraction of lowest and highest values to winsorize
 *
 *  @return		Winsorized sample
 */
template<typename T> std::vector<T> Winsorized(std::vector<T> vx,
			       	    	       const double& frac);

/** @brief Randomly select a single element from an array
 *
 *  @param	begin	Iterator to the first element
 *  @param	end	Iterator to the last element
 *
 *  @return		Iterator to the randomly selected element
 */
template<typename I> I RandomElement(I begin,
				     I end);

/** @brief Randomly resample n data points into m data points with or without replacement
 *
 *  @param	vx	Input sample
 *  @param	m	Amount of points in the subsample
 *  @param	rep	With replacement if true
 *
 *  @return		Subsample of m data points
 */
template<typename T> std::vector<T> Resample(const std::vector<T>& vx,
					     size_t m=0,
					     const bool rep=true);

#include "Statistics-inl.hh"

} // namespace Math

#endif
