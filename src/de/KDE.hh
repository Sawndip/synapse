// C++ includes
#include <iostream>
#include <ctime>
#include <cmath>
#include <vector>

// Additional includes
#include "Matrix.hh"
#include "Statistics.hh"
#include "DGaus.hh"

/** @brief Computes the Kernel Density Estimator (KDE) of a set of points
 */
class KDE {
 public:
  /** @brief Default constructor, initilizes everything to 0 */
  KDE();

  /** @brief Proper constructor, initilizes everything */
  KDE(std::vector<std::vector<double>> mat);

  /** @brief Copy constructor */
  KDE(const KDE& kde);

  /** @brief Equality operator */
  KDE& operator=(const KDE& kde);

  /** @brief Overloaded function operator, calls Evaluate
   *
   *  @param	v	Vector of coordinates
   */
  double operator()(const std::vector<double>& v) const		{ return Evaluate(v); }

  /** @brief Overloaded function operator, calls Evaluate
   *
   *  @param	v	One dimensional value
   */
  double operator()(const double& v) const			{ return Evaluate(v); }

  /** @brief Overloaded function operator, calls Evaluate
   *
   *  @param	v	Pointer to an array of coordinates
   */
  double operator()(const double* v) const 			{ return Evaluate(v); }

  /** @brief Destructor */
  ~KDE();

  /** @brief Performs the density estimation at an n-point v
   *
   *  @param	v	Vector of coordinates of the n-point
   */
  double Evaluate(const std::vector<double>& v) const;

  /** @brief Performs the density estimation at an n-point v
   *
   *  @param	v	One dimensional value of the point
   */
  double Evaluate(const double& v) const;

  /** @brief Performs the density estimation at an n-point v
   *
   *  @param	v	Pointer to an array of coordinates of the n-point
   */
  double Evaluate(const double* v) const;

  /** @brief Returns the value of the KDE in the d-dimensional point (fast algorithm)
   *  @param	points	dxN matrix of points in which to evaluate the KDE
   *  @param	eps	Tolarable error with respect to the exact KDE
   *  @return		Approximated value of the KDE in all of the points (N-vector)
   */
  std::vector<double> FastFunction(std::vector<std::vector<double>> points,
				    double eps = 1e-3);

  /** @brief Returns the value of the KDE in the d-dimensional point (ROOT format)
   *  @param	x	d-array of variables
   *  @param	par	Array of parameters
   *  @return		Value of the KDE in x
   */
  double FunctionROOT(double *x, double *par);

  /** @brief Initializes the KDE parameters, bandwidth */
  void Initialize();

 private:

  /** @brief Returns the Silverman bandwidth matrix of order n 
   *  @param	vsig	d-vector of standard deviations of the variables
   *  @param	n	Number of samples
   *  @return		Diagonal Silverman bandwitdh
   */
  Matrix<double> SilvermanBandwidth(std::vector<double> vsig, size_t n);

  /** @brief Returns the Scott bandwidth matrix of order n 
   *  @param	vsig	Vector of standard deviations on the variables
   *  @param	n	Number of samples
   *  @return		Diagonal Scott bandwitdh
   */
  Matrix<double> ScottBandwidth(std::vector<double> vsig, size_t n);

  std::string 				_btype;		///< Bandwidth algorithm name
  Matrix<double>			_H;		///< nxn bandwidth matrix
  std::vector<std::vector<double>>	_samples;	///< nxN sample matrix
  DGaus					_kernel;	///< Multivariate Gaussian Kernel
};
