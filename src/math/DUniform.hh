#ifndef DUNIFORM_HH
#define DUNIFORM_HH

// C++ includes
#include <iostream>
#include <vector>
#include <cmath>

// ROOT includes
#include "TBox.h"

// Other includes
#include "DFunction.hh"

/** @brief Defines the methods associated with the Uniform distribution in any dimension
 */
class DUniform : public DFunction {
 public:

  /** @brief Basic constructor, initilizes a uniform distribution of dimension n
   *
   *  @param 	n	Dimension of the space
   */
  DUniform(const size_t n=1);

  /** @brief 1D constructor, initilizes a uniform distribution of boundaries l and u
   *
   *  @param 	l	Lower bound
   *  @param 	u	Upper bound
   */
  DUniform(const double l,
	   const double u);

  /** @brief nD constructor, initilizes a uniform distibution of given boundaries
   *
   *  @param 	l	Vector of lower bounds
   *  @param 	u	Vector of upper bounds
   */
  DUniform(const std::vector<double>& l,
	   const std::vector<double>& u);

  /** @brief Copy constructor */
  DUniform(const DUniform& df);

  /** @brief Equality operator */
  DUniform& operator=(const DUniform& df);

  /** @brief Destructor */
  ~DUniform();

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

  /** @brief Square that encompasses a certain p-value in 2 dimensions
   *
   *  @param 	alpha	P-value
   */
  std::vector<TObject*> Contour2D(const double alpha) const;

  /** @brief Cube that encompasses a certain p-value in 3 dimensions
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

  /** @brief Return the value of the Uniform distribution for an R following the Linf symmetry
   *
   *  @param	R	Radius to evaluate at (Linf symmetry)
   */
  double Radial(const double R) const;

  /** @brief Returns the radius of the alpha-contour, following an Linf symmetry
   *
   *  @param 	alpha	P-value
   */
  double Radius(const double alpha) const;

  /** @brief Returns a random value sampled from the uniform distribution */
  double Random();

  /** @brief Returns a random vector sampled from the uniform distribution */
  std::vector<double> RandomVector();

 private:

  /** @brief Initializes the distribution */
  bool Initialize();

  /** @brief Returns the parameters necessary to produce a ROOT 3D contour
   *
   *  @param	C	Level of the contour
   */
  std::vector<double> OffsetParameters(const double C) const;

  std::vector<double>	_l;		///< Lower bounds
  std::vector<double>	_u;		///< Upper bounds
  double		_vol;		///< Volume of the n-orthotope
};

#endif
