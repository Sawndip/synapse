#ifndef DFUNCTION_HH
#define DFUNCTION_HH

// C++ includes
#include <iostream>

// ROOT includes
#include "TStyle.h"
#include "TObject.h"
#include "TLine.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TF3.h"
#include "TRandom3.h"
#include "TCanvas.h"

// Other includes
#include "Pitch.hh"
#include "Exception.hh"
#include "Assert.hh"

/** @brief Parent class for all the probability density functions D*.
 *
 *          Sets up common and virtual functions and members that may be called
 *	    by all of the child function classes
 */
class DFunction {
 public:
  /** @brief Default constructor, initilizes everything to 0 */
  DFunction();

  /** @brief Copy constructor */
  DFunction(const DFunction& df);

  /** @brief Equality operator */
  DFunction& operator=(const DFunction& df);

  /** @brief Overloaded function operator, calls Evaluate
   *
   *  @param	v	Vector of coordinates
   */
  double operator()(const std::vector<double>& v) const	{ return Evaluate(v); }

  /** @brief Overloaded function operator, calls Evaluate
   *
   *  @param	v	One dimensional value
   */

  double operator()(const double& v) const		{ return Evaluate(v); }

  /** @brief Overloaded function operator, calls Evaluate
   *
   *  @param	v	Pointer to an array of coordinates
   */
  double operator()(const double* v) const		{ return Evaluate(v); }

  /** @brief Destructor */
  ~DFunction();

  /** @brief Cumulative distribution function at v
   *
   *  @param 	v	Value in which to evaluate the CDF
   */
  virtual double CDF(const double& v) const = 0;

  /** @brief If the function has an L1, L2 or Linf symmetry, returns its CDF at radius R
   *
   *  @param	R	Radius to evaluate at
   */
  virtual double CDFRadial(const double& R) const = 0;

  /** @brief Segments that encompass a certain p-value in 1 dimension
   *
   *  @param 	alpha	P-value
   */
  virtual std::vector<TLine*> Contour1D(const double alpha) const = 0;

  /** @brief Geometrical objects that encompass a certain p-value in 2 dimensions
   *
   *  @param 	alpha	P-value
   */
  virtual std::vector<TObject*> Contour2D(const double alpha) const = 0;

  /** @brief TF3 contour that encompasses a certain p-value in 3 dimensions
   *
   *  @param 	alpha	P-value
   */
  virtual TF3* Contour3D(const double alpha) const = 0;

  /** @brief Segments that encompass a certain p-value in radius
   *
   *  @param 	alpha	P-value
   */
  virtual std::vector<TLine*> ContourRadial(const double alpha) const = 0;

  /** @brief Volume of the smallest contour of a given p-value
   *
   *  @param 	alpha	P-value
   */
  virtual double ContourVolume(const double alpha) const = 0;

  /** @brief Returns the dimension of the space */
  size_t Dimension() const				{ return _dim; }

  /** @brief Returns the value of the PDF
   *
   *  @param	v	Vector of coordinates
   */
  virtual double Evaluate(const std::vector<double>& v) const = 0;

  /** @brief Returns the value of the PDF
   *
   *  @param	v	One dimensional value
   */
  double Evaluate(const double& v) const;

  /** @brief Returns the value of the PDF
   *
   *  @param	v	Pointer to an array of coordinates
   */
  double Evaluate(const double* v) const;

  /** @brief Returns the level of the alpha-contour
   *
   *  @param 	alpha	P-value
   */
  virtual double Level(const double alpha) const = 0;

  /** @brief Sets the lower bound of axis i */
  void SetLowerBound(const double lower, const size_t i=0);

  /** @brief Returns the lower bound of axis i */
  double LowerBound(const size_t i=0) const		{ return _lower[i]; }

  /** @brief Sets the vector of lower bounds */
  void SetLowerBounds(std::vector<double> lower)	{ _lower = lower; }

  /** @brief Returns the array of lower bounds */
  std::vector<double> LowerBounds() const		{ return _lower; }

  /** @brief Returns the name of the function */
  void SetName(const std::string& name) 		{ _name = name; }

  /** @brief Returns the name of the function */
  std::string Name() const				{ return _name; }

  /** @brief Returns the probability content within the defined interval */
  virtual double Norm() const = 0;

  /** @brief If the function has an L1, L2 or Linf symmetry, returns its value at radius R
   *
   *  @param	R	Radius to evaluate at
   */
  virtual double Radial(const double R) const = 0;

  /** @brief Returns a random value sampled from the distribution */
  virtual double Random() = 0;

  /** @brief Returns a random vector sampled from the distribution */
  virtual std::vector<double> RandomVector() = 0;

  /** @brief Returns the upper bound of axis i */
  void SetRange(const double lower,
		const double upper,
		const size_t i=0);

  /** @brief Returns the upper bound of axis i */
  void SetRange(const std::vector<double>& lower,
		const std::vector<double>& upper);

  /** @brief Returns the upper bound of axis i */
  void Range(std::vector<double>& lower,
	     std::vector<double>& upper) const;

  /** @brief Returns the name of the function */
  void SetTitle(const std::string& title) 		{ _title = title; }

  /** @brief Returns the name of the function */
  std::string Title() const				{ return _title; }

  /** @brief Sets the upper bound of axis i */
  void SetUpperBound(const double upper, const size_t i=0);

  /** @brief Returns the upper bound of axis i */
  double UpperBound(const size_t i=0) const		{ return _upper[i]; }

  /** @brief Sets the vector of upper bounds */
  void SetUpperBounds(std::vector<double> upper)	{ _upper = upper; }

  /** @brief Returns the array of upper bounds */
  std::vector<double> UpperBounds() const		{ return _upper; }

  /** @brief Returns a 1D graph of the function in the requested axis
   *
   *  @param	idx	ID of the first coordinate to represent
   *  @param	x	Point included in the line (fixes the remaining coordinates)
   */
  TGraph* Graph(size_t idx=0, std::vector<double> x=std::vector<double>()) const;

  /** @brief Returns a radial graph of the function */
  TGraph* GraphRadial() const;

  /** @brief Returns a 2D histogram of the function in the requested axis
   *
   *  @param	idx	ID of the first coordinate to represent
   *  @param	idy	ID of the second coordinate to represent
   *  @param	x	Point included in the plane (fixes the remaining coordinates)
   */
  TH2F* Graph2D(size_t idx=0, size_t idy=1, std::vector<double> x=std::vector<double>()) const;

  /** @brief Returns a radial graph of the cummulative density function */
  TGraph* GraphCDFRadial() const;

  /** @brief Draws the function on whichever TCanvas is currently being used (ROOT) 
   *
   *  @param	opt	Drawing options
   *  @param	idx	ID of the first coordinate to represent
   *  @param	idy	ID of the second coordinate to represent
   *  @param	x	Point that the projection intersects	
   */
  void Draw(const std::string opt="", int idx=0, int idy=1,
		std::vector<double> x=std::vector<double>()) const;

  /** @brief Draws the radial function on whichever TCanvas is currently being used (ROOT)
   *
   *  @param	opt	Drawing options
   */
  void DrawRadial(const std::string opt="") const;

  /** @brief Creates a canvas, draws the function on it and save it to a file
   *
   *  @param	opt	Drawing options
   *  @param	idx	ID of the first coordinate to represent
   *  @param	idy	ID of the second coordinate to represent
   *  @param	x	Point that the projection intersects	
   */
  void Print(const std::string opt="", int idx=0, int idy=1,
		std::vector<double> x=std::vector<double>()) const;

  /** @brief Creates a canvas, draws the radial function on it and save it to a file
   *
   *  @param	opt	Drawing options
   */
  void PrintRadial(const std::string opt="") const;

 protected:

  size_t		_dim;		///< Dimension of the space the function lives in
  std::vector<double> 	_lower;		///< Array of lower bounds
  std::vector<double> 	_upper;		///< Array of upper bounds
  std::string		_name;		///< Name of the function
  std::string		_title;		///< Title of the function
  TRandom3		_rdmzer;	///< Pseudorandom number generator
};

#endif
