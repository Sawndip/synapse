// C++ includes
#include <iostream>
#include <cmath>

/** @brief Environment of beam related objects and methods */
namespace Beam {

/** @brief Contains the value, measurement uncertainty and statistical uncertainty of a variable.
 */
class Variable {
 public:
  /** @brief Default constructor, initilizes everything to 0 */
  Variable();

  /** @brief Normal constructor, sets the value and uncertainties
   *
   *  @param	val	Value of the variable
   *  @param	merr	Measurement error on the variable
   *  @param	serr	Statistical error on the variable
   */
  Variable(const double& val,
	   const double merr=0.,
	   const double serr=0.);

  /** @brief Copy constructor */
  Variable(const Variable& var);

  /** @brief Equality operator */
  Variable& operator=(const Variable& var);

  /** @brief Destructor */
  ~Variable();

  /** @brief Overloads the double cast operator */
  explicit operator double() const;

  /** @brief Overloads the ostream operator */
  friend std::ostream& operator<<(std::ostream& os, const Variable& var);

  /** @brief Returns the value of the variable */
  const double& GetValue() const		{ return _val; }

  /** @brief Sets the value of the variable */
  void SetValue(const double& val)		{ _val = val; }

  /** @brief Returns the measurement uncertainty of the variable */
  const double& GetMError() const		{ return _merr; }

  /** @brief Sets the measurement uncertainty of the variable */
  void SetMError(const double& merr)		{ _merr = merr; }

  /** @brief Returns the value of the variable */
  const double& GetSError() const		{ return _serr; }

  /** @brief Returns the reference to the statistical error on the variable to set it */
  void SetSError(const double& serr)		{ _serr = serr; }

  /** @brief Returns the total quadratic error on the variable */
  double GetError() const			{ return sqrt(_merr*_merr+_serr*_serr); }

 private:
  double	_val;		///< Value of the variable
  double	_merr;		///< Value of the measurement uncertainty
  double 	_serr;		///< Value of the statistical uncertainty
};
} // namespace Beam
