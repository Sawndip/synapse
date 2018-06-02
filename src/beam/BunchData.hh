#ifndef BUNCHDATA_HH
#define BUNCHDATA_HH

// C++ includes
#include <map>
#include <vector>

// Other includes
#include "Statistics.hh"
#include "Definitions.hh"

/** @brief Environment of beam related objects and methods */
namespace Beam {

/** @brief Strutcure of a bunch map (maps variable tags to vectors) */
typedef std::map<std::string, std::vector<double>> 	BunchMap;

/** @brief Holds a the underlying phase space vectors that form a Bunch
 *
 *	   It stores the array of phase space vectors \f$(x, p_x, y, p_y, p_z)\f$
 *	   of all the particle that are contained in the bunch.
 */
class BunchData {
 public:
  /** @brief Default constructor, initilizes everything to 0 */
  BunchData();

  /** @brief Map constructor, sets the particle sample from maps
   *
   *  @param	samples		Dictionary that maps each variable onto a vector of values
   *  @param	errors		Dictionary that maps each variable onto a vector of uncertainties
   */
  BunchData(const BunchMap& samples,
            const BunchMap& errors);

  /** @brief Copy constructor */
  BunchData(const BunchData& bunchdata);

  /** @brief Equality operator */
  BunchData& operator=(const BunchData& bunchdata);

  /** @brief Destructor */
  ~BunchData();

  /** @brief Overload the subscript operator, return the samples for the requested variable */
  const std::vector<double>& operator[](const std::string& var) const	{ return Samples(var); }

  /** @brief Returns the matrix of transverse samples */
  const Matrix<double>& PhaseSpace() const				{ return _trans; }

  /** @brief Returns the matrix of uncertainties on the transverse samples */
  const Matrix<double>& PhaseSpaceErrors() const			{ return _etrans; }

  /** @brief Returns the array of samples of variable var */
  const std::vector<double>& Samples(const std::string var) const;

  /** @brief Returns the array of uncertainties on the samples of variable var */
  const std::vector<double>& Errors(const std::string var) const;

  /** @brief Returns one element of the covariance matrix of the bunch */
  const Matrix<double>& S() const					{ return _S; }

  /** @brief Returns the phase space covariance matrix of the bunch */
  const double& S(const std::string vara, const std::string varb) const;

  /** @brief Returns the matrix of transverse optical samples */
  const Matrix<double>& TraceSpace() const				{ return _opt; }

  /** @brief Returns the matrix of uncertainties on the transverse optical samples */
  const Matrix<double>& TraceSpaceErrors() const			{ return _eopt; }

  /** @brief Returns the array of optical samples of variable var */
  const std::vector<double>& Opticals(const std::string var) const;

  /** @brief Returns the array of uncertainties on the optical samples of variable var */
  const std::vector<double>& OpticalErrors(const std::string var) const;

  /** @brief Returns one element of the covariance matrix of the bunch */
  const Matrix<double>& T() const					{ return _T; }

  /** @brief Returns the phase space covariance matrix of the bunch */
  const double& T(const std::string vara, const std::string varb) const;

  /** @brief Returns the dimension of space in which the bunch lives */
  const size_t Size() const						{ return _trans.Ncols(); }

  /** @brief Returns the list of transverse axes */
  const std::vector<std::string>& Axes() const				{ return _axes; }

  /** @brief Returns the list of BunchMap dictionary keys */
  const std::vector<std::string>& Keys() const				{ return _keys; }

  /** @brief Returns the position of the requested key */
  const size_t& Id(const std::string& key) const			{ return _ids.at(key); }

  /** @brief Checks if the variable is provided */
  bool Contains(const std::string& var) const;

  /** @brief Checks if the variable is provided, throws if not */
  void AssertContains(const std::string& var) const;

 private:

  /** @brief Normal initializer, sets the particle sample and bunch variables
   *
   *  @param	samples		Dictionary that maps each variable onto a vector of values
   *  @param	errors		Dictionary that maps each variable onto a vector of errors
   */
  void Initialize(const BunchMap& samples,
		  const BunchMap& errors);

  Matrix<double>		_trans;		///< Transverse phase space
  Matrix<double>		_etrans;	///< Transverse phase space uncertainties
  Matrix<double>		_S;		///< Transverse phase space covariance matrix
  Matrix<double>		_opt;		///< Transverse trace space
  Matrix<double>		_eopt;		///< Transverse trace space uncertainties
  Matrix<double>		_T;		///< Transverse trace space covariance matrix
  std::vector<double>		_pz;		///< Vector of longitudinal momenta
  std::vector<double>		_epz;		///< Vector of longitudinal momenta uncertainties
  std::vector<std::string>	_axes;		///< List of transverse axes
  std::vector<std::string>	_keys;		///< List of transverse variables
  std::map<std::string, size_t>	_ids;		///< Reciprocal of BunchMap keys, return id of var
};
} // namespace Beam

#endif
