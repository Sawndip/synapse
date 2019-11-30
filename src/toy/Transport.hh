#ifndef TRANSPORT_HH
#define TRANSPORT_HH

// C++ includes
#include <iostream>
#include <ctime>
#include <cmath>
#include <map>

// Additional includes
#include "Matrix.hh"

/** @brief Performs the linear transport of a beam through beam line elements
 *
 *	  First build the object by providing it with the dimensionality of the
 *	  space to transport, then add elements to the beam line*
 *
 *	  The transport is operated by applying transfer matrices in succession:
 *	   - \f$X' = MX\f$;
 *	   - \f$\Sigma' = M\Sigma M^{T}\f$;
 *	   - \f$M = M_n\cdots M_1\f$.
 */
class Transport {
 public:
  /** @brief Default constructor, initilizes everything to 0 */
  Transport();

  /** @brief Dimensional constructor, sets the dimension of the space
   *
   *  @param	dim	Dimension of the space
   */
  Transport(const size_t dim);

  /** @brief Copy constructor */
  Transport(const Transport& transport);

  /** @brief Equality operator */
  Transport& operator=(const Transport& transport);

  /** @brief Destructor */
  ~Transport();

  /** @brief Performs the Transport on a particle 
   *
   *  @param	particle	Input particle: x, (y), px, (py), pz [mm, MeV/c]
   */
  void TransportParticle(Matrix<double>& particle);

  /** @brief Performs the transport on a beam of particles 
   *
   *  @param	beam	Input beam: x, (y), px, (py), pz [mm, MeV/c]
   */
  void TransportBunch(std::map<std::string, std::vector<double>>& beam);

  /** @brief Adds a drift space to the chain of transport
   *
   *  @param	L	Length of the drift space
   *  @param	p	Nominal longitudinal momentum
   */
  void AddDriftSpace(const double& L, const double& p);

  /** @brief Adds a constant solenoid to the chain of transport
   *
   *  @param	L	Length of the drift space
   *  @param	K	Ratio B_0/(2Brho_0) of the nominal field of the magnetic rigity of the ref.
   */
  void AddSolenoid(const double& L, const double& K);

  /** @brief Adds a constant non-linear solenoid to the chain of transport
   *
   *  @param	L	Length of the drift space
   *  @param	K	Ratio B_0/(2Brho_0) of the nominal field of the magnetic rigity of the ref.
   */
  void AddSolenoidNL(const double& L, const double& K);

  /** @brief Gets the total transfer matrix resulting from the 
   *
   *  @return		Transfer matrix M
   */
  const Matrix<double>& TransferMatrix() const	{ return _M; }

 private:

  size_t			_dim;	///< Dimension of the space
  Matrix<double> 		_M;	///< Transfer matrix
  std::vector<Matrix<double>>	_T;	///< Second order transfer 3-tensor
};

#endif
