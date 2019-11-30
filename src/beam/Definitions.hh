#ifndef DEFINITIONS_HH
#define DEFINITIONS_HH

// Cpp includes
#include <iostream>
#include <string>
#include <map>
#include <cmath>

/** @brief Environment of beam related objects and methods */
namespace Beam {

/** @brief Enumerates the phase space variables 
 *
 *	   The axes are set up so that the \f$\vec{z}\f$ axis points downstream, the
 *	   \f$\vec{y}\f$ axis point upwards from the ground and the \f$\vec{x}\f$ axis
 *	   follows the right hand convetion, i.e. \f$ \vec{x} = \vec{y}\times\vec{z}\f$.
 *	   Coordinates :
 *	     - x: position along the transverse \f$\vec{x}\f$ axis expressed in [mm]
 *	     - x: position along the transverse \f$\vec{y}\f$ axis expressed in [mm]
 *	     - px: momentum along the transverse \f$\vec{x}\f$ axis expressed in [MeV/c]
 *	     - py: momentum along the transverse \f$\vec{y}\f$ axis expressed in [MeV/c]
 *	     - pz: momentum along the longitudinal \f$\vec{z}\f$ axis expressed in [MeV/c]
 */
enum 	PhaseSpaceVariable {x, px, y, py, pz};

/** @brief Enumerates the types of summary statistics of a beam bunch
 *
 *  	   An enumeration allows to distinguish between different summary statistics
 *   	    - alpha: Twiss parameter alpha
 *   	    - beta: Twiss parameter beta
 *   	    - gamma: Twiss parameter gamma
 *   	    - mecl: Mechanical angular momentum
 *   	    - eps: Geometric emittance
 *   	    - neps: Normalised emittance
 *   	    - trans: Tansmission
 *   	    - subeps: &alpha;-subemittance
 *   	    - amp: &alpha;-amplitude
 *   	    - vol: &alpha;-fractional emittance
 *	    - mode: amplitude mode
 *	    - mom: mean momentum
 *	    - trans: transmission
 *	    - disp: dispersion
 */
enum SumStat {alpha, beta, gamma, mecl, eps, neps, amp, subeps, vol, den, mom, trans, disp};

/** @brief Structure that stores the characteristics of each summary statistic */
struct SumStatStruct {
  std::string name;	///< Name tag of the statistic
  std::string title;	///< Long name of the statistic
  std::string label;	///< Symbolic representation of the statistic (TLatex)
  std::string unit;	///< Units in which the statistic is expressed
  bool frac;		///< Whether or not it is a fractional quantity
};

/** @brief Dictionary that maps each summary statistic to its characteristics **/
extern std::map<SumStat, SumStatStruct> SumStatDict;

} // namesapce Beam

#endif
