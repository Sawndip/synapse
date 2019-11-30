#ifndef STREAM_HH
#define STREAM_HH

// C++ includes
#include <stdlib.h>
#include <map>
#include <vector>

// Other includes
#include "Statistics.hh"
#include "DGaus.hh"
#include "DensityEstimator.hh"
#include "ProbabilityContour.hh"
#include "Drawer.hh"
#include "Bunch.hh"

/** @brief Environment of beam related objects and methods */
namespace Beam {

/** @brief Holds a set of particle bunches
 *
 *	   Holds particle bunches at an array of z positions, or planes, to evaluate global
 * 	   characteristics of a stream of particles along a the full extent of a beam line.
 */
class Stream {
 public:
  /** @brief Default constructor, initilizes everything to 0 */
  Stream();

  /** @brief Map constructor, sets the particle stream from a map of bunches
   *
   *  @param	bunches		Dictionary that maps each plane id onto a Bunch
   *  @param	reaches		Reach of each particle in the beam (max plane id)
   */
  Stream(const std::map<size_t, Bunch>& bunches,
	 const std::vector<size_t> reaches = std::vector<size_t>());

  /** @brief Copy constructor */
  Stream(const Stream& stream);

  /** @brief Equality operator */
  Stream& operator=(const Stream& stream);

  /** @brief Destructor */
  ~Stream();

  /** @brief Returns the number of bunches */
  const size_t Size() const					{ return _bunches.size(); }

  /** @brief Returns the first iterator of the list of bunches */
  const std::map<size_t, Bunch>::const_iterator begin() const	{ return _bunches.begin(); }

  /** @brief Returns the last iterator of the list of bunches */
  const std::map<size_t, Bunch>::const_iterator end() const	{ return _bunches.end(); }

  /** @brief Returns the most upstream bunch in the stream */
  const Bunch& Front() const;

  /** @brief Returns the most downstream bunch in the stream */
  const Bunch& Back() const;

  /** @brief Overload the subscript operator, return a bunch */
  const Bunch& operator[](const size_t plane_id) const		{ return _bunches.at(plane_id); }

  /** @brief Overload the subscript operator, return a bunch, allows mods */
  Bunch& operator[](const size_t plane_id)			{ return _bunches[plane_id]; }

  /** @brief Returns the list of plane ids that compose the stream */
  const std::vector<size_t>& GetPlaneIds() const		{ return _plane_ids; }

  /** @brief Returns the reaches of all the particles that compose the stream */
  const std::vector<size_t>& GetReaches() const			{ return _reaches; }

  /** @brief Returns the reach of the particle of requested id */
  const size_t& GetReach(const size_t id) const			{ return _reaches[id]; }

  /** @brief Returns an the transmission in the requested plane id */
  double Transmission(const size_t plane_id) const;

  /** @brief Sets the core fraction of all the bunches in the stream */
  void SetCoreFraction(const double frac,
		       const std::string fom="amp");

  /** @brief Returns the change in the summary statistic between two planes
   *
   *  @param 	stat	Tag of the summary statistic
   *  @param 	idu	Id of the upstream plane
   *  @param 	idd	Id of the downstream plane
   *
   *  @return		Change and its uncertainty (variable)
   */
  Variable Change(const SumStat& stat,
		  const size_t idu,
		  const size_t idd) const;

  /** @brief Returns an evolution graph of the requested statistic 
   *
   *  @param 	stat	Tag of the summary statistics to be drawn
   *
   *  @return		Pointer to a TGraphErrors
   */
  TGraphErrors* EvolutionGraph(const SumStat& stat) const;

  /** @brief Returns an evolution graph of the beam transmission
   *
   *  @return		Pointer to a TGraphErrors
   */
  TGraphErrors* TransmissionGraph() const;

  /** @brief Returns the relative change in a fractional quantities as a function of the fraction
   *
   *  @brief 	idu	Id of the upstream reference plane
   *  @brief 	idd	Id of the downstream reference plane
   *
   *  @return		Map of TGraphErrors
   */
  std::map<SumStat, TGraphErrors*> FractionalGraphs(const size_t idu,
						    const size_t idd);

  /** @brief Returns a function that expresses the fractional quantify as a function of emittance */
  TF1* FractionalFunction(const SumStat& stat) const;

 private:

  /** @brief Normal initializer, sets the list global variables */
  void Initialize();

  /** @brief Sets the number of particles contained in the core of all the bunches in the stream */
  void SetCoreSize(const size_t size,
		   const std::string fom="amp");

  /** @brief Returns true if the stream contains the requested plane id */
  bool Contains(const size_t plane_id) const;

  std::map<size_t, Bunch>	_bunches;	///< Map of beam bunches that compose the stream
  std::vector<size_t>		_plane_ids;	///< List of plane ids that compose the stream
  std::vector<size_t>		_reaches;	///< Reach of each particle in the beam
  double			_fraction;	///< Fraction of the first bunch in the core
};
} // namespace Beam

#endif
