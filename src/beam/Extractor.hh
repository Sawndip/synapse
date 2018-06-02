#ifndef EXTRACTOR_HH
#define EXTRACTOR_HH

// C++ includes
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <regex>
#include <sys/stat.h>
#include <algorithm>

// ROOT includes
#include "TFile.h"
#include "TNtuple.h"
#include "TKey.h"

// Additional includes
#include "Exception.hh"
#include "Pitch.hh"
#include "ProgressBar.hh"
#include "Bunch.hh"
#include "Stream.hh"

/** @brief Environment of beam related objects and methods */
namespace Beam {

/** @brief Name of the requested data type in the imported ROOT files */
static 	std::map<std::string, std::string> 	DataTypeName =  {{"data", "Data"},
								 {"recmc", "RecMC"},
							 	 {"truth", "Truth"},
							 	 {"utruth", "UncutTruth"}};

/** @brief Extracts the data from an array of imported ROOT file
 *
 *	   Stores the list of data files to be analyzed by phase space evolution algorithms.
 *	   Upon request, loops over the list the list of files provided and looks for
 *	   the requested data type. 
 *
 *	   May be requested to limit the extraction to a set of virtual planes in
 *	   the simulation or a set of stations in the data.
 *
 *
 */
class Extractor {
 public:
  /** @brief Default constructor, initilizes everything to 0 */
  Extractor();

  /** @brief Data extractor constructor, sets the list of data files. 
   *
   *  @param	data_files	Path to the ROOT imported data files
   */
  Extractor(const std::vector<std::string>& data_files);

  /** @brief Copy constructor */
  Extractor(const Extractor& ext);

  /** @brief Equality operator */
  Extractor& operator=(const Extractor& ext);

  /** @brief Destructor */
  ~Extractor();

  /** @brief Gets the value of the requested dictionary entry */
  const std::string& Get(const std::string& key) const	{ return _dictionary.at(key); }

  /** @brief Gets the list of keys in the dictionary */
  const std::vector<std::string> GetListOfKeys() const;

  /** @brief Extracts the beam at a specific z position
   *
   *  @param	data_type	Type of data
   *  @param	plane_id	Requested plane ID
   *  @param	measerr		Extract the measurement errors
   *
   *  @return 			Bunch object at the specified plane ID
   */
  Bunch GetBunch(const std::string& data_type,
	         const size_t plane_id,
	         const bool measerr=false) const;

  /** @brief Extracts beam at every requested z position
   *
   *  @param	data_type	Type of data
   *  @pram	plane_ids	Requested plane IDs. if empty, extracts everything
   *  @param	measerr		Extract the measurement errors
   *
   *  @return 			Stream object for the selected plane IDs
   */
  Stream GetStream(const std::string& data_type,
		   const std::vector<size_t> plane_ids=std::vector<size_t>(),
		   const bool measerr=false) const;

  /** @brief Extracts beam at every requested z position and for the requested data types
   *
   *  @param	data_types	Types of data
   *  @pram	plane_ids	Map of requested plane IDs. If empty, extracts everything
   *  @param	measerr		Extract the measurement errors
   *
   *  @return 			Stream object for the selected plane IDs
   */
  std::map<std::string, Stream>
	GetStreams(const std::vector<std::string>& data_types,
		   const std::map<std::string, std::vector<size_t>> plane_ids =
			 std::map<std::string, std::vector<size_t>>(),
		   const bool measerr=false) const;

  /** @brief Sets the run name */
  void SetRunName(const std::string& run_name)	{ _run_name = run_name; }

  /** @brief Returns the run name */
  std::string GetRunName() const		{ return _run_name; }

 private:

  /** @brief Initializes the extractor 
   *
   *	     Check that they exist and are of the right type.
   */
  void Initialize();

  /** @brief Checks if file exists at the specified path */
  bool Exists(const std::string& file) const;

  std::vector<std::string>		_data_files;	///< List of data files
  std::string				_run_name;	///< Name of the run
  std::map<std::string, std::string>	_dictionary;	///< Dictionary of run characteristics
};
} // namespace Beam

#endif
