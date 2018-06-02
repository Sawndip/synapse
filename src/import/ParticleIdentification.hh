// C++ includes
#include <iostream>
#include <vector>
#include <string>

// ROOT includes
#include "TStyle.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TSpectrum.h"

// MAUS includes
#include "DataStructure/Spill.hh"
#include "DataStructure/Data.hh"
#include "DataStructure/TOFEvent.hh"

// Additional modules
#include "ProgressBar.hh"
#include "Statistics.hh"
#include "Exception.hh"
#include "Pitch.hh"

/** @brief Performs particle identification.
 *
 *	   Currently only one implemented identification method (TOF).
 *
 *	   Time-of-flight particle identification:
 *	    -# fills a histogram with the time-of-flights between TOF stations i and j;
 *	    -# fits the muon and pion peaks with a two-peak Gaussian;
 *	    -# determines the time limits between the e, mu and pi peaks.
 *	    -# when called, checks in which interval the time-of-flight falls.
 *
 *	    Returns the particle species as a Geant4 particle ID:
 *	     - Electron/positron : 11;
 *	     - Muon/Antimuon: 13;
 *	     - Charged pion: 211.
 **/
class ParticleIdentification {
 public:
  /** @brief Default constructor, initilizes everything to 0 **/
  ParticleIdentification();

  /** @brief TOF PID constructor, sets the distance between used TOF stations and their ids
   *
   *  @param	data_files	List of file names to be processed
   *  @param	ids		IDs of the TOF stations to take into account
   *  @param	nentries	Number of entries required to make the fit
   **/
  ParticleIdentification(const std::vector<std::string>& data_files,
			 const std::vector<size_t>& ids,
			 const size_t nentries = 1e4);

  /** @brief Manual TOF PID constructor, sets the boundaries of the muon peak manually
   *
   *  @param	ids		IDs of the TOF stations to take into account
   *  @param	mumin		Lower boundary of the muon peak
   *  @param	mumax		Upper boundary of the muon peak
   **/
  ParticleIdentification(const std::vector<size_t>& ids,
			 const double& mumin,
			 const double& mumax);

  /** @brief Copy constructor **/
  ParticleIdentification(const ParticleIdentification& pid);

  /** @brief Equality operator **/
  ParticleIdentification& operator=(const ParticleIdentification& pid);

  /** @brief Destructor **/
  ~ParticleIdentification();

  /** @brief Initializes the TOF particle identificator
   *
   *  @param	data_files	List of file names to be processed
   *  @param	nentries	Number of entries required to make the fit
   **/
  void InitializeTOF(std::vector<std::string> data_files, size_t nentries);

  /** @brief Returns the ID of the particle using whichever method was initialized
   *
   *  @param	event		MAUS reconstructed event (single trigger)
   *
   *  @return			Geant4 particle ID
   **/
  int GetID(MAUS::ReconEvent* event) const;

  /** @brief Returns the ID of the particle using the TOF identification method
   *
   *  @param	event		MAUS reconstructed event (single trigger)
   *
   *  @return			Geant4 particle ID
   **/
  int GetTofID(MAUS::ReconEvent* event) const;

 private:

  std::string		_method;	///< Name of the method used for PID
  std::vector<size_t>	_tof_ids;	///< IDs of the TOFs used for PID
  std::vector<double>	_tof_params;	///< Parameters of the TOF PID
  double		_tof_dz;	///< Distance between the TOFs
};
