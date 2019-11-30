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

// Additional modules
#include "ProgressBar.hh"
#include "Globals.hh"
#include "Statistics.hh"
#include "Exception.hh"
#include "Pitch.hh"
#include "MiceTrack.hh"

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

  /** @brief Electron peak TOF PID constructor, sets the boundaries w.r.t the electron peak
   *
   *  @param	data_files	List of file names to be processed
   *  @param	ids		IDs of the TOF stations to take into account
   *  @param	mumin		Lower boundary of the muon peak w.r.t to the electron peak
   *  @param	mumax		Upper boundary of the muon peak w.r.t to the electron peak
   *  @param	nentries	Number of entries required to fit the peak
   **/
  ParticleIdentification(const std::vector<std::string>& data_files,
			 const std::vector<size_t>& ids,
			 const double& mumin,
			 const double& mumax,
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
  void InitializeTOF(std::vector<std::string> data_files,
		     size_t nentries);

  /** @brief Initializes the TOF particle identificator manually
   *
   *  @param	data_files	List of file names to be processed
   *  @param	mumin		Lower boundary of the muon peak
   *  @param	mumax		Upper boundary of the muon peak
   *  @param	nentries	Number of entries required to fit the peak
   **/
  void InitializeTOFPeak(std::vector<std::string> data_files,
			 const double& mumin,
			 const double& mumax,
			 size_t nentries);

  /** @brief Returns the ID of the particle using whichever method was initialized
   *
   *  @param	fom		Figure of merit for pid
   *
   *  @return			Geant4 particle ID
   **/
  int GetID(const double& fom) const;

  /** @brief Returns the ID of the particle using the TOF identification method
   *
   *  @param	tof		Time-of-flight
   *
   *  @return			Geant4 particle ID
   **/
  int GetTofID(const double& tof) const;

  /** @brief Returns the distance between the two time-of-flight counters **/
  const double& GetTofDistance() const			{ return _tof_dz; }

 private:

  /** @brief Returns the ID of the particle using the TOF identification method
   *
   *  @param	data_files	List of file names to be processed
   *  @param	min		Lower boundary of the TOFs
   *  @param	max		Upper boundary of the TOFs
   *  @param	nentries	Number of entries requested
   *
   *  @return			Array of time-of-flights
   **/
  std::vector<double> TOFs(const std::vector<std::string>& data_files,
			   const double& min,
			   const double& max,
			   const size_t nentries);

  std::string		_method;	///< Name of the method used for PID
  std::vector<size_t>	_tof_ids;	///< IDs of the TOFs used for PID
  std::vector<double>	_tof_params;	///< Parameters of the TOF PID
  double		_tof_dz;	///< Distance between the TOFs
};
