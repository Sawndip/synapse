#ifndef MICETRACK_HH
#define MICETRACK_HH

// c++ includes
#include <map>

// ROOT includes
#include "TVector3.h"

typedef std::map<size_t, TVector3>	SPArray;
typedef std::map<size_t, int>		PIDArray;

class MiceTrack {
 public:

  // Track parameters
  size_t		spillid = 0;	///< Spill ID
  size_t		eventid = 0;	///< Event ID

  // Detector counters
  std::vector<size_t>	tof_nsp = std::vector<size_t>(2, 0);	///< Number of space points in TOFs
  std::vector<size_t>	tk_nt = std::vector<size_t>(2, 0);	///< Number of tracks in trackers

  // TOF times
  std::vector<double> 	t = std::vector<double>(2, -1.);	///< Time-of-flight times

  // Quality criteria
  std::vector<double> 	tk_chi2 = std::vector<double>(2, -1.);	///< Trackers chi squared over ndf
  std::vector<double> 	tk_maxr = std::vector<double>(2, -1.);	///< Maximum radius in the trackers

  // Phase space
  SPArray		pos;		///< Array of positions in the sample planes
  SPArray		pose;		///< Array of position uncertainties in the sample planes
  SPArray		mom;		///< Array of 3-momenta in the sample planes
  SPArray		mome;		///< Array of 3-momentum uncertainties in the sample planes

  // Particle id
  PIDArray		pid;		///< Particle ID if known in the sample planes
};

#endif
