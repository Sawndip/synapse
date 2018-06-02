// Cpp includes
#include <ctime>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>

// Root includes
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TNtuple.h"
#include "THStack.h"

// MAUS includes
#include "DataStructure/JobHeader.hh"
#include "DataStructure/JobHeaderData.hh"
#include "DataStructure/Spill.hh"
#include "DataStructure/Data.hh"
#include "DataStructure/MCEvent.hh"
#include "DataStructure/VirtualHit.hh"
#include "DataStructure/TOFEvent.hh"
#include "DataStructure/SciFiEvent.hh"

// Additional modules
#include "ProgressBar.hh"
#include "ParticleIdentification.hh"
#include "Aperture.hh"
#include "GeometryHandler.hh"
#include "Pitch.hh"
#include "Globals.hh"
#include "InfoBox.hh"
#include "DGaus.hh"
#include "Statistics.hh"

/** @file  ImportMAUSMinimal.cc
 *
 *  @brief Imports the minimal MAUS simulation
 *
 *	   Algorithm that imports the minimal MAUS simulation output to a single TNtuple:
 *	    - Truth: true phase space information of each track at each virtual plane.
 *	   This is used for simulations that only include virtual planes and do not require cuts.
 **/

/** @brief	Main function
 *
 *  @param	argc		Number of command line arguments
 *  @param	argv		List of command line aguments
 */
int main(int argc, char ** argv) {

  // Reroot the info messages to the terminal for now
  Pitch::setAnOutput(Pitch::info, std::cerr);

  // Importer global parameters (stored in globals)
  Globals &globals = Globals::GetInstance(argc, argv);
  if ( !globals.GetDataFiles().size() ) {
    Pitch::print(Pitch::error, "No data file specified, please execute in the following fashion:");
    Pitch::print(Pitch::error, "./import_minimal [options] data0.root [... dataN.root]");
    return 2;
  } else {
    for (const std::string& file : globals.GetDataFiles())
      if ( std::string(file).find(".root") > std::string(file).size() ) {
        Pitch::print(Pitch::error, "Invalid data file type: "+std::string(file));
	return 2;
      }
  }

  // Initialize the geometry handler with the relevant detectors
  std::vector<std::string> dets = {"tku", "tkd"};
  GeometryHandler geoh(globals["geometry_filename"], dets);

  // Information relevant to the importer
  std::string maus_version = "3.0.2";

  // Loop over the input data files, extract the relevant to the beam evolution (x, y, px, py, pz)
  std::string out_file = globals["import_filename"];
  TFile *out = new TFile(out_file.c_str(), "RECREATE");
  TNtuple* truth_samples = new TNtuple("truth_ntuple", "",
	"SpillID:EventID:VirtualPlaneID:x:y:z:px:py:pz");

  for (const std::string& file : globals.GetDataFiles()) {
    // Set up the ROOT file
    TFile data_file(file.c_str());		// Load the MAUS output file
    if ( !data_file.IsOpen() ) {
      Pitch::print(Pitch::error, "Data file not found: "+file);
      return 2;
    }
    Pitch::print(Pitch::info, "Processing "+file);

    TTree *T0 = (TTree*)data_file.Get("JobHeader");		// Pull out the TTree
    if ( !T0 ) {
      Pitch::print(Pitch::warning, "Could not find job header, cannot get MAUS version");
    } else {
      MAUS::JobHeaderData *job_ptr = new MAUS::JobHeaderData();	// A variable to store the ptr
      T0->SetBranchAddress("job_header", &job_ptr);  		// Set the address of job_ptr
      T0->GetEntry();
      if ( !maus_version.size() ) {
        maus_version = job_ptr->GetJobHeader()->GetMausVersionNumber();
        Pitch::print(Pitch::info, maus_version);
      } else if ( job_ptr->GetJobHeader()->GetMausVersionNumber() != maus_version ) {
        Pitch::print(Pitch::warning, "MAUS versions don't match, shouldn't proceed");
      }
    }

    // Set up the data pointer
    TTree *T = (TTree*)data_file.Get("Spill");	// Pull out the TTree
    MAUS::Data *data_ptr = new MAUS::Data(); 	// A variable to store the ptr to each spill
    T->SetBranchAddress("data", &data_ptr);  	// Set the address of data_ptr

    // Loop over the spills, get the mc events
    ProgressBar pbar;
    for (size_t i = 0; i < (size_t)T->GetEntries(); ++i) {

      // Display the progress in %
      pbar.GetProgress(i, (size_t)T->GetEntries());

      // Update the spill pointed to by data_ptr
      T->GetEntry(i);  
      MAUS::Spill* spill = data_ptr->GetSpill();  
      if (spill == NULL || !(spill->GetDaqEventType() == "physics_event"))
	  continue;

      // Loop over MC events in spill if provided and requested 
        std::vector<MAUS::MCEvent*>* mcevts = spill->GetMCEvents();
        for ( size_t ev = 0; ev < mcevts->size(); ++ev ) {
          if ( !mcevts->at(ev) )
	      continue;

          // Loop over the virtual hits, identify the virtual plane ids
          std::vector<MAUS::VirtualHit>* vhits = mcevts->at(ev)->GetVirtualHits();
	  MAUS::VirtualHit vhit;
	  size_t vid;
	  MAUS::ThreeVector pos, mom;
          for ( size_t vh = 0; vh < std::min((size_t)globals["nvirtuals"], vhits->size()); vh++ ) {
	    vhit = vhits->at(vh);
	    vid = vhit.GetStationId()-1;
	    pos = vhit.GetPosition();
	    mom = vhit.GetMomentum();

	    // Fill the virtual samples
	    truth_samples->Fill(i, ev, vid, pos.x(), pos.y(), pos.z(), mom.x(), mom.y(), mom.z());
	  }
        } // End of the list of MCEvents

    } // End of the list of TTree entries
    data_file.Close();
  } // End of the list of files

  // Fill the output file with the Ntuples and information about the sample
  out->cd();
  if ( truth_samples->GetEntries() )
      truth_samples->Write("Truth");
  delete truth_samples;

  TNamed("MausVersion", maus_version.c_str()).Write();
  out->Close();
}
