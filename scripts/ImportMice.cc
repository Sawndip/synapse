// Cpp includes
#include <iostream>
#include <vector>
#include <map>

// Root includes
#include "TString.h"
#include "TFile.h"
#include "TTree.h"

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
#include "Pitch.hh"
#include "ProgressBar.hh"
#include "Globals.hh"
#include "Aperture.hh"
#include "MiceTrack.hh"
#include "CircleFitter.hh"

/** @file  ImportMice.cc
 *
 *  @brief Imports the MICE data tree structure into a set of MICE tracks
 *
 *	   Algorithm that imports the MICE data and/or reconstructed simulation along with
 *	   the Monte Carlo truth to a set of two TTree branches made of MiceTrack objects:
 *	    - Truth: Monte Carlo truth;
 *	    - Recon: Reconstructed Monte Carlo OR Data;
 */

/** @brief	Main function
 *
 *  @param	argc		Number of command line arguments
 *  @param	argv		List of command line aguments
 */
int main(int argc, char ** argv) {

  // Reroot the info messages to the terminal for now
  Pitch::setAnOutput(Pitch::info, std::cerr);

  // Importer global parameters (stored in globals)
  // TODO make the file check a function of Globals!
  Globals &globals = Globals::GetInstance(argc, argv);
  if ( !globals.GetDataFiles().size() ) {
    Pitch::print(Pitch::error, "No data file specified, please execute in the following fashion:");
    Pitch::print(Pitch::error, "./import_data [options] data0.root [... dataN.root]");
    return 2;
  } else {
    for (const std::string& file : globals.GetDataFiles()) {
      if ( std::string(file).find(".root") > std::string(file).size() ) {
        Pitch::print(Pitch::error, "Invalid data file type: "+std::string(file));
	return 2;
      }
    }
  }

  // Information relevant to the importer
  std::string maus_version = globals["maus_version"];

  // Load mice Track library to ROOT, ignore version errors
  gSystem->Load("../lib/libMiceTrack.so");
  gErrorIgnoreLevel = kBreak;

  // Initialize a MiceTrack object for the MC truth and the reconstruction
  MiceTrack* track_tru = new MiceTrack;
  MiceTrack* track_rec = new MiceTrack;

  // Loop over the input data files, extract the data relevant to the beam evolution
  std::string out_file = globals["import_filename"];
  TFile *out = new TFile(out_file.c_str(), "RECREATE");

  size_t spilloff = 0;
  size_t tkid, planeid;
  double x0, y0, rad;

  TTree tree("Track", "Main MICE data tree");
  TBranch *branch_tru(NULL), *branch_rec(NULL);

  for (const std::string& file : globals.GetDataFiles()) {
    // Set up the ROOT file
    TFile data_file(file.c_str());
    if ( !data_file.IsOpen() ) {
      Pitch::print(Pitch::error, "Data file not found: "+file);
      return 2;
    }
    Pitch::print(Pitch::info, "Processing "+file);

    // If chunked, get the chunk number from the file name
    if ( globals["chunked"] ) {
      std::string chunk_number = std::regex_replace(file,
				  std::regex("import_|(.*/)|(_sim.root)"), std::string(""));
      spilloff = 1e3*atoi(chunk_number.c_str());
      std::cerr << spilloff << std::endl;
    }

    // Get the MAUS version from the loaded file, flag inconsistencies
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
        Pitch::print(Pitch::warning, "MAUS versions do not match, shouldn't proceed");
      }
    }

    // Set up the data pointer
    TTree *T = (TTree*)data_file.Get("Spill");	// Pull out the TTree
    MAUS::Data *data_ptr = new MAUS::Data(); 	// A variable to store the ptr to each spill
    T->SetBranchAddress("data", &data_ptr);  	// Set the address of data_ptr

    // Identify whether it is a simulation or real data
    bool simulated(false), recon(false);
    for (size_t i = 0; i < std::min((size_t)100, (size_t)T->GetEntries()); i++) {
      T->GetEntry(i);
      if ( data_ptr->GetSpill()->GetMCEvents()->size() ) {
	simulated = true;
	if ( !branch_tru )
	    branch_tru = tree.Branch("Truth", track_tru);
      }
      if ( data_ptr->GetSpill()->GetReconEvents()->size() ) {
	recon = true;
  	if ( !branch_rec )
	    branch_rec = tree.Branch("Recon", track_rec);
      }
      if ( simulated && recon )
  	  break;
    }
    if ( !simulated && !recon ) {
      Pitch::print(Pitch::warning, "Nothing found in this file, skipping");
      continue;
    }

    // Loop over the spills, get the recon and mc events
    ProgressBar pbar;
    for (size_t i = 0; i < (size_t)T->GetEntries(); ++i) {

      // Display the progress in %
      pbar.GetProgress(i, (size_t)T->GetEntries());

      // Update the spill pointed to by data_ptr
      T->GetEntry(i);  
      MAUS::Spill* spill = data_ptr->GetSpill();  
      if (spill == NULL || !(spill->GetDaqEventType() == "physics_event"))
	  continue;

      // Loop over events in the spill if provided
      std::vector<MAUS::ReconEvent*>* revts = spill->GetReconEvents();
      std::vector<MAUS::MCEvent*>* mcevts = spill->GetMCEvents();
      for ( size_t ev = 0; ev < revts->size(); ++ev ) {
        if ( recon && revts->at(ev) ) {

	  // Reset the track information
	  *track_rec = MiceTrack();

	  // Set the spill and event
	  track_rec->spillid = spilloff+i;
	  track_rec->eventid = ev;

	  // If the information is present, get the TOF01 time-of-flight
	  if ( !globals["notofs"] ) {

	    // Get the TOF space points
	    MAUS::TOFEvent* tofevent = revts->at(ev)->GetTOFEvent();
            MAUS::TOFEventSpacePoint tofsp = tofevent->GetTOFEventSpacePoint();
            std::vector<std::vector<MAUS::TOFSpacePoint>> tofsps = {tofsp.GetTOF0SpacePointArray(),
							   	    tofsp.GetTOF1SpacePointArray()};

  	    // Get the number of space points and the times in each TOF station
  	    for (size_t j = 0 ; j < 2; j++) {
	      track_rec->tof_nsp[j] = tofsps[j].size();
              track_rec->t[j] = -1.;
              if ( tofsps[j].size() == 1 )
                  track_rec->t[j] = tofsps[j][0].GetTime();
            }

	    // If there is no space point in TOF1, skip as the experiment is not triggered
	    if ( !track_rec->tof_nsp[0] )
	        continue;
	  }

	  // Fetch the kalman tracks and extract the phase space coordinates
          MAUS::SciFiEvent* sfevt = revts->at(ev)->GetSciFiEvent(); // Pull out scifi event
          MAUS::SciFiTrackPArray sfts = sfevt->scifitracks();

	  // Loop over Kalman tracks (helical only)
          for ( const MAUS::SciFiTrack *sft : sfts ) {
	    // If it is not fitted with a helix, no momentum info, not a valid track
	    if ( !sft->GetAlgorithmUsed() || !sft->pr_track_pointer_helical() )
	        continue;

	    // Get the tracker id
	    tkid = sft->tracker();

	    // Increment the counter
	    track_rec->tk_nt[tkid]++;

	    // Save the chi^2/ndf
	    track_rec->tk_chi2[tkid] = sft->chi2()/sft->ndf();

	    // Get the underlying helix from PR, record the maximum radius
	    MAUS::SciFiHelicalPRTrack* prtrack = sft->pr_track_pointer_helical();
	    x0 = prtrack->get_circle_x0();
	    y0 = prtrack->get_circle_y0();
	    rad = fabs(prtrack->get_R());
	    track_rec->tk_maxr[tkid] = sqrt(x0*x0+y0*y0)+rad;

	    // Loop over the Kalman track points
	    MAUS::SciFiTrackPointPArray tps = sft->scifitrackpoints();
	    for (const MAUS::SciFiTrackPoint *tpoint : tps ) {

	      if ( !tpoint->plane() ) {
		planeid = tkid*5 + tpoint->station()-1;
	      	track_rec->pos[planeid] = tpoint->pos();
	      	track_rec->pose[planeid] = tpoint->pos_error();
	      	track_rec->mom[planeid] = tpoint->mom();
	      	track_rec->mome[planeid] = tpoint->mom_error();
	      }
	    }
	  }
        }

	// If simulated, get the truth that corresponds to the reconstruction
	if ( simulated && mcevts->at(ev) ) {

	  // Reset the track information
	  *track_tru = MiceTrack();

	  // Set the spill and event
	  track_tru->spillid = spilloff+i;
	  track_tru->eventid = ev;

	  // Find the hit ID of the relevant beam line elements
	  std::map<std::string, int> vids =
		{{"tku",-1},{"tkd",-1},{"tkd3",-1},{"tof0",-1},{"tof1",-1},{"min",-1},{"max",-1}};
          std::vector<MAUS::VirtualHit>* vhits = mcevts->at(ev)->GetVirtualHits();
  	  MAUS::VirtualHit vhit;
	  size_t vid;
          for ( size_t vh = 0; vh < vhits->size(); vh++ ) {
	    vhit = vhits->at(vh);
	    vid = vhit.GetStationId()-1;

	    for (const std::string& det : {"tku", "tkd", "tkd3", "tof0", "tof1", "min", "max"})
	      if ( vid == (size_t)globals[det+"_vid"] )
	          vids[det] = vh;
	  }

	  // Store the true time-of-flight information
	  if ( !globals["notofs"] ) {
	    if ( vids["tof0"] > 0 )
	      track_tru->t[0] = vhits->at(vids["tof0"]).GetTime();
	    if ( vids["tof1"] > 0 )
	      track_tru->t[1] = vhits->at(vids["tof1"]).GetTime();
          }

          // Measure the true maximal radius by fitting the true tracker station space points
	  std::map<std::string, bool> fiducial;
	  MAUS::ThreeVector pos;
	  if ( vids["tku"] > 0 ) {
	    track_tru->tk_nt[0]++;
	    std::vector<double> xx, yy;
            for (size_t vh = vids["min"]; vh <= (size_t)vids["tku"]; vh++) {
  	      vhit = vhits->at(vh);
	      pos = vhit.GetPosition();
	      xx.push_back(pos.x());
	      yy.push_back(pos.y());
	    }
	    FitCircle(xx, yy, x0, y0, rad);
	    track_tru->tk_maxr[0] = sqrt(x0*x0+y0*y0)+rad;
	  }
	  if ( vids["tkd3"] > 0 ) {
	    track_tru->tk_nt[1]++;
	    size_t limit = vids["max"] > 0 ? vids["max"]+1 : vhits->size();
	    std::vector<double> xx, yy;
            for (size_t vh = vids["tkd"]; vh < limit; vh++) {
  	      vhit = vhits->at(vh);
	      pos = vhit.GetPosition();
	      xx.push_back(pos.x());
	      yy.push_back(pos.y());
	    }
	    FitCircle(xx, yy, x0, y0, rad);
	    track_tru->tk_maxr[1] = sqrt(x0*x0+y0*y0)+rad;
	  }

          for ( size_t vh = 0; vh < vhits->size(); vh++ ) {
	    vid = vhits->at(vh).GetStationId()-1;
	    track_tru->pid[vid] = vhits->at(vh).GetParticleId();
	    track_tru->pos[vid] = vhits->at(vh).GetPosition();
	    track_tru->mom[vid] = vhits->at(vh).GetMomentum();
	  }
        }

	// Fill the tree
	tree.Fill();
      } // End of the list of events
    } // End of the list of TTree entries

    // Close the current file
    data_file.Close();
  } // End of the list of files

  // Save the MAUS version, write the tree
  out->cd();
  tree.Write();
  TNamed("MausVersion", maus_version.c_str()).Write();
  out->Close();
}
