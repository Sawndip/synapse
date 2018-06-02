// Cpp includes
#include <iostream>
#include <vector>
#include <map>

// Root includes
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TLegend.h"
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
#include "Pitch.hh"
#include "ProgressBar.hh"
#include "InfoBox.hh"
#include "Globals.hh"
#include "GeometryHandler.hh"
#include "DGaus.hh"
#include "Statistics.hh"
#include "Aperture.hh"
#include "CircleFitter.hh"
#include "ParticleIdentification.hh"

/** @file  ImportMAUSData.cc
 *
 *  @brief Imports the MICE reconstructed data
 *
 *	   Algorithm that imports the MICE data and/or reconstructed simulation along with
 *	   the Monte Carlo truth to a set of 4 simple TNtuples:
 *	    - Data: real data;
 *	    - RecMC: reconstructed Monte Carlo simulation (digitized);
 *	    - Truth: truth that directly corresponds to the selected RecMC;
 *	    - UncutTruth: truth that includes inefficiencies and impurities.
 *	   The code applies a series of cuts on the sample, as defined in the datacards.
 *	   The code produces a plethora of diagnostics, if requested:
 *	    - x, y, px, py, pz, tof profiles (for each cut, each data type);
 *	    - Poincaré sections (2D correlations for each cut, each data type);
 *	    - Counters: records the amount of particles that makes each cut individually.
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
  Globals &globals = Globals::GetInstance(argc, argv);
  if ( !globals.GetDataFiles().size() ) {
    Pitch::print(Pitch::error, "No data file specified, please execute in the following fashion:");
    Pitch::print(Pitch::error, "./import_data [options] data0.root [... dataN.root]");
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

  // Set up and draw the apertures according to the requested model
  Beam::Aperture apertures;
  if ( (double)globals["aperture"] > 0 ) {
    apertures.Add(geoh["tku"].z()-1e3, geoh["tkd"].z()+1e3, globals["aperture"]);
  } else {
    apertures.SetMICEDefault(globals["geometry_filename"]);
  }
  apertures.Draw();

  // If real data is included, set up particle identification to select muons
  ParticleIdentification particle_id;

  // Variables involved in the diagnostics
  std::vector<std::string> vars = {"x", "y", "px", "py", "pz"};		    // Emittance pars tags
  std::vector<std::string> labels = {"x", "y", "p_{x}", "p_{y}", "p_{z}"};  // Emittance pars names
  std::vector<std::string> units = {"mm", "mm", "MeV/c", "MeV/c", "MeV/c"}; // Emittance units
  std::vector<double> llims = {-250, -250, -100, -100, 0};	            // Range lower limits
  std::vector<double> ulims = {250, 250, 100, 100, 350};		    // Range upper limits

  // Diagnostics plots (if requested)
  std::map<std::string, std::string> type_names =
	{{"truth","Truth"}, {"recmc","Simulation"}, {"data","Data"}};
  std::vector<std::string> cut_types = {"none", "pid", "ref", "qual", "fid", "mom", "all"};
  std::map<std::string, std::string> cut_names =
	{{"none","None"}, {"pid", "Particle ID"}, {"ref","Reference"}, {"qual","Quality"},
	 {"fid","Fiducial"}, {"mom","Momentum"}, {"all","All"}};
  std::vector<int> colors = {kBlack, kMagenta+1, kBlue+1, kCyan+1, kGreen+1, kOrange+1, kRed+1};
  std::vector<int> dcolors = {kBlue+2, kRed+2, kGreen+2};
  std::vector<std::string> types;
  std::map<std::string, bool> cuts, dcuts;
  bool tkdevt;
  std::map<std::string, std::map<std::string, std::map<std::string, TH1F*>>> hdiag;
  std::map<std::string, std::map<std::string, THStack*>> hsdiag;
  std::map<std::string, std::map<std::string, size_t>> counters;
  std::map<std::string, std::map<std::string, std::vector<double>>> dmeans;
  std::map<std::string, std::map<std::string, Matrix<double>>> dcovmats;
  std::map<std::string, std::map<std::string, size_t>> dn;
  std::map<std::string, size_t> ids = {{"x",0}, {"y",1}, {"px",2}, {"py",3}};
  TLegend *diag_leg = new TLegend(.75, .6, .89, .89);
  if ( globals["diagnostics"] ) {
    // Set up the legend
    diag_leg->SetLineColorAlpha(0, 0);
    diag_leg->SetFillStyle(0);
    bool leg_set(false);

    // For each category of data, set up histograms and counters
    for (const std::string& type : {"truth", "recmc", "data"})
      if ( globals[type] ) {
        types.push_back(type);

	// Initialize the counters
	counters[type]["evt"] = 0;
	counters[type]["tof0"] = 0;
	counters[type]["tof1"] = 0;
	counters[type]["tkuevt"] = 0;
	for (const std::string& cut : cut_types)
	    counters[type][cut] = 0;
	counters[type]["tkdevt"] = 0;
	for (const std::string& cut : cut_types)
	    counters[type]["tkd"+cut] = 0;

        // Initialize the histogram stacks
        hsdiag[type]["tof"] = new THStack(TString::Format("tof01_%s", type.c_str()), ";t_{01}  [ns]");
	for (const std::string& det : {"tku", "tkd"}) {
          hsdiag[type][det+"ptot"] = new THStack(TString::Format("%s_ptot_%s",
							det.c_str(), type.c_str()), ";p [MeV/c]");
	  for (size_t i = 0; i < vars.size(); i++)
	      hsdiag[type][det+vars[i]] = new THStack(TString::Format("%s_%s_%s", det.c_str(),
		vars[i].c_str(), type.c_str()), TString::Format(";%s [%s]",
		labels[i].c_str(), units[i].c_str()));
	}

        // Intitialize the individual histograms
        size_t id = 0;
        for (const std::string& cut : cut_types) {
	  // Initialize the TOF01 diagnostics histogram
          hdiag[type]["tof"][cut] = new TH1F(TString::Format("tof01_%s_%s", type.c_str(),
		cut.c_str()), TString::Format("%s", cut.c_str()), 100, 24, 40);

	  // Initialize the tracker diagnostics histograms
	  for (const std::string& det : {"tku", "tkd"}) {
            hdiag[type][det+"ptot"][cut] = new TH1F(TString::Format("%s_ptot_%s_%s", det.c_str(),
		type.c_str(), cut.c_str()), TString::Format("%s", cut.c_str()),
		100, llims[4], ulims[4]);
	    for (size_t i = 0; i < vars.size(); i++)
                hdiag[type][det+vars[i]][cut] = new TH1F(TString::Format("%s_%s_%s_%s", det.c_str(),
			vars[i].c_str(), type.c_str(), cut.c_str()), TString::Format("%s",
			cut.c_str()), 100, llims[i], ulims[i]);
	  }

	  // Set the histogram styles
          for (const std::pair<std::string, THStack*>& hist : hsdiag[type]) {
	    hdiag[type][hist.first][cut]->SetLineColor(colors[id]);
	    hdiag[type][hist.first][cut]->SetLineWidth(2);
	    hist.second->Add(hdiag[type][hist.first][cut]);
	  }
          ++id;

	  if ( !leg_set )
		diag_leg->AddEntry(hdiag[type]["tof"][cut], cut_names[cut].c_str(), "l");

	  // Initialize the means and covariance matrices
          dmeans[type][cut] = std::vector<double>(4, 0.);
	  dcovmats[type][cut] = Matrix<double>(4, 4, 0.);
	  dn[type][cut] = 0;
        }

	leg_set = true;
      }
  }

  // Information relevant to the importer
  std::string maus_version = "3.1.2";

  // Loop over the input data files, extract the relevant to the beam evolution (x, y, px, py, pz)
  size_t spilloff = 0;
  std::string out_file = globals["import_filename"];
  TFile *out = new TFile(out_file.c_str(), "RECREATE");

  TNtuple* truth_samples = new TNtuple("truth_ntuple", "",
	"SpillID:EventID:VirtualPlaneID:x:y:z:px:py:pz");
  TNtuple* uncut_truth_samples = new TNtuple("uncut_truth_ntuple", "",
	"SpillID:EventID:VirtualPlaneID:x:y:z:px:py:pz");
  TNtuple* recmc_samples = new TNtuple("recmc_ntuple", "",
	"SpillID:EventID:TrackerID:StationID:x:y:z:px:py:pz:ex:ey:ez:epx:epy:epz");
  TNtuple* data_samples = new TNtuple("data_ntuple", "",
	"SpillID:EventID:TrackerID:StationID:x:y:z:px:py:pz:ex:ey:ez:epx:epy:epz");

  for (const std::string& file : globals.GetDataFiles()) {
    // Set up the ROOT file
    TFile data_file(file.c_str());		// Load the MAUS output file
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
        Pitch::print(Pitch::warning, "MAUS versions don't match, shouldn't proceed");
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
      if ( data_ptr->GetSpill()->GetMCEvents()->size() )
	  simulated = true;
      if ( data_ptr->GetSpill()->GetReconEvents()->size() && (globals["recmc"] || globals["data"]) )
	  recon = true;
    }
    if ( !simulated && !recon ) {
      Pitch::print(Pitch::warning, "Nothing found in this file, skipping");
      continue;
    }

    // If there is reconstructed data in the file, set up the identifier
    if ( recon && globals["do_pid"] ) {
      if ( globals["fit_pid"] ) {
        std::vector<std::string> data_files;
        data_files.push_back(file);
        particle_id = ParticleIdentification(data_files, {0, 1});
      } else {
        particle_id = ParticleIdentification({0, 1}, globals["tofmin"], globals["tofmax"]);
      }
    }

    // Loop over the spills, get the recon events
    double tkfid = (double)globals["tracker_fiducial"];
    ProgressBar pbar;
    for (size_t i = 0; i < (size_t)T->GetEntries(); ++i) {

      // Display the progress in %
      pbar.GetProgress(i, (size_t)T->GetEntries());

      // Update the spill pointed to by data_ptr
      T->GetEntry(i);  
      MAUS::Spill* spill = data_ptr->GetSpill();  
      if (spill == NULL || !(spill->GetDaqEventType() == "physics_event"))
	  continue;

      // Loop over recon events in spill if provided and requested
      // Get the truth stricltly if the recon exists
      if ( recon ) {
        std::vector<MAUS::ReconEvent*>* revts = spill->GetReconEvents();
        std::vector<MAUS::MCEvent*>* mcevts = spill->GetMCEvents();
        for ( size_t ev = 0; ev < revts->size(); ++ev ) {
          if ( !revts->at(ev) )
	      continue;

	  // Reset all the cuts
	  std::string type = simulated ? "recmc" : "data";
	  for (const std::string cut : cut_types) {
	    cuts[cut] = true;
	    dcuts[cut] = true;
	  }

	  // If diagnostics are requested, get the TOF01 time-of-flight
	  double t01(0);
	  if ( globals["diagnostics"] && !globals["notofs"] ) {

	    // Get the TOF space points
	    MAUS::TOFEvent* tofevent = revts->at(ev)->GetTOFEvent();
            MAUS::TOFEventSpacePoint tofsp = tofevent->GetTOFEventSpacePoint();
            std::vector<std::vector<MAUS::TOFSpacePoint>> tofsps = {tofsp.GetTOF0SpacePointArray(),
							   	    tofsp.GetTOF1SpacePointArray()};

  	    // Get the times and positions in each TOF, return if not a single SP
  	    std::vector<double> t(2);
	    double both(true);
  	    size_t j;
	    counters[type]["evt"]++;
  	    for (j = 0 ; j < 2; j++) {
              if ( tofsps[j].size() == 1 ) {
	        counters[type]["tof"+std::to_string((int)j)]++;
                t[j] = tofsps[j][0].GetTime();
              } else {
                both = false;
              }
            }
	    if ( !both )
		continue;

	    t01 = t[1]-t[0];
	  }

	  // If PID was requested, reject the pions and electrons from the sample, keep only muons
	  if ( globals["do_pid"] && !globals["notofs"] ) {
	    int pid = particle_id.GetID(revts->at(ev));
	    if ( !pid ) {
	      continue;
	    } else if ( pid != 13 ) {
	      cuts["pid"] = false;
	    }
	  } else {
	    cuts["pid"] = true;
	  }

          //////////////////////////////
          ////// TKU&D COORDINATES /////
          //////////////////////////////
	  // Fetch the kalman tracks and extract the phase space coordinates
          MAUS::SciFiEvent* sfevt = revts->at(ev)->GetSciFiEvent(); // Pull out scifi event
          MAUS::SciFiTrackPArray sfts = sfevt->scifitracks();

	  // Track point and error container
          std::map<std::string, bool> has_track = {{"tku",false}, {"tkd",false}};
          std::map<std::string, std::vector<std::map<std::string, double>>> vpoint, verror, tvpoint;
	  for ( const std::string &det : dets ) {
	    vpoint[det].resize(5);
	    verror[det].resize(5);
	    tvpoint[det].resize(5);
	  }

	  // Loop over Kalman tracks (helical)
 	  tkdevt = false;
          for ( const MAUS::SciFiTrack *track : sfts ) {
	    // If it is not fitted with a helix, no momentum info, not a valid track
	    std::string det = track->tracker() ? "tkd" : "tku";
	    if ( det == "tku" ) {
	      counters[type]["tkuevt"]++;
	    } else {
	      tkdevt = true;
	    }
	    if ( !track->pr_track_pointer_helical() || !track->GetAlgorithmUsed() )
	        continue;

	    // Only keep the tracks that have a decent chi^2/ndf
	    if ( track->chi2()/track->ndf() > 10) {
	      if ( det == "tku" ) {
		cuts["qual"] = false;
	      } else {
		dcuts["qual"] = false;
	      }
	    }

	    // Get the underlying helix from PR, check that it stays with within fiducial
	    MAUS::SciFiHelicalPRTrack* prtrack = track->pr_track_pointer_helical();
	    double x0 = prtrack->get_circle_x0();
	    double y0 = prtrack->get_circle_y0();
	    double rad = fabs(prtrack->get_R());
	    double maxrad = sqrt(x0*x0+y0*y0)+rad;
	    if ( maxrad > tkfid ) {
	      has_track[det] = false;
	      if ( det == "tku" ) {
		cuts["fid"] = false;
	      } else {
		dcuts["fid"] = false;
	      }
	    }

	    // Loop over the Kalman track points
	    has_track[det] = true;
	    MAUS::SciFiTrackPointPArray tps = track->scifitrackpoints();
            MAUS::ThreeVector pos, pos_error, mom, mom_error;
	    for (const MAUS::SciFiTrackPoint *tpoint : tps ) {

	      if ( !tpoint->plane() ) {
	      	pos = tpoint->pos();
	      	pos_error = tpoint->pos_error();
	      	mom = tpoint->mom();
	      	mom_error = tpoint->mom_error();
	      	size_t st = tpoint->station()-1;

		// Fill the containers
	        vpoint[det][st]["x"] = pos.x();
	        vpoint[det][st]["y"] = pos.y();
	        vpoint[det][st]["z"] = pos.z();
	        vpoint[det][st]["px"] = mom.x();
	        vpoint[det][st]["py"] = mom.y();
	        vpoint[det][st]["pz"] = mom.z();

	        verror[det][st]["x"] = pos_error.x();
	        verror[det][st]["y"] = pos_error.y();
	        verror[det][st]["z"] = pos_error.z();
	        verror[det][st]["px"] = mom_error.x();
	        verror[det][st]["py"] = mom_error.y();
	        verror[det][st]["pz"] = mom_error.z();

	        // Check that there are no NaN values, should be removed when moving to 2.9.0 (TODO)
	        for ( const std::string &var : {"x", "y", "z", "px", "py", "pz"} ) {
		  if ( std::isnan(vpoint[det][st][var]) || std::isnan(verror[det][st][var]) ) {
		    has_track[det] = false;
		    break;
		  } 
		}
		if ( !has_track[det] )
		    break;
	      }
	    }
	  }

	  // If a tracker was not hit, set the cuts to false
	  if ( !has_track["tku"] ) {
	    cuts["ref"] = false;
	    cuts["fid"] = false;
	    cuts["qual"] = false;
          }
	  if ( !has_track["tkd"] ) {
	    dcuts["ref"] = false;
	    dcuts["fid"] = false;
	    dcuts["qual"] = false;
          }

	  // Check that both trackers have a full track if through particles are requested
	  if ( globals["through_part"] && !has_track["tkd"] )
	      continue;

	  // Only keep the momenta selected at the reference plane
	  cuts["mom"] = false;
	  if ( globals["total_momentum_cut"] ) {
	    MAUS::ThreeVector upv
		(vpoint["tku"][0]["px"], vpoint["tku"][0]["py"], vpoint["tku"][0]["pz"]);
	    if ( upv.mag() > (double)globals["min_momentum"] &&
	         upv.mag() < (double)globals["max_momentum"] )
	        cuts["mom"] = true;
          } else {
	    if ( vpoint["tku"][0]["pz"] > (double)globals["min_momentum"] &&
	         vpoint["tku"][0]["pz"] < (double)globals["max_momentum"] )
	        cuts["mom"] = true;
	  }

	  // Check that the downstream track does not have an absurd momentum
	  dcuts["mom"] = false;
	  if ( has_track["tkd"] ) {
	    MAUS::ThreeVector downv
		(vpoint["tkd"][0]["px"], vpoint["tkd"][0]["py"], vpoint["tkd"][0]["pz"]);
	    if ( downv.mag() > (double)globals["tkd_min_momentum"] &&
	         downv.mag() < (double)globals["tkd_max_momentum"] ) {
	      dcuts["mom"] = true;
	    }
          }

	  // For now, test all the cuts and only keep events that pass them all
	  for (const std::string& cut : cut_types)
	    if ( !cuts[cut] ) {
	      cuts["all"] = false;
	      break;
	    }
	  for (const std::string& cut : cut_types)
	    if ( !dcuts[cut] ) {
	      dcuts["all"] = false;
	      break;
	    }

	  // Fill the diagnostics graphs
	  if ( globals["diagnostics"] ) {
	    for (const std::string& cut : cut_types)
	      if ( cuts["all"] && dcuts[cut] )
		  counters[type]["tkd"+cut]++;
	    for (const std::string& cut : cut_types) {
	      if ( cuts[cut] ) {
		counters[type][cut]++;
		if ( !globals["notofs"] )
		    hdiag[type]["tof"][cut]->Fill(t01);
		for (const std::string& det : {"tku", "tkd"}) {
		  if ( has_track[det] ) {
		    if ( det == "tkd" && !dcuts["all"] )
			continue;
		    double ptot = 0.;
		    for (const std::string& var : {"px", "py", "pz"})
			ptot += vpoint[det][0][var]*vpoint[det][0][var];
		    hdiag[type][det+"ptot"][cut]->Fill(sqrt(ptot));

		    for (const std::string& var : vars)
		        hdiag[type][det+var][cut]->Fill(vpoint[det][0][var]);
		    if ( det == "tku" ) {
			std::vector<double> vec(4);
			for (size_t k = 0; k < 4; k++)
			    vec[k] = vpoint[det][0][vars[k]];
			Math::IncrementCovarianceMatrix(dn[type][cut],
				dcovmats[type][cut], dmeans[type][cut], vec);
			for (size_t k = 0; k < 4; k++)
			    Math::IncrementMean(dn[type][cut], dmeans[type][cut][k], vec[k]);
			++dn[type][cut];
		    }
		  }
		}
	      }
	    }
	  }


	  if ( !cuts["all"] ) {
	    continue;
	  } else if ( tkdevt ) {
	    counters[type]["tkdevt"]++;
          }
	  if ( !dcuts["all"] )
	      has_track["tkd"] = false;

	  // If simulated, get the truth that corresponds to the reconstruction
	  // If it cannot be found do not save the reconstruction either
	  std::map<std::string, int> vids =
		  {{"tku",-1},{"tkd",-1},{"tkd3",-1},{"tof0",-1},{"tof1",-1},{"min",-1},{"max",-1}};
	  if ( simulated ) {

            std::vector<MAUS::VirtualHit>* vhits = mcevts->at(ev)->GetVirtualHits();
  	    MAUS::VirtualHit vhit;
	    MAUS::ThreeVector pos, mom;
	    size_t vid;
            for ( size_t vh = 0; vh < std::min((size_t)globals["nvirtuals"], vhits->size()); vh++ ) {
	      vhit = vhits->at(vh);
	      vid = vhit.GetStationId()-1;
	      pos = vhit.GetPosition();

	      for (const std::string& det : {"tku", "tkd", "tkd3", "tof0", "tof1", "min", "max"})
	        if ( vid == (size_t)globals[det+"_vid"] )
	            vids[det] = vh;
	    }

	    if ( vids["tof1"] > 0 || globals["notofs"] ) {
	      counters["truth"]["evt"]++;
	      counters["truth"]["tof1"]++;
	      counters["truth"]["tof0"]++;
	      if ( !globals["notofs"] )
	          t01 = vhits->at(vids["tof1"]).GetTime() - vhits->at(vids["tof0"]).GetTime();
            }

	    if ( globals["diagnostics"] ) {
	      for (const std::string& cut : cut_types) {
	        if ( cuts[cut] ) {
		  counters["truth"][cut]++;
		  if ( vids["tof1"] > 0 && !globals["notofs"] )
		      hdiag["truth"]["tof"][cut]->Fill(t01);

		  for (const std::string& det : {"tku", "tkd"}) {
		    if ( has_track[det] && vids[det] > 0 ) {
		      vhit = vhits->at(vids[det]);
	      	      pos = vhit.GetPosition();
	      	      mom = vhit.GetMomentum();
		      tvpoint[det][0]["x"] = pos.x();
		      tvpoint[det][0]["y"] = pos.y();
		      tvpoint[det][0]["px"] = mom.x();
		      tvpoint[det][0]["py"] = mom.y();
		      tvpoint[det][0]["pz"] = mom.z();

		      double ptot = 0.;
		      for (const std::string& var : {"px", "py", "pz"})
			  ptot += tvpoint[det][0][var]*tvpoint[det][0][var];
		      hdiag["truth"][det+"ptot"][cut]->Fill(sqrt(ptot));

		      for (const std::string& var : vars)
		          hdiag["truth"][det+var][cut]->Fill(tvpoint[det][0][var]);
		      if ( det == "tku" ) {
			  std::vector<double> vec(4);
			  for (size_t k = 0; k < 4; k++)
			      vec[k] = tvpoint[det][0][vars[k]];
			  Math::IncrementCovarianceMatrix(dn["truth"][cut],
				dcovmats["truth"][cut], dmeans["truth"][cut], vec);
			  for (size_t k = 0; k < 4; k++)
			      Math::IncrementMean(dn["truth"][cut], dmeans["truth"][cut][k], vec[k]);
			  ++dn["truth"][cut];
		      }
		    }
		  }
	        }
	      }
	    }

            // Add the truth for each particle that makes it into the reconstructed sample
            for ( size_t vh = 0; vh < std::min((size_t)globals["nvirtuals"], vhits->size()); vh++ ) {
  	      vhit = vhits->at(vh);
	      vid = vhit.GetStationId()-1;
	      pos = vhit.GetPosition();
	      mom = vhit.GetMomentum();

	      // In between the two trackers, apply the default apertures
	      if ( vid > (size_t)globals["tku_vid"] && vid < (size_t)globals["tkd_vid"]
		   && !apertures.IsIn(pos.x(), pos.y(), pos.z()) )
		  continue;

	      // When the downtream tracker is reached, use it as a criterion
	      if ( vid >= (size_t)globals["tkd_vid"] && !has_track["tkd"]  )
	          break; 

	      truth_samples->Fill(
		spilloff+i, ev, vid, pos.x(), pos.y(), pos.z(), mom.x(), mom.y(), mom.z());
            }

            // Add the uncut truth for each upstream particle that makes the cut, regardless of
	    // whether the downstream tracker is present. First check the true fiducial
	    std::map<std::string, bool> fiducial;
	    fiducial["tku"] = false;
	    if ( vids["tku"] > 0 ) {
	      fiducial["tku"] = true;
	      std::vector<double> xx, yy;
              for ( size_t vh = vids["min"]; vh <= (size_t)vids["tku"]; vh++ ) {
  	        vhit = vhits->at(vh);
		pos = vhit.GetPosition();
		xx.push_back(pos.x());
		yy.push_back(pos.y());
	      }
	      double tx0, ty0, trad;
	      FitCircle(xx, yy, tx0, ty0, trad);
	      double maxrad = sqrt(tx0*tx0+ty0*ty0)+trad;
	      if ( maxrad > tkfid )
		  fiducial["tku"] = false;
	    }
	    fiducial["tkd"] = false;
	    if ( vids["tkd3"] > 0 ) { // Must at least make it to station 3 to form a track
	      fiducial["tkd"] = true;
	      size_t limit = vids["max"] > 0 ? vids["max"] : vhits->size();
	      std::vector<double> xx, yy;
              for ( size_t vh = vids["tkd"]; vh < limit; vh++ ) {
  	        vhit = vhits->at(vh);
		pos = vhit.GetPosition();
		xx.push_back(pos.x());
		yy.push_back(pos.y());
	      }
	      double tx0, ty0, trad;
	      FitCircle(xx, yy, tx0, ty0, trad);
	      double maxrad = sqrt(tx0*tx0+ty0*ty0)+trad;
	      if ( maxrad > tkfid )
		  fiducial["tkd"] = false;
	    }

            for ( size_t vh = 0; vh < std::min((size_t)globals["nvirtuals"], vhits->size()); vh++ ) {
  	      vhit = vhits->at(vh);
	      vid = vhit.GetStationId()-1;
	      pos = vhit.GetPosition();
	      mom = vhit.GetMomentum();

	      // Do not bother if TKU ref. plane is not hit (should be impossible...)
	      if ( vids["tku"] < 0 )
		  break;

	      // In the upstream tracker, check that the particle is indeed a muon and
	      // that it indeed lays within the fiducial volume of the upstream tracker
	      if ( vid <= (size_t)globals["tku_vid"] ) {
		if ( vhits->at(vids["tku"]).GetParticleId() != -13 )
		    break; // Will not turn back into a muon, can interupt
		if ( !fiducial["tku"] )
		    continue;
	      }

	      // In between the two trackers, apply the default apertures
	      if ( vid > (size_t)globals["tku_vid"] && vid < (size_t)globals["tkd_vid"] )
		if ( !apertures.IsIn(pos.x(), pos.y(), pos.z()) )
		    continue;

	      // When the downtream tracker is reached, check that the track is within the fiducial
	      // throughout the tracker, and that it is a muon. Do not request a TKD recon track
	      if ( vid >= (size_t)globals["tkd_vid"] ) {
		if ( vids["tkd"] < 0 ) // Downstream ref. plane missing
		    break;
		if ( vhits->at(vids["tkd"]).GetParticleId() != -13 )
		    break; // Will not turn back into a muon, can interupt
		if ( !fiducial["tkd"] )
		    continue;
	      }

	      uncut_truth_samples->Fill(
		spilloff+i, ev, vid, pos.x(), pos.y(), pos.z(), mom.x(), mom.y(), mom.z());
            }
	  }

	  // Fill ther reconstructed samples
	  std::vector<float> array = {spilloff+(float)i, (float)ev};
	  for ( const std::string &det : dets ) {
	    if ( !has_track[det] || ( simulated && vids[det] < 0 ) )
		continue;

	    for ( size_t st = 0; st < 5; st++ ) {
	      array.resize(2);
	      array.push_back((det == "tku") ? 0 : 1);
	      array.push_back(st);

	      for ( const std::string &var : {"x", "y", "z", "px", "py", "pz"} )
	          array.push_back(vpoint[det][st][var]);
	      for ( const std::string &var : {"x", "y", "z", "px", "py", "pz"} )
	          array.push_back(verror[det][st][var]);

	      if ( simulated ) {
	          recmc_samples->Fill(&(array[0]));
	      } else {
	          data_samples->Fill(&(array[0]));
              }

	    }
	  }

        } // End of the list of recEvents
      }

    } // End of the list of TTree entries
    data_file.Close();
  } // End of the list of files

  // Fill the output file with the Ntuples and information about the sample
  out->cd();
  if ( truth_samples->GetEntries() )
      truth_samples->Write("Truth");
  delete truth_samples;

  if ( uncut_truth_samples->GetEntries() )
      uncut_truth_samples->Write("UncutTruth");
  delete uncut_truth_samples;

  if ( recmc_samples->GetEntries() )
      recmc_samples->Write("RecMC");
  delete recmc_samples;

  if ( data_samples->GetEntries() )
      data_samples->Write("Data");
  delete data_samples;

  TNamed("MausVersion", maus_version.c_str()).Write();
  out->Close();

  // If the diagnostics plots are requested, print them
  if ( globals["diagnostics"] ) {

    // Print the counters to a data file
    ofstream fcounters;
    fcounters.open("diag_counters.txt");
    for (const std::string type : types) {
      fcounters << type_names[type] << "\n";
      fcounters << "Event\t" << counters[type]["evt"] << "\n";
      fcounters << "TOF0 SP\t" << counters[type]["tof0"] << "\n";
      fcounters << "TOF1 SP\t" << counters[type]["tof1"] << "\n";
      fcounters << "Both TOFs\t" << counters[type]["none"] << "\n";
      fcounters << "Muon\t" << counters[type]["pid"] << "\n";
      fcounters << "\n";
      fcounters << "TKU Event\t" << counters[type]["tkuevt"] << "\n";
      fcounters << "TKU Reference\t" << counters[type]["ref"] << "\n";
      fcounters << "TKU Quality\t" << counters[type]["qual"] << "\n";
      fcounters << "TKU Fiducial\t" << counters[type]["fid"] << "\n";
      fcounters << "TKU Momentum\t" << counters[type]["mom"] << "\n";
      fcounters << "TKU All\t" << counters[type]["all"] << "\n";
      fcounters << "\n";
      fcounters << "TKD Event\t" << counters[type]["tkdevt"] << "\n";
      fcounters << "TKD Reference\t" << counters[type]["tkdref"] << "\n";
      fcounters << "TKD Quality\t" << counters[type]["tkdqual"] << "\n";
      fcounters << "TKD Fiducial\t" << counters[type]["tkdfid"] << "\n";
      fcounters << "TKD Momentum\t" << counters[type]["tkdmom"] << "\n";
      fcounters << "TKD All\t" << counters[type]["tkdall"] << "\n";
      fcounters << "\n";
    }
    fcounters.close();

    // Initialize the RMS ellipses from the recorded means and covariances (+ frames)
    TFile outfile("diag_histograms.root", "RECREATE");
    std::map<std::string, std::map<std::string, std::map<std::string, TEllipse*>>> ellipses;
    std::map<std::string, std::map<std::string, std::map<std::string, TGraph*>>> centres;
    std::map<std::string, std::map<std::string, TH2F*>> frames;
    std::map<std::string, std::vector<double>> mean;
    std::map<std::string, Matrix<double>> mat;
    DGaus gaus;
    double xmin, xmax, ymin, ymax;
    double margin = .5;
    size_t i, j, id;
    for (const std::string type : types)
      for (i = 0; i < 4; i++)
	for (j = i+1; j < 4; j++) {
	  id = 0;
          for (const std::string cut : cut_types) {
	    mean[cut] = {dmeans[type][cut][i], dmeans[type][cut][j]};
	    mat[cut] = Matrix<double>({{dcovmats[type][cut][i][i], dcovmats[type][cut][i][j]},
		     		      {dcovmats[type][cut][j][i], dcovmats[type][cut][j][j]}});
	    gaus = DGaus(mean[cut], mat[cut]);
	    ellipses[type][cut][vars[i]+vars[j]] = (TEllipse*)gaus.Contour2D(.39347)[0];
	    ellipses[type][cut][vars[i]+vars[j]]->SetLineColor(colors[id]);

	    centres[type][cut][vars[i]+vars[j]] = new TGraph();
	    centres[type][cut][vars[i]+vars[j]]->SetPoint(0, mean[cut][0], mean[cut][1]);
	    centres[type][cut][vars[i]+vars[j]]->SetMarkerStyle(kFullCross);
	    centres[type][cut][vars[i]+vars[j]]->SetMarkerColor(colors[id]);
	    ++id;
	  }

	  xmin = mean["none"][0]-(1.+margin)*sqrt(mat["none"][0][0]);
	  xmax = mean["none"][0]+(1.+margin)*sqrt(mat["none"][0][0]);
	  ymin = mean["none"][1]-(1.+margin)*sqrt(mat["none"][1][1]);
	  ymax = mean["none"][1]+(1.+margin)*sqrt(mat["none"][1][1]);
	  frames[type][vars[i]+vars[j]] = new TH2F(TString::Format("frame_%s_%s", 
		type.c_str(), (vars[i]+vars[j]).c_str()), 
		TString::Format(";%s [%s];%s [%s]", labels[i].c_str(), units[i].c_str(),
		labels[j].c_str(), units[j].c_str()), 1, xmin, xmax, 1, ymin, ymax);
	}

    // Print the comparision between the different levels of cut
    for (const std::string& type : types) {

      // Initialize the info box
      std::string data_type = (type == "data") ? "Preliminary" : "[simulation]";
      InfoBox info(data_type, maus_version,
  	  globals["run_name"].AsString(), globals["user_cycle"].AsString());
      info.SetPosition("tl");

      // Draw the stacks
      for (const std::pair<std::string, THStack*>& hist : hsdiag[type]) {
	TCanvas* c = new TCanvas("c", "c", 1200, 800);
	hist.second->Draw("NOSTACK");
	hist.second->Write(hist.second->GetName());
        gPad->SetLogy();
	diag_leg->Draw("SAME");
	info.Draw();
	c->SaveAs(TString::Format("diag_%s_%s.pdf", hist.first.c_str(), type.c_str()));
	delete c;
      }

      // Draw the RMS ellipses
      gStyle->SetOptStat(0);
      for (const std::pair<std::string, TEllipse*>& ell : ellipses[type]["none"]) {
        TCanvas* c = new TCanvas("c", "c", 1200, 800);
        frames[type][ell.first]->Draw();
        frames[type][ell.first]->Write(TString::Format("frame_%s_%s",
						type.c_str(), ell.first.c_str()));
      
        for (const std::string cut : cut_types) {
	  centres[type][cut][ell.first]->Draw("PSAME");
	  centres[type][cut][ell.first]->Write(TString::Format("centroid_%s_%s_%s",
						type.c_str(), cut.c_str(), ell.first.c_str()));
	  ellipses[type][cut][ell.first]->Draw("SAME");
	  ellipses[type][cut][ell.first]->Write(TString::Format("ellipse_%s_%s_%s",
						type.c_str(), cut.c_str(), ell.first.c_str()));
	}
        diag_leg->Draw("SAME");
	info.Draw();
        c->SaveAs(TString::Format("diag_%s_%s.pdf", ell.first.c_str(), type.c_str()));
        delete c;
      }
    }

    // Initialize the info box
    std::string data_type = globals["data"] ? "Preliminary" : "[simulation]";
    InfoBox info(data_type, maus_version,
  	globals["run_name"].AsString(), globals["user_cycle"].AsString());
    info.SetPosition("tl");

    // Print the comparision between the different types of data histograms
    std::string cut = globals["compare_cut"];
    THStack* hscomp;
    TLegend *comp_leg = new TLegend(.7, .75, .89, .89);
    comp_leg->SetFillStyle(0);
    comp_leg->SetLineColorAlpha(0, 0);
    bool leg_set(false);

    for (const std::pair<std::string, THStack*>& hist : hsdiag[types[0]]) {
      std::string title = std::string(hsdiag[types[0]][hist.first]->GetTitle());
      hscomp = new THStack(TString::Format("%s_comp", hist.first.c_str()), title.c_str());

      id = 0;
      size_t ref_n = hdiag[types[0]][hist.first][cut]->GetEntries();
      for (const std::string& type : types) {
	hdiag[type][hist.first][cut]->SetLineColor(dcolors[id]);
	hdiag[type][hist.first][cut]->Scale((double)ref_n/hdiag[type][hist.first][cut]->GetEntries());
	hscomp->Add(hdiag[type][hist.first][cut]);
	if ( !leg_set )
	    comp_leg->AddEntry(hdiag[type][hist.first][cut], type_names[type].c_str(), "l");
	++id;
      }
      leg_set = true;

      TCanvas* c = new TCanvas("c", "c", 1200, 800);
      gPad->SetLogy();
      hscomp->Draw("NOSTACK");
      hscomp->Write(hscomp->GetName());
      comp_leg->Draw("SAME");
      info.Draw();
      c->SaveAs(TString::Format("diag_%s_comp.pdf", hist.first.c_str()));
      delete c;
    }

    // Print the comparision between the different types of data ellipses
    for (const std::pair<std::string, TEllipse*>& ell : ellipses[types[0]][cut]) {
      TCanvas* c = new TCanvas("c", "c", 1200, 800);
      frames[types[0]][ell.first]->Draw();
      id = 0;
      for (const std::string& type : types) {
	centres[type][cut][ell.first]->SetMarkerColor(dcolors[id]);
	centres[type][cut][ell.first]->Draw("PSAME");
	ellipses[type][cut][ell.first]->SetLineColor(dcolors[id]);
	ellipses[type][cut][ell.first]->Draw("SAME");
	++id;
      }
      comp_leg->Draw("SAME");
      info.Draw();
      c->SaveAs(TString::Format("diag_%s_comp.pdf", ell.first.c_str()));
      delete c;
    }

    // Move all the output to an appropriate directory
    outfile.Close();
    std::string uc_name = globals["user_cycle"];
    uc_name.replace(uc_name.find("/"), 1, "-");
    std::string run_name = globals["run_name"];
    if ( run_name.find(".") != std::string::npos )
        run_name.replace(run_name.find("."), 1, "-");
    std::string dir_name = uc_name+"_"+run_name;

    Pitch::print(Pitch::info, "Moving the diagnostics graphs to "+dir_name+"/diag");
    std::string sysCmd = "mkdir -p "+dir_name+"/diag; mv diag_* "+dir_name+"/diag";
    if ( std::system(sysCmd.c_str()) )
        Pitch::print(Pitch::error, "Couldn't move diagnostics plots");
  }
}
