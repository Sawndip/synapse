// Cpp includes
#include <iostream>
#include <vector>
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
#include "Pitch.hh"
#include "ProgressBar.hh"
#include "InfoBox.hh"
#include "Globals.hh"
#include "GeometryHandler.hh"
#include "DGaus.hh"
#include "Statistics.hh"
#include "Aperture.hh"
#include "CircleFitter.hh"

/** @file  ImportMAUSSimulation.cc
 *
 *  @brief Imports the MAUS simulation truth standalone
 *
 *         Algorithm that imports the MAUS simulation truth standalone, without the
 * 	   corresponding digitized information. It stores it in a single TNtuple:
 *	    - Truth: true phase space information of each track at each virtual plane.
 *	   The code applies a series of cuts on the sample, as defined in the datacards.
 *	   The code produces a plethora of diagnostics, if requested:
 *	    - x, y, px, py, pz, tof profiles (for each cut, each data type);
 *	    - Poincaré sections (2D correlations for each cut, each data type);
 *	    - Counters: records the amount of particles that makes each cut individually.
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
    Pitch::print(Pitch::error, "./import_sim [options] dat0.root [... dataN.root]");
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

  // Variables involved in the diagnostics
  std::vector<std::string> vars = {"x", "y", "px", "py", "pz"};		    // Emittance pars tags
  std::vector<std::string> labels = {"x", "y", "p_{x}", "p_{y}", "p_{z}"};  // Emittance pars names
  std::vector<std::string> units = {"mm", "mm", "MeV/c", "MeV/c", "MeV/c"}; // Emittance units
  std::vector<double> llims = {-250, -250, -100, -100, 0};	            // Range lower limits
  std::vector<double> ulims = {250, 250, 100, 100, 350};		    // Range upper limits

  // Diagnostics plots (if requested)
  std::vector<std::string> cut_types = {"none", "pid", "ref", "qual", "fid", "mom", "all"};
  std::map<std::string, std::string> cut_names =
	{{"none","None"}, {"pid", "Particle ID"}, {"ref","Reference"}, {"qual","Quality"},
	 {"fid","Fiducial"}, {"mom","Momentum"}, {"all","All"}};
  std::vector<int> colors = {kBlack, kMagenta+1, kBlue+1, kCyan+1, kGreen+1, kOrange+1, kRed+1};
  std::vector<int> dcolors = {kBlue+2, kRed+2, kGreen+2};
  std::map<std::string, bool> cuts;
  std::map<std::string, std::map<std::string, TH1F*>> hdiag;
  std::map<std::string, THStack*> hsdiag;
  std::map<std::string, size_t> counters;
  std::map<std::string, std::vector<double>> dmeans;
  std::map<std::string, Matrix<double>> dcovmats;
  std::map<std::string, size_t> dn;
  std::map<std::string, size_t> ids = {{"x",0}, {"y",1}, {"px",2}, {"py",3}};
  TLegend *diag_leg = new TLegend(.75, .6, .89, .89);
  if ( globals["diagnostics"] ) {
    // Set up the legend
    diag_leg->SetLineColorAlpha(0, 0);
    diag_leg->SetFillStyle(0);

    // Initialize the counters
    counters["evt"] = 0;
    counters["tof0"] = 0;
    counters["tof1"] = 0;
    counters["tkuevt"] = 0;
    for (const std::string& cut : cut_types)
	counters[cut] = 0;
    counters["tkdevt"] = 0;
    for (const std::string& cut : cut_types)
	counters["tkd"+cut] = 0;

    // Initialize the histogram stacks
    hsdiag["tof"] = new THStack("tof01_truth", ";t_{01}  [ns]");
    for (const std::string& det : {"tku", "tkd"}) {
      hsdiag[det+"ptot"] = new THStack(TString::Format("%s_ptot_truth", det.c_str()), ";p [MeV/c]");
      for (size_t i = 0; i < vars.size(); i++)
	  hsdiag[det+vars[i]] = new THStack(TString::Format("%s_%s_truth", det.c_str(),
		vars[i].c_str()), TString::Format(";%s [%s]",
		labels[i].c_str(), units[i].c_str()));
    }

    // Intitialize the individual histograms
    size_t id = 0;
    for (const std::string& cut : cut_types) {
      // Initialize the TOF01 diagnostics histogram
      hdiag["tof"][cut] = new TH1F(TString::Format("tof01_truth_%s",
	cut.c_str()), TString::Format("%s", cut.c_str()), 100, 24, 40);

      // Initialize the tracker diagnostics histograms
      for (const std::string& det : {"tku", "tkd"}) {
        hdiag[det+"ptot"][cut] = new TH1F(TString::Format("%s_ptot_truth_%s", det.c_str(),
		cut.c_str()), TString::Format("%s", cut.c_str()),
		100, llims[4], ulims[4]);
      	for (size_t i = 0; i < vars.size(); i++)
	    hdiag[det+vars[i]][cut] = new TH1F(TString::Format("%s_%s_truth_%s", det.c_str(),
		vars[i].c_str(), cut.c_str()), TString::Format("%s",
		cut.c_str()), 100, llims[i], ulims[i]);
      }

      // Set the histogram styles
      for (const std::pair<std::string, THStack*>& hist : hsdiag) {
	hdiag[hist.first][cut]->SetLineColor(colors[id]);
	hdiag[hist.first][cut]->SetLineWidth(2);
	hist.second->Add(hdiag[hist.first][cut]);
      }
      ++id;

      // Add an entry to the TLegend
      diag_leg->AddEntry(hdiag["tof"][cut], cut_names[cut].c_str(), "l");

      // Initialize the means and covariance matrices
      dmeans[cut] = std::vector<double>(4, 0.);
      dcovmats[cut] = Matrix<double>(4, 4, 0.);
      dn[cut] = 0;
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
    bool simulated(false);
    for (size_t i = 0; i < std::min((size_t)20, (size_t)T->GetEntries()); i++) {
      T->GetEntry(i);
      if ( data_ptr->GetSpill()->GetMCEvents()->size() )
	  simulated = true;
    }
    if ( !simulated ) {
      Pitch::print(Pitch::warning, "No MC truth found in this file, skipping");
      continue;
    }

    // Loop over the spills, get the recon events
    ProgressBar pbar;
    for (size_t i = 0; i < (size_t)T->GetEntries(); ++i) {

      // Display the progress in %
      pbar.GetProgress(i, (size_t)T->GetEntries());

      // Update the spill pointed to by data_ptr
      T->GetEntry(i);  
      MAUS::Spill* spill = data_ptr->GetSpill();  
      if (spill == NULL || !(spill->GetDaqEventType() == "physics_event"))
	  continue;

      // Loop over MC events in spill if provided
      if ( simulated ) { 
        std::vector<MAUS::MCEvent*>* mcevts = spill->GetMCEvents();
        for ( size_t ev = 0; ev < mcevts->size(); ++ev ) {
          if ( !mcevts->at(ev) )
	      continue;

	  // Reset all the cuts
	  for (const std::string cut : cut_types)
	      cuts[cut] = true;

          //////////////////////////////
          ////// TKU&D COORDINATES /////
          //////////////////////////////
          // Loop over the virtual hits, identify the virtual plane ids
          std::vector<MAUS::VirtualHit>* vhits = mcevts->at(ev)->GetVirtualHits();
	  bool lost(false);
	  MAUS::VirtualHit vhit;
	  MAUS::ThreeVector pos;
	  size_t vid;
	  std::map<std::string, int> vids =
		{{"tku",-1}, {"tkd",-1}, {"tof0",-1}, {"tof1",-1}, {"max", -1}};
          for ( size_t vh = 0; vh < std::min((size_t)globals["nvirtuals"], vhits->size()); vh++ ) {
	    vhit = vhits->at(vh);
	    vid = vhit.GetStationId()-1;
	    pos = vhit.GetPosition();

	    if ( !apertures.IsIn(pos.x(), pos.y(), pos.z()) ) {
	      lost = true;
	      break;
	    }

	    for (const std::string& det : {"tku", "tkd", "tof0", "tof1", "max"})
	      if ( vid == (size_t)globals[det+"_vid"] )
	          vids[det] = vh;

	    if ( vid == (size_t)globals["tof1_vid"] ) {
	      counters["evt"]++;
	      counters["tof1"]++;
	      counters["tof0"]++;
            }
	  }
	  if ( globals["through_part"] && (lost || vhits->size() < (size_t)globals["nvirtuals"]) )
	      continue;

	  // Get the time-of-flight, reject tracks that didn't get to TOF1
	  if ( vids["tof0"] < 0 || vids["tof1"] < 0 )
	      continue;
	  double t01(0);
	  if ( globals["diagnostics"] )
	      t01 = vhits->at(vids["tof1"]).GetTime()-vhits->at(vids["tof0"]).GetTime();

	  // Check that the particle is a muon at TOF1
	  if ( fabs(vhits->at(vids["tof1"]).GetParticleId()) != 13 )
	      cuts["pid"] = false;

  	  // Check the relevant selection criteria
	  if ( vids["tku"] < 0 ) {
	    for (const std::string& cut : {"ref", "qual", "fid", "mom"})
	        cuts[cut] = false;
	  } else {
	    vhit = vhits->at(vids["tku"]);
	    if ( globals["total_momentum_cut"] ) {
	      if ( vhit.GetMomentum().mag() < (double)globals["min_momentum"] ||
		   vhit.GetMomentum().mag() > (double)globals["max_momentum"] )
		  cuts["mom"] = false;
	    } else {
	      if ( vhit.GetMomentum().z() < (double)globals["min_momentum"] ||
		   vhit.GetMomentum().z() > (double)globals["max_momentum"] )
		  cuts["mom"] = false;
	    }
	  }

	  // For now, test all the cuts and only keep events that pass them all
	  for (const std::string& cut : cut_types)
	    if ( !cuts[cut] ) {
	      cuts["all"] = false;
	      break;
	    }

	  if ( globals["diagnostics"] ) {
	    MAUS::VirtualHit vhit;
	    MAUS::ThreeVector pos, mom;
	    std::map<std::string, std::map<std::string, double>> tpoint;
	    for (const std::string& det : {"tku", "tkd"}) {
	      if ( vids[det] < 0 )
		  continue;

	      vhit = vhits->at(vids[det]);
	      pos = vhit.GetPosition();
	      mom = vhit.GetMomentum();
	      tpoint[det] = {{"x", pos.x()}, {"y", pos.y()},
		{"px", mom.x()}, {"py", mom.y()}, {"pz", mom.z()}};
	    }

	    for (const std::string& cut : cut_types)
	      if ( cuts[cut] ) {
		counters[cut]++;
		hdiag["tof"][cut]->Fill(t01);
		for (const std::string& det : {"tku", "tkd"})
		  if ( vids[det] >= 0 && (det == "tkd" ? vids["max"] >= 0 : true) ) {
		    for (const std::string& var : vars)
		        hdiag[det+var][cut]->Fill(tpoint[det][var]);
		    if ( det == "tku" ) {
			std::vector<double> vec(4);
			for (size_t k = 0; k < 4; k++)
			    vec[k] = tpoint[det][vars[k]];
			Math::IncrementCovarianceMatrix(dn[cut],
				dcovmats[cut], dmeans[cut], vec);
			for (size_t k = 0; k < 4; k++)
			    Math::IncrementMean(dn[cut], dmeans[cut][k], vec[k]);
			++dn[cut];
		    }
		  }
	      }
	  }

	  if ( !cuts["all"] )
	      continue;

	  // Fetch the virutal hits
          for ( size_t vh = 0; vh < vhits->size(); vh++ ) {

	    // Reject the particles outside of the fiducial surface
	    MAUS::VirtualHit vhit = vhits->at(vh);
	    size_t vid = vhit.GetStationId()-1;
	    MAUS::ThreeVector pos = vhit.GetPosition();
	    if ( !apertures.IsIn(pos.x(), pos.y(), pos.z()) )
		  break;

	    // Reject the particles in TKD that don't make it to station 5 (reflect data)
	    if ( vids["tkd"] > 0 && (int)vid >= vids["tkd"] && vids["max"] < 0 )
		break;

  	    // Reject the particles with a pz below 50, the trackers break down there (TODO)
	    MAUS::ThreeVector mom = vhit.GetMomentum();
	    if ( mom.z() < 50 )
		break;

	    // Fill the virtual samples
	    truth_samples->Fill(spilloff+i, ev, vid,
		pos.x(), pos.y(), pos.z(), mom.x(), mom.y(), mom.z());
	  }
        } // End of the list of MCEvents
      }

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

  // If the diagnostics plots are requested, print them
  if ( globals["diagnostics"] ) {

    // Print the counters to a data file
    ofstream fcounters;
    fcounters.open("diag_counters.txt");
    fcounters << std::left << std::setw(20) << "Truth" << "\n";
    fcounters << std::left << std::setw(20) << "Event" << counters["evt"] << "\n";
    fcounters << std::left << std::setw(20) << "TOF0 SP" << counters["tof0"] << "\n";
    fcounters << std::left << std::setw(20) << "TOF1 SP" << counters["tof1"] << "\n";
    fcounters << std::left << std::setw(20) << "Both TOFs" << counters["none"] << "\n";
    fcounters << std::left << std::setw(20) << "Muon" << counters["pid"] << "\n";
    fcounters << "\n";
    fcounters << std::left << std::setw(20) << "TKU Event" << counters["tkuevt"] << "\n";
    fcounters << std::left << std::setw(20) << "TKU Reference" << counters["ref"] << "\n";
    fcounters << std::left << std::setw(20) << "TKU Quality" << counters["qual"] << "\n";
    fcounters << std::left << std::setw(20) << "TKU Fiducial" << counters["fid"] << "\n";
    fcounters << std::left << std::setw(20) << "TKU Momentum" << counters["mom"] << "\n";
    fcounters << std::left << std::setw(20) << "TKU All" << counters["all"] << "\n";
    fcounters << "\n";
    fcounters << std::left << std::setw(20) << "TKD Event" << counters["tkdevt"] << "\n";
    fcounters << std::left << std::setw(20) << "TKD Reference" << counters["tkdref"] << "\n";
    fcounters << std::left << std::setw(20) << "TKD Quality" << counters["tkdqual"] << "\n";
    fcounters << std::left << std::setw(20) << "TKD Fiducial" << counters["tkdfid"] << "\n";
    fcounters << std::left << std::setw(20) << "TKD Momentum" << counters["tkdmom"] << "\n";
    fcounters << std::left << std::setw(20) << "TKD All" << counters["tkdall"] << "\n";
    fcounters.close();

    // Initialize the RMS ellipses from the recorded means and covariances (+ frames)
    TFile outfile("diag_histograms.root", "RECREATE");
    std::map<std::string, std::map<std::string, TEllipse*>> ellipses;
    std::map<std::string, std::map<std::string, TGraph*>> centres;
    std::map<std::string, TH2F*> frames;
    std::map<std::string, std::vector<double>> mean;
    std::map<std::string, Matrix<double>> mat;
    DGaus gaus;
    double xmin, xmax, ymin, ymax;
    double margin = .5;
    size_t i, j, id;
    for (i = 0; i < 4; i++)
      for (j = i+1; j < 4; j++) {
	id = 0;
        for (const std::string cut : cut_types) {
	  mean[cut] = {dmeans[cut][i], dmeans[cut][j]};
	  mat[cut] = Matrix<double>({{dcovmats[cut][i][i], dcovmats[cut][i][j]},
		     		    {dcovmats[cut][j][i], dcovmats[cut][j][j]}});
	  gaus = DGaus(mean[cut], mat[cut]);
	  ellipses[cut][vars[i]+vars[j]] = (TEllipse*)gaus.Contour2D(.39347)[0];
	  ellipses[cut][vars[i]+vars[j]]->SetLineColor(colors[id]);

	  centres[cut][vars[i]+vars[j]] = new TGraph();
	  centres[cut][vars[i]+vars[j]]->SetPoint(0, mean[cut][0], mean[cut][1]);
	  centres[cut][vars[i]+vars[j]]->SetMarkerStyle(kFullCross);
	  centres[cut][vars[i]+vars[j]]->SetMarkerColor(colors[id]);
	  ++id;
	}

	xmin = mean["none"][0]-(1.+margin)*sqrt(mat["none"][0][0]);
	xmax = mean["none"][0]+(1.+margin)*sqrt(mat["none"][0][0]);
	ymin = mean["none"][1]-(1.+margin)*sqrt(mat["none"][1][1]);
	ymax = mean["none"][1]+(1.+margin)*sqrt(mat["none"][1][1]);
	frames[vars[i]+vars[j]] =
		new TH2F(TString::Format("frame_truth_%s", (vars[i]+vars[j]).c_str()), 
		TString::Format(";%s [%s];%s [%s]", labels[i].c_str(), units[i].c_str(),
		labels[j].c_str(), units[j].c_str()), 1, xmin, xmax, 1, ymin, ymax);
	}

    // Initialize the info box
    std::string data_type = globals["data"] ? "Preliminary" : "[simulation]";
    InfoBox info(data_type, maus_version,
  	globals["run_name"].AsString(), globals["user_cycle"].AsString());
    info.SetPosition("tl");


    // Print the comparision between the different levels of cut
    for (const std::pair<std::string, THStack*>& hist : hsdiag) {
      TCanvas* c = new TCanvas("c", "c", 1200, 800);
      hist.second->Draw("NOSTACK");
      hist.second->Write(hist.second->GetName());
      gPad->SetLogy();
      diag_leg->Draw("SAME");
      info.Draw();
      c->SaveAs(TString::Format("diag_%s_truth.pdf", hist.first.c_str()));
      delete c;
    }

    // Draw the RMS ellipses
    gStyle->SetOptStat(0);
    for (const std::pair<std::string, TEllipse*>& ell : ellipses["none"]) {
      TCanvas* c = new TCanvas("c", "c", 1200, 800);
      frames[ell.first]->Draw();
      frames[ell.first]->Write(TString::Format("frame_truth_%s", ell.first.c_str()));
      
      for (const std::string cut : cut_types) {
	centres[cut][ell.first]->Draw("PSAME");
	centres[cut][ell.first]->Write(TString::Format("centroid_truth_%s_%s",
						cut.c_str(), ell.first.c_str()));
	ellipses[cut][ell.first]->Draw("SAME");
	ellipses[cut][ell.first]->Write(TString::Format("ellipse_truth_%s_%s",
						cut.c_str(), ell.first.c_str()));
      }
      diag_leg->Draw("SAME");
      info.Draw();
      c->SaveAs(TString::Format("diag_%s_truth.pdf", ell.first.c_str()));
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
