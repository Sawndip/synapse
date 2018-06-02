// Cpp includes
#include <ctime>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <regex>
#include <map>

// Root includes
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TObject.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TF1.h"
#include "TVector3.h"
#include "TMath.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TRandom3.h"
#include "TColor.h"

// Additional modules
#include "ProgressBar.hh"
#include "InfoBox.hh"
#include "Bunch.hh"
#include "GeometryHandler.hh"
#include "Vector.hh"
#include "Voronoi.hh"
#include "Globals.hh"
#include "DSpiral.hh"

/** @file  SystematicCorrections.cc
 *
 *  @brief Produces systematics coorections to the MICE data.
 *
 *	   Algorithm that produces systematic corrections to the MICE data provided with
 *	   with a MAUS simulation of the beam.
 **/

/** @brief	Main function
 *
 *  @param	argc		Number of command line arguments
 *  @param	argv		List of command line aguments
 */
int main(int argc, char ** argv) {

  // Reroot the info messages to the terminal for now
  Pitch::setAnOutput(Pitch::info, std::cerr);

  // Systematics algorithm global parameters (stored in globals)
  Globals &globals = Globals::GetInstance(argc, argv);
  if ( !globals.GetDataFiles().size() ) {
    Pitch::print(Pitch::error, "No data file specified, please execute in the following fashion:");
    Pitch::print(Pitch::error, "./syst [options] import0.root [... importN.root]");
    return 2;
  } else {
    for (const std::string& file : globals.GetDataFiles())
      if ( std::string(file).find(".root") > std::string(file).size() ) {
        Pitch::print(Pitch::error, "Invalid data file type: "+std::string(file));
	return 2;
      }
  }

  // Get the run name from the first file name (remove superfluous with regex)
  std::string run_name = std::regex_replace(globals.GetDataFiles()[0],
				std::regex("import_|(.*/)|(.root)"), std::string(""));
  Pitch::print(Pitch::info, "Run name: "+run_name);

  // Initialize the geometry handler with the relevant detectors
  std::vector<std::string> dets = {"tku", "tkd"};	// Detectors
  std::map<std::string, std::vector<double>> tk_st_z;
  std::vector<double> stations_z = {-550., -350., -100., 200., 550.};
  std::map<std::string, DetectorGlobals> tk_globals;
  try {
    for (const std::string& det : dets) {
      tk_globals[det] = (det == "tku") ? DetectorGlobals(globals["geometry_filename"], "Tracker0")
				       : DetectorGlobals(globals["geometry_filename"], "Tracker1");
      tk_st_z[det].resize(5);
      for (size_t st = 0; st < 5; st++)
          tk_st_z[det][st] = (det == "tku") ? tk_globals[det].z() - stations_z[st]
					    : tk_globals[det].z() + stations_z[st];
    }
  } catch ( Exceptions::Exception& e ) {
    if ( globals["recmc"] || globals["data"] ) {
      Pitch::print(Pitch::error, "Could not get the tracker positions\n"+std::string(e.what()));
      return 2;    
    } 
  }

  // Set the Stat box to be along the top right corner of the TPad in the TCanvas
  gSystem->Load("libTree");
  gStyle->SetOptStat(0);
  gStyle->SetStatX(.9);
  gStyle->SetStatY(.9);
  gStyle->SetOptStat(2210);
  gStyle->SetOptFit(1);

  std::vector<double> stops = {0., 0.125, .25, .375, .5, .625, .75, .875, 1.};
  std::vector<double> red = {0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764};
  std::vector<double> green = {0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832};
  std::vector<double> blue = {0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539};
  TColor::CreateGradientColorTable(9, &stops[0], &red[0], &green[0], &blue[0], 255, 1);

  std::map<std::string, int> colors = {{"truth",kBlue+2}, {"recmc",kRed+2}, {"data",kGreen+2}};
  std::map<std::string, int> markers = {{"truth",33}, {"recmc",20}, {"data",21}};

  // Variables involved in the analysis
  std::vector<std::string> vars = {"x", "y", "px", "py", "pz"};		    // Emittance pars tags
  std::vector<std::string> labels = {"x", "y", "p_{x}", "p_{y}", "p_{z}"};  // Emittance pars names
  std::vector<std::string> units = {"mm", "mm", "MeV/c", "MeV/c", "MeV/c"}; // Emittance units
  std::vector<double> llims = {-250, -250, -100, -100, 0};	            // Range lower limits
  std::vector<double> ulims = {250, 250, 100, 100, 250};		    // Range upper limits

  // Containers for the true samples, one per virtual plane
  size_t nvirtuals = (size_t)globals["nvirtuals"];
  std::vector<std::map<std::string,
	std::vector<double>>> tsamples(nvirtuals);	// True samples (x, y, px, py, pz)
  std::vector<std::map<std::string,
	std::vector<double>>> utsamples(nvirtuals);	// True samples (x, y, px, py, pz)
  std::vector<size_t> ttransmissions(nvirtuals); 	// True transmission
  std::vector<double> virtuals_z(nvirtuals);		// Position of the virtual planes
  std::map<std::string, std::vector<int>> tk_vids;	// Virtual IDs of tracker stations
  for (const std::string &det : dets)
      tk_vids[det] = std::vector<int>(5, -1);

  // Containers for the reconstructed MC samples and real data, one per tracker station
  std::vector<std::string> types;			// Types of reconstructed data
  std::map<std::string,std::map<std::string,
	std::vector<std::map<std::string, std::vector<double>>>>> rsamples, rerrors;
  std::map<std::string, std::map<std::string, std::vector<size_t>>> rtransmissions;
  for (const std::string &type : {"recmc", "data"}) {
    for (const std::string &det : dets) {
      rsamples[type][det].resize(5);
      rerrors[type][det].resize(5);
      rtransmissions[type][det].resize(5);
    }
  }

  // Fill the sample maps from the imported ROOT file
  std::string maus_version;
  for (const std::string& file : globals.GetDataFiles()) {
    // Set up the ROOT file and data pointer
    TFile data_file(file.c_str());		// Load the MAUS output file
    if ( !data_file.IsOpen() ) {
      Pitch::print(Pitch::error, "Data file not found: "+file);
      return 2;
    }
    Pitch::print(Pitch::info, "Processing "+file);

    // Get the MAUS version if the file was produced in MAUS
    if ( data_file.GetListOfKeys()->Contains("MausVersion") ) {
      maus_version = data_file.Get("MausVersion")->GetTitle();
      maus_version = maus_version.substr(maus_version.size()-5);
      Pitch::print(Pitch::info, "MAUS version: "+maus_version);
    }

    // Import the truth if it is provided
    if ( data_file.GetListOfKeys()->Contains("Truth") && globals["truth"] ) {
      TNtuple* truth_samples = (TNtuple*)data_file.Get("Truth");
      float* ntuple;
      size_t vid;
      TVector3 mom, pos;

      Pitch::print(Pitch::info, "Processing the truth");
      ProgressBar pbar;
      for (size_t i = 0; i < (size_t)truth_samples->GetEntries(); ++i) {

	// Fetch the variables from the Ntuple
	truth_samples->GetEntry(i);
	ntuple = truth_samples->GetArgs();
	vid = ntuple[2];
	pos = TVector3(ntuple[3], ntuple[4], ntuple[5]);
	mom = TVector3(ntuple[6], ntuple[7], ntuple[8]);
	
	// Fill the true samples
        tsamples[vid]["x"].push_back(pos.x());
        tsamples[vid]["y"].push_back(pos.y());
        tsamples[vid]["px"].push_back(mom.x());
        tsamples[vid]["py"].push_back(mom.y());
        tsamples[vid]["pz"].push_back(mom.z());
	ttransmissions[vid]++;

	// If it's the first time this V plane is encountered, add it
	// Retain the virtual IDs of the plane corresponding to the tracker stations
	if ( !virtuals_z[vid] ) {
	  virtuals_z[vid] = pos.z();
	  if ( globals["mice"] ) {
	    for (const std::string& det : dets)
	      for (size_t st = 0; st < 5; st++)
	        if ( fabs(tk_st_z[det][st] - pos.z()) < 0.5 ) {
		    tk_vids[det][st] = vid;
		    break;
	        }
	  }
	}

        // Display the progress in %
        pbar.GetProgress(i, (size_t)truth_samples->GetEntries());
      }
    }

    // Import the uncut truth if it is provided
    if ( data_file.GetListOfKeys()->Contains("UncutTruth") && globals["truth"] ) {
      TNtuple* truth_samples = (TNtuple*)data_file.Get("UncutTruth");
      float* ntuple;
      size_t vid;
      TVector3 mom, pos;

      Pitch::print(Pitch::info, "Processing the uncut truth");
      ProgressBar pbar;
      for (size_t i = 0; i < (size_t)truth_samples->GetEntries(); ++i) {

	// Fetch the variables from the Ntuple
	truth_samples->GetEntry(i);
	ntuple = truth_samples->GetArgs();
	vid = ntuple[2];
	pos = TVector3(ntuple[3], ntuple[4], ntuple[5]);
	mom = TVector3(ntuple[6], ntuple[7], ntuple[8]);
	
	// Fill the true samples
        utsamples[vid]["x"].push_back(pos.x());
        utsamples[vid]["y"].push_back(pos.y());
        utsamples[vid]["px"].push_back(mom.x());
        utsamples[vid]["py"].push_back(mom.y());
        utsamples[vid]["pz"].push_back(mom.z());

        // Display the progress in %
        pbar.GetProgress(i, (size_t)truth_samples->GetEntries());
      }
    }

    // Import the reconstructed MC/data if it is provided
    for (const std::string& type_uc : {"RecMC", "Data"}) {
      std::string type = type_uc;
      std::transform(type.begin(), type.end(), type.begin(), ::tolower);
      if ( data_file.GetListOfKeys()->Contains(type_uc.c_str()) && globals[type] ) {
        types.push_back(type);
        TNtuple* rec_samples = (TNtuple*)data_file.Get(type_uc.c_str());
        float* ntuple;
	std::string det;
        size_t st;
        TVector3 mom, pos, mome, pose;

        Pitch::print(Pitch::info, "Processing the "+type);
        ProgressBar pbar;
        for (size_t i = 0; i < (size_t)rec_samples->GetEntries(); ++i) {

	  // Fetch the variables from the Ntuple
	  rec_samples->GetEntry(i);
	  ntuple = rec_samples->GetArgs();
	  det = ntuple[2] ? "tkd" : "tku";
	  st = ntuple[3];
	  pos = TVector3(ntuple[4], ntuple[5], ntuple[6]);
	  mom = TVector3(ntuple[7], ntuple[8], ntuple[9]);
	  pose = TVector3(ntuple[10], ntuple[11], ntuple[12]);
	  mome = TVector3(ntuple[13], ntuple[14], ntuple[15]);

	  // Fill the reconstructed samples
          rsamples[type][det][st]["x"].push_back(pos.x());
          rsamples[type][det][st]["y"].push_back(pos.y());
          rsamples[type][det][st]["px"].push_back(mom.x());
          rsamples[type][det][st]["py"].push_back(mom.y());
          rsamples[type][det][st]["pz"].push_back(mom.z());

          rerrors[type][det][st]["x"].push_back(pose.x());
          rerrors[type][det][st]["y"].push_back(pose.y());
          rerrors[type][det][st]["px"].push_back(mome.x());
          rerrors[type][det][st]["py"].push_back(mome.y());
          rerrors[type][det][st]["pz"].push_back(mome.z());
	  rtransmissions[type][det][st]++;

          // Display the progress in %
          pbar.GetProgress(i, (size_t)rec_samples->GetEntries());
        }
      }
    }

  } // End of the list of files

  // Import the set of amplitudes upstream and downstream in the simulatio and the recmc
  std::map<std::string, std::vector<double>> tamps, utamps, ramps;
  Beam::Bunch tbeam, utbeam, rbeam;
  for (const std::string set : {"up", "down"}) {
    tbeam = (set == "up") ? Beam::Bunch(tsamples[(size_t)globals["tku_vid"]]) :
			    Beam::Bunch(tsamples[(size_t)globals["tkd_vid"]]);
    utbeam = (set == "up") ? Beam::Bunch(utsamples[(size_t)globals["tku_vid"]]) :
			     Beam::Bunch(utsamples[(size_t)globals["tkd_vid"]]);
    rbeam = (set == "up") ? Beam::Bunch(rsamples["recmc"]["tku"][0]) :
			    Beam::Bunch(rsamples["recmc"]["tkd"][0]);
    if ( globals["corrected"] ) {
      tbeam.SetCorrectedAmplitudes();
      utbeam.SetCorrectedAmplitudes();
      rbeam.SetCorrectedAmplitudes();
    }
    tamps[set] = tbeam.Amplitudes();
    utamps[set] = utbeam.Amplitudes();
    ramps[set] = rbeam.Amplitudes();
  }

  // Draw and save migration matrix for the upstream and the downstream samples
  TFile *outfile = new TFile(TString::Format("%s_syst.root", run_name.c_str()), "RECREATE");
  for (const std::string set : {"up", "down"}) {
    // Set the migration matrix
    TH2F* hresol = new TH2F(TString::Format("syst_res_%s", set.c_str()),
			    ";Reconstructed A_{#perp}  [mm];True A_{#perp}  [mm]",
	  		    20, 0, 100, 20, 0, 100);
    hresol->FillN(tamps[set].size(), &ramps[set][0], &tamps[set][0], NULL, 1);

    // Normalize the bins to 1
    size_t nentries;
    for (size_t i = 0; i < (size_t)hresol->GetNbinsX(); i++) {
      nentries = 0;
      for (size_t j = 0; j < (size_t)hresol->GetNbinsY(); j++)
          nentries += (size_t)hresol->GetBinContent(i+1, j+1);

      if ( nentries )
        for (size_t j = 0; j < (size_t)hresol->GetNbinsY(); j++)
            hresol->SetBinContent(i+1, j+1, hresol->GetBinContent(i+1, j+1)/nentries);
    }

    // Draw and save
    TCanvas* c = new TCanvas("c", "c", 1200, 800);
    hresol->Draw("COLZ");
    hresol->Write(hresol->GetName());
    c->SaveAs(TString::Format("syst_res_%s.pdf", set.c_str()));
    delete c;
  }

  // Draw the efficiency as a function of the true amplitude
  for (const std::string set : {"up", "down"}) {
    // Set the ammplitude histograms
    TH1F* heff_rec = new TH1F(TString::Format("syst_eff_rec_%s", set.c_str()),
			  	";True A_{#perp}  [mm];A_{#perp}/A_{#perp}", 20, 0, 100);
    TH1F* heff = new TH1F(TString::Format("syst_eff_%s", set.c_str()),
			  	";True A_{#perp}  [mm];A_{#perp}/A_{#perp}", 20, 0, 100);
    heff_rec->FillN(tamps[set].size(), &tamps[set][0], NULL);
    heff->FillN(utamps[set].size(), &utamps[set][0], NULL);

    // Compute the ratio
    heff_rec->Sumw2();
    heff->Sumw2();
    heff->Divide(heff_rec);
    heff->SetMarkerStyle(20);
    heff->SetLineWidth(2);
    heff->SetLineColor(kBlack);

    // Draw and save
    TLine *line = new TLine(0, 1, 100, 1);
    line->SetLineColorAlpha(kBlack, .5);

    TCanvas* c = new TCanvas("c", "c", 1200, 800);
    heff->Draw("EP");
    heff->SetMaximum(2); // TODO
    line->Draw("SAME");
    heff->Write(heff->GetName());
    c->SaveAs(TString::Format("syst_eff_%s.pdf", set.c_str()));
    delete c;
  }

  // Draw the uncertainty on the amount of particles in each bin due to the field
  double field_err = 0.001; // Field uncertainty from the hall probes (0.1%)
  for (const std::string set : {"up", "down"}) {
    // Set up the amplitudes of the data
    std::string det = (set == "up") ? "tku" : "tkd";
    rbeam = Beam::Bunch(rsamples["data"][det][0]);
    if ( globals["corrected"] )
        rbeam.SetCorrectedAmplitudes();
    ramps[set] = rbeam.Amplitudes();

    // Set up the baseline histogram
    TH1F* hbase = new TH1F(TString::Format("syst_base_%s", set.c_str()),
			  	";A_{#perp}  [mm]", 20, 0, 100);
    hbase->FillN(ramps[set].size(), &ramps[set][0], NULL);

    // Scale all the transverse momentum up by field_err
    for (size_t i = 0; i < rsamples["data"][det][0]["px"].size(); i++) {
	rsamples["data"][det][0]["px"][i] *= (1.+field_err);
	rsamples["data"][det][0]["py"][i] *= (1.+field_err);
    }
    rbeam = Beam::Bunch(rsamples["data"][det][0]);
    if ( globals["corrected"] )
        rbeam.SetCorrectedAmplitudes();
    ramps[set] = rbeam.Amplitudes();
    
    // Set up the +1sigma histogram
    TH1F* hplus = new TH1F(TString::Format("syst_plus_%s", set.c_str()),
			  	";A_{#perp}  [mm]", 20, 0, 100);
    hplus->FillN(ramps[set].size(), &ramps[set][0], NULL);

    // Scale all the transverse momentum down by field_err
    for (size_t i = 0; i < rsamples["data"][det][0]["px"].size(); i++) {
	rsamples["data"][det][0]["px"][i] *= (1.-field_err)/(1.+field_err);
	rsamples["data"][det][0]["py"][i] *= (1.-field_err)/(1.+field_err);
    }
    rbeam = Beam::Bunch(rsamples["data"][det][0]);
    if ( globals["corrected"] )
        rbeam.SetCorrectedAmplitudes();
    ramps[set] = rbeam.Amplitudes();
    
    // Set up the -1sigma histogram
    TH1F* hminus = new TH1F(TString::Format("syst_minus_%s", set.c_str()),
			  	";A_{#perp}  [mm]", 20, 0, 100);
    hminus->FillN(ramps[set].size(), &ramps[set][0], NULL);

    // Compute the ratios
    hbase->Sumw2();
    hplus->Sumw2();
    hminus->Sumw2();
    hplus->Divide(hbase);
    hminus->Divide(hbase);

    hplus->SetMarkerStyle(20);
    hplus->SetLineWidth(2);
    hplus->SetLineColor(kRed+1);
    hminus->SetMarkerStyle(20);
    hminus->SetLineWidth(2);
    hminus->SetLineColor(kBlue+1);

    // Draw and save
    TLine *line = new TLine(0, 1, 100, 1);
    line->SetLineColorAlpha(kBlack, .5);

    TCanvas* c = new TCanvas("c", "c", 1200, 800);
    hplus->Draw("");
    hminus->Draw("SAME");
    hplus->Write(hplus->GetName());
    hminus->Write(hminus->GetName());
    c->SaveAs(TString::Format("syst_scale_%s.pdf", set.c_str()));
    delete c;
  }
  outfile->Close();

  Pitch::print(Pitch::info, "Moving the amplitude graphs to "+run_name+"/amps");
  std::string sysCmd = "mkdir -p "+run_name+"/amps; mv syst*.pdf *_syst.root "+run_name+"/amps; done";
  if ( std::system(sysCmd.c_str()) )
      Pitch::print(Pitch::error, "Couldn't move systematics plots");
}
