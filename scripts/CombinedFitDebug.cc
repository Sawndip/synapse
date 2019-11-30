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

// Additional modules
#include "Pitch.hh"
#include "ProgressBar.hh"
#include "InfoBox.hh"
#include "Globals.hh"
#include "GeometryHandler.hh"
#include "DGaus.hh"
#include "Statistics.hh"
#include "Aperture.hh"
#include "MiceTrack.hh"
#include "ParticleIdentification.hh"
#include "EnergyLoss.hh"
#include "Drawer.hh"

enum Cut {None, TOF0SP, TOF1SP, TOF01, Track, Chi2, Fiducial, Momentum, EnergyLoss, All};

TH1F* ProjectionRMSX(TH2F* hist, bool fit)  {

  TH1F* h1d = new TH1F(TString::Format("%s_1d_rms", hist->GetName()),
		       TString::Format("%s;%s;%s", hist->GetTitle(), 
		       hist->GetXaxis()->GetTitle(), hist->GetYaxis()->GetTitle()),
		       hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());

  TF1* fgaus;
  if ( fit )
      fgaus = new TF1("fgaus", "[0]*TMath::Gaus(x, [1], [2])", -1e3, 1e3);

  size_t k;
  TH1D* hproj;
  double rms, error;
  for (k = 0; k < (size_t)hist->GetNbinsX(); k++) {
    hproj = hist->ProjectionY("bin", k+1, k+1);

    if ( !fit ) {
      rms = hproj->GetRMS();
      error = hproj->GetRMSError();
    } else {
      fgaus->SetRange(hproj->GetMean()-hproj->GetRMS(), hproj->GetMean()+hproj->GetRMS());
      fgaus->SetParameters(hproj->GetMaximum(), hproj->GetMean(), hproj->GetRMS());
      hproj->Fit(fgaus, "RQ");
      rms = fabs(fgaus->GetParameter(2));
      error = fgaus->GetParError(2);
    }

//    if ( error > rms )
//	continue;

    h1d->SetBinContent(k+1, rms);
    h1d->SetBinError(k+1, error);
    delete hproj;
  }

  return h1d;
}

std::map<Cut, std::string> CutName = {{None, 		"None"},
				      {TOF0SP, 		"TOF0 SP"},
				      {TOF1SP, 		"TOF1 SP"},
				      {TOF01, 		"Time-of-flight"},
				      {Track, 		"Track"},
				      {Chi2, 		"Chi squared"},
				      {Fiducial,	"Fiducial"},
				      {Momentum, 	"Momentum"},
				      {EnergyLoss, 	"Energy Loss"},
				      {All, 		"All"}};

bool CheckCut(MiceTrack* track,
	      const Cut& cut,
	      const int& planeid) {

  Globals &globals = Globals::GetInstance();
  bool upstream = !planeid || planeid == (int)globals["tku_vid"];
  int tkid = upstream ? 0 : 1;
  switch(cut) {
      case None:
        return true;
      case TOF0SP:
        return track->tof_nsp[0] == 1;
      case TOF1SP:
        return track->tof_nsp[1] == 1;
      case TOF01:
        return (track->tof_nsp[0] == 1 && track->tof_nsp[0] == 1 &&
		track->t[1]-track->t[0] > (double)globals["tof_min"] &&
		track->t[1]-track->t[0] < (double)globals["tof_max"]);
      case Track:
        return track->tk_nt[tkid] == 1;
      case Chi2:
        return (track->tk_nt[tkid] == 1 &&
	        track->tk_chi2[tkid] < (double)globals["tk_chi2"]);
      case Fiducial:
        return (track->tk_nt[tkid] == 1 &&
		track->tk_maxr[tkid] < (double)globals["tk_fiducial"]);
      case Momentum:
	if ( upstream ) {
	  return ( track->tk_nt[0] == 1 &&
		  track->mom[planeid].Mag() > (double)globals["tku_mom_min"] &&
	          track->mom[planeid].Mag() < (double)globals["tku_mom_max"]);
	} else {
	  return (track->tk_nt[1] == 1 &&
		  track->mom[planeid].Mag() > (double)globals["tkd_mom_min"] &&
	          track->mom[planeid].Mag() < (double)globals["tkd_mom_max"]);
	}
      case EnergyLoss:
	return true;
      case All:
	return true;
  }

  return false;
}

void DoEnergyLossMuon(double& p, const double& cst) {

  // Convert momentum to energy then C factor
  double m = 105.66;
  double E = sqrt(p*p+m*m);
  double C = (E*E+m*m)/E;

  // If the particle is stopped, set to 0. and return
  double Estar = (C-cst)/2.;
  if ( m > Estar ) {
    p = 0.;
    return;
  }

  // Find the energy after a thickness par[1] and a stopping constant par[2]
  double Ex = Estar*(1.+sqrt(1.-pow(m/Estar, 2)));

  // Return the momentum
  p = sqrt(Ex*Ex-m*m);
}

double TOFMomentum(const double& tof, const double& dist, const int& diff) {

  // First compute the time-of-flight momentum just before TOF1, assume muon
  double mom = 105.66/sqrt(pow(299792458*tof*1e-9/dist, 2)-1);

  // Apply the energy loss with the appropriate stopping constant
  double cst(0.);
  if ( diff == 0 ) { 		// 3 mm beam
//      cst = 11.5598; 		// 11.5598 +/- 0.00044644 MeV
     cst = 12.0159;
  } else if ( diff == 4 ) { 	// 6 mm beam
//      cst = 17.6816; 		// 17.6816 +/- 0.000368964 MeV
      cst = 18.1207;
  } else if ( diff == 15 ) { 	// 10 mm beam
//      cst = 37.7302; 		// +/- 0.00032925 MeV
      cst = 39.7302; 		// +/- 0.00032925 MeV
  } else {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "This diffuser thickness is not calibrated: "+std::to_string(diff),
	  "TOFMomentum"));
  }

  DoEnergyLossMuon(mom, cst);
  return mom;
}

double DownstreamMomentum(double mom, const int& abs) {

  // Apply the energy loss with the appropriate stopping constant
  double cst(0.);
  if ( abs == 0 ) { 		// No absorber
     cst = 0.700;
  } else if ( abs == 1 ) { 	// LiH absorber
     cst = 8.500;
  } else {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "This absorber setting is not recognized: "+std::to_string(abs),
	  "DownstreamMomentum"));
  }

  DoEnergyLossMuon(mom, cst);
  return mom;
}

double Tolerance(const double& diff) {

  // Return the energy loss tolerance
  if ( diff == 0 ) { 		// 3 mm beam
      return 15.;
  } else if ( diff == 4 ) { 	// 6 mm beam
      return 20.;
  } else if ( diff == 15 ) { 	// 10 mm beam
      return 25.;
  } else {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "This diffuser thickness is not calibrated: "+std::to_string(diff),
	  "Tolerance"));
  }
  return 0.;
}

double Uncertainty(const double& diff) {

  // Return the energy loss tolerance
  if ( diff == 0 ) { 		// 3 mm beam
      return 2.;
  } else if ( diff == 4 ) { 	// 6 mm beam
      return 3.;
  } else if ( diff == 15 ) { 	// 10 mm beam
      return 6.;
  } else {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "This diffuser thickness is not calibrated: "+std::to_string(diff),
	  "Uncertainty"));
  }
  return 0.;
}

struct HistDef {
  std::string name;
  std::string label;
  std::string unit;
  double nbins;
  double llim;
  double ulim;
};

/** @file  ApplyCuts.cc
 *
 *  @brief Applys cuts to the MICE data
 *
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

  // Importer global parameters, check that ROOT data files are provided
  Globals &globals = Globals::GetInstance(argc, argv);
  globals.CheckDataFiles("root");

  // Load TTree and Track library to ROOT
  gSystem->Load("libTree.so");
  gSystem->Load("libMausCpp.so");
  gSystem->Load("../lib/libMiceTrack.so");

  /*// If reconstructed data is included, set up particle identification to select muons
  ParticleIdentification particle_id;
  double tofdz = 0.;
  if ( (globals["data"] || globals["recmc"]) && !globals["notofs"] ) {
    particle_id = ParticleIdentification(globals.GetDataFiles(),
						{0, 1}, globals["tof_min"], globals["tof_max"]);
    tofdz = particle_id.GetTofDistance();
    Pitch::print(Pitch::info, "Distance between TOFs: "+std::to_string(tofdz)+" m");
  }*/

  // Histograms
  std::map<std::string, std::map<std::string, std::vector<double>>> data;
  std::map<std::string, std::map<std::string, TH1F*>> hists;
  std::map<std::string, THStack*> stacks;
  std::vector<std::string> types = {"truth", "recmc", "data"};
  std::map<std::string, int> colors = {{"truth",kBlue}, {"recmc",kRed}, {"data",kGreen+2}};

  stacks["tof"] = new THStack("tof_stack", ";t_{01}  [ns]");
  stacks["tof_mom"] = new THStack("tof_mom_stack", ";p_{01}  [MeV/c]");
  for (const std::string& type : types) {
    hists["tof"][type] = new TH1F(("tof_"+type).c_str(), ";t_{01}  [ns]", 50, 25, 26.5);
    hists["tof"][type]->SetLineColor(colors[type]);
    hists["tof"][type]->SetLineWidth(2);
    stacks["tof"]->Add(hists["tof"][type]);

    hists["tof_mom"][type] = new TH1F(("tof_mom_"+type).c_str(), ";p_{01}  [MeV/c]", 50, 100, 300);
    hists["tof_mom"][type]->SetLineColor(colors[type]);
    hists["tof_mom"][type]->SetLineWidth(2);
    stacks["tof_mom"]->Add(hists["tof_mom"][type]);
  }

  std::map<std::string, TH1F*> residuals;
  residuals["tof"] = new TH1F("res_tof", ";t_{01}^{rec}-t_{01}^{tru}  [ns]", 50, -1, 1);
  residuals["tof_mom"] = new TH1F("res_tof_mom", ";p_{01}^{rec}-p_{01}^{tru}  [MeV/c]", 50, -10, 10);

  TH2F* dpz_dpz = new TH2F("dpz_dpz", ";Recon #sigma_{p_{z}} [MeV/c]; p_{z}^{rec} - p_{z}^{true} [MeV/c]", 
			   20, 0, 10, 50, -20, 20);

  // Loop over the entries and check that the combined fitter is doing its job
  std::string maus_version = globals["maus_version"];
  std::vector<float> array_tru(9), array_rec(16);
  double tofmom(0.), tkumom(0.), tofmom_err(0.), tkumom_err(0.),
	 tkdmom(0.), tkdmom_err(0.), scale(1.), meanmom(0.);
  size_t nplanes, tkuid, tkdid;
  TRandom3 rdmzer(time(NULL));
  for (const std::string& file : globals.GetDataFiles()) {

    // Set up the ROOT file
    TFile data_file(file.c_str());		// Load the MAUS output file
    if ( !data_file.IsOpen() ) {
      Pitch::print(Pitch::error, "Data file not found: "+file);
      return 2;
    }
    Pitch::print(Pitch::info, "Processing "+file);

    TTree *T = (TTree*)data_file.Get("Track");	// Pull out the TTree
    MiceTrack *track_tru = NULL;	 	// A variable to store the truth from each track
    MiceTrack *track_rec = NULL;	 	// A variable to store the rec data from each track
    if ( globals["truth"] )
        T->SetBranchAddress("Truth", &track_tru);  	// Set the address of data_ptr
    if ( globals["recmc"] || globals["data"] )
        T->SetBranchAddress("Recon", &track_rec);  	// Set the address of data_ptr

    // Get the type of reconstructed data
    std::string type;
    if ( track_tru && !track_rec) {
      type = "truth";
    } else {
      type = track_tru ? "recmc" : "data";
    }

    // Get the upstream and downstream reference plane IDs
    tkuid = track_rec ? 0 : globals["tku_vid"];
    tkdid = track_rec ? 5 : globals["tkd_vid"];

    // Loop over the tracks
    ProgressBar pbar;
    for (size_t i = 0; i < .05*(size_t)T->GetEntries(); ++i) {

      // Display the progress in %
      pbar.GetProgress(i, (size_t)T->GetEntries());

      // Update the track pointed to by *track_tru and *track_rec
      T->GetEntry(i);

      // Check that the particle is a muon
//      if ( track_tru->t[1] - track_tru->t[0] < 28 || track_tru->t[1] - track_tru->t[0] > 33 )
//      if ( track_tru->t[1] - track_tru->t[0] < 28 || track_tru->t[1] - track_tru->t[0] > 31.5 )
//      if ( track_tru->t[1] - track_tru->t[0] < 28 || track_tru->t[1] - track_tru->t[0] > 30 )
//      if ( track_tru->t[1] - track_tru->t[0] > 28 )
//	    continue;      

      // Compare the true time-of-flight with the measured time-of-flight
      for ( const std::string& itype : types ) {
	MiceTrack* track = (itype == "truth") ? track_tru : track_rec;
        if ( !track )
	    continue;
	if ( itype != "truth" && type != itype )
	    continue;

	// Fill the time-of-flight information
	double tof = track->t[1] - track->t[0];

        data["tof"][itype].push_back(tof);

	// Fill the time-of-flight muon momentum information
	double beta = 7642.3/(tof*299.792458);
	double gamma = 1./sqrt(1-beta*beta);
        double tof_mom = beta*gamma*105.66;
//	if ( itype == "truth" )
//	    tof_mom = track_tru->mom.at(32).Mag();
        data["tof_mom"][itype].push_back(tof_mom);
      }

/*        if ( track_rec->mom.find(0) != track_rec->mom.end() &&
	     track_tru->mom.find(64) != track_tru->mom.end() ) {
	    dpz_dpz->Fill(track_rec->mome.at(0).z(), track_rec->mom.at(0).z()-track_tru->mom.at(64).z());
        }*/

      // If they the energy lost between TOF1 and TKU is not consistent with a muon seeing
      // the expected amount of material given the diffuser setting, cut out
      /*if ( (upassed[TOF01] || globals["notofs"]) && upassed[Track] ) {

	// Get the reconstructed TOF momentum and the tracker momentum
        // If the simulation does not include the TOFs, smear the truth
	if ( globals["notofs"] ) {
	  if ( track_tru->mom.find(globals["tku_vid"]) != track_tru->mom.end() ) {
	    tofmom = track_tru->mom.at(globals["tku_vid"]).Mag()
		  +rdmzer.Gaus(0., Uncertainty(globals["diffuser_thickness"]));
	  } else {
	    continue;
	  }
        } else {
	  tofmom = TOFMomentum(track_cut->t[1]-track_cut->t[0], tofdz, globals["diffuser_thickness"]);
	}
	tkumom = track_cut->mom[tkuid].Mag();

	// Take the weighted average momentum estimate
        tkumom_err = track_cut->mome[tkuid].Mag();
        tofmom_err = Uncertainty(globals["diffuser_thickness"]);
	meanmom = tkumom/pow(tkumom_err, 2) + tofmom/pow(tofmom_err, 2);
        meanmom /= (1./pow(tkumom_err, 2) + 1./pow(tofmom_err, 2));
      }*/


    }
  }

  for (const std::pair<std::string, THStack*>& stack_pair : stacks) {
    for (const std::string& type : types) {
      hists[stack_pair.first][type]->FillN(data[stack_pair.first][type].size(),
						&data[stack_pair.first][type][0], NULL);
      hists[stack_pair.first][type]->Sumw2();
      hists[stack_pair.first][type]->Scale(1./hists[stack_pair.first][type]->GetEntries());
      hists[stack_pair.first][type]->SetTitle(TString::Format("%s: %0.3f #pm %0.3f ns", type.c_str(), hists[stack_pair.first][type]->GetMean(), hists[stack_pair.first][type]->GetMeanError()));
    }

    TCanvas* c = new TCanvas("c", "c", 1200, 800);
    stack_pair.second->Draw("NOSTACKHIST");
    std::cerr << stack_pair.first << "   " << hists[stack_pair.first]["truth"]->GetMean() << "  " << hists[stack_pair.first]["truth"]->GetMeanError() << std::endl;
    gPad->BuildLegend();
//    gPad->BuildLegend(.125, .675, .45, .875);
    c->SaveAs(("stack_"+stack_pair.first+".pdf").c_str());
    delete c;
  }

  for (const std::pair<std::string, TH1F*>& res_pair : residuals) {
    std::vector<double> diffs;
    for (size_t i = 0; i < data[res_pair.first]["truth"].size(); i++)
	diffs.push_back(data[res_pair.first]["recmc"][i]-data[res_pair.first]["truth"][i]);

    TCanvas* c = new TCanvas("c", "c", 1200, 800);
    res_pair.second->FillN(diffs.size(), &diffs[0], NULL);
    res_pair.second->Draw("");
    c->SaveAs(("res_"+res_pair.first+".pdf").c_str());
    delete c;
  }

  TCanvas* c = new TCanvas("c", "c", 1200, 800);
  gPad->SetLogz();
//  dpz_dpz->Draw("COLZ");
  ProjectionRMSX(dpz_dpz, false)->Draw("E");
  (new TLine(1, 1, 7.5, 7.5))->Draw("SAME");
  c->SaveAs("dpz_dpz.pdf");
  delete c;
}
