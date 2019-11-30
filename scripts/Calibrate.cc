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

// par[0] = mass, par[1] = thickness*Kstar
double EnergyLossFunction(double *x, double *par) {

  // Convert momentum to energy then C factor
  double E = sqrt(x[0]*x[0]+par[0]*par[0]);
  double C = (E*E+par[0]*par[0])/E;

  // If the particle is stopped, return the particle momentum
  double Estar = (C-par[1])/2.;
  if ( par[0] > Estar )
      return x[0];

  // Find the energy after a thickness par[1] and a stopping constant par[2]
  double Ex = Estar*(1.+sqrt(1.-pow(par[0]/Estar, 2)));

  // Return the momentum
  return x[0]-sqrt(Ex*Ex-par[0]*par[0]);
}

void Limits(const std::vector<double>& samples, double& llim, double& ulim) {

  double mean = Math::Mean(samples);
  double rms = Math::RMS(samples);
  llim = mean-5*rms;
  ulim = mean+5*rms;
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

/** @file  Calibrate.cc
 *
 *  @brief Calibrate to the MICE data
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

  // Load TTree and Track library to ROOT
  gSystem->Load("libTree.so");
  gSystem->Load("../lib/libMiceTrack.so");

  // Extract the time-of-flight and tracker momentum information
  std::vector<size_t> ids;
  std::vector<double> etofs, tofs, moms, momes, dmoms;
  double t01;
  double tofmom, tkumom, meanmom, tofmom_err, tkumom_err;
  for (const std::string& file : globals.GetDataFiles()) {

    // Set up the ROOT file
    TFile data_file(file.c_str());		// Load the MAUS output file
    if ( !data_file.IsOpen() ) {
      Pitch::print(Pitch::error, "Data file not found: "+file);
      return 2;
    }
    Pitch::print(Pitch::info, "Processing "+file);

    TTree *T = (TTree*)data_file.Get("Track");	// Pull out the TTree
    MiceTrack *track_rec = NULL;	 	// A variable to store the data from each track
    T->SetBranchAddress("Recon", &track_rec);  	// Set the address of data_ptr

    // Loop over the tracks
    ProgressBar pbar;
    for (size_t i = 0; i < (size_t)T->GetEntries(); ++i) {

      // Display the progress in %
      pbar.GetProgress(i, (size_t)T->GetEntries());

      // Update the track pointed to by *track	
      T->GetEntry(i);

      // Only apply the calibration to particles that have a time-of-flight
      if ( track_rec->tof_nsp[0] != 1 || track_rec->tof_nsp[1] != 1 )
	  continue;
      t01 = track_rec->t[1]-track_rec->t[0];
      if ( t01 > 24 && t01 < 27 )
          etofs.push_back(t01);

      // Only apply the calibration to particles that have a single track in TKU
      if ( track_rec->tk_nt[0] != 1 ||
	   track_rec->tk_chi2[0] > (double)globals["tk_chi2"]  ||
	   track_rec->tk_maxr[0] > (double)globals["tk_fiducial"] )
	  continue;

      // Only fit the muon peak
      if ( t01 < (double)globals["tof_min"] || t01 > (double)globals["tof_max"] )
	  continue;

      // Skip abhorrent values of momentum
      if ( track_rec->mom[0].Mag() > 300 )
	  continue;

      // Fill the time-of-flight and momentum arrays
      tofs.push_back(t01);    
      moms.push_back(track_rec->mom[0].Mag());
      momes.push_back(track_rec->mome[0].Mag());

      // If there is no downstream track, abort here
      if ( track_rec->tk_nt[1] != 1 ||
	   track_rec->tk_chi2[1] > (double)globals["tk_chi2"]  ||
	   track_rec->tk_maxr[1] > (double)globals["tk_fiducial"] )
	  continue;

      // Fill the downstream momentum array
      ids.push_back(moms.size()-1);
      dmoms.push_back(track_rec->mom[5].Mag());
    }
  }

  // Fit the electron peak, get its position and deduce the distance between the TOFs
  double clight = 299792458; // Speed of light [m/s]
  double llim, ulim;
  Limits(etofs, llim, ulim);
  
  TCanvas* c = new TCanvas("c", "c", 1200, 800);
  TH1F* hpeak = new TH1F("hpeak", ";t_{01}  [ns]", 100, llim, ulim);
  hpeak->FillN(etofs.size(), &etofs[0], NULL);
  hpeak->Draw("");
  c->SaveAs("tofpeak.pdf");
  delete c;

  double peakpos = hpeak->GetMean();
  double dist = peakpos*1e-9*clight;
  Pitch::print(Pitch::info, "Distance between TOF0 and TOF1: "+std::to_string(dist));
//  double peakres = hpeak->GetMeanError();
//  double distres = peakres*1e-9*clight;

  // Compute the estimated TOF momenta
  double mass = 105.66;
  std::vector<double> tofmoms(moms.size()), plosses(moms.size());
  for (size_t i = 0; i < tofs.size(); i++) {
    tofmoms[i] = 105.66/sqrt(pow(clight*tofs[i]*1e-9/dist, 2)-1);
    plosses[i] = tofmoms[i]-moms[i];
  }
  
  // Initialize the fit
  TF1* fitloss = new TF1("fitloss", EnergyLossFunction, 0, 300, 2);
  fitloss->FixParameter(0, mass);
  fitloss->SetParameter(1, 10);
  fitloss->SetParLimits(1, 0, 1000);
  gStyle->SetOptStat(0);

  // Histograms
  TFile *outfile = new TFile("cal_histograms.root", "RECREATE");

  Limits(tofmoms, llim, ulim);
  TH2F* hist = new TH2F("cal_ploss", ";p_{01}  [MeV/c]; #Deltap  [MeV/c]", 100, 175, 225, 40, 40, 90);
  hist->FillN(tofs.size(), &tofmoms[0], &plosses[0], NULL, 1);
  hist->Sumw2();
  hist->Scale(1./hist->GetEntries());

  c = new TCanvas("c", "c", 1200, 800);
  hist->Draw("COLZ");
  hist->Fit(fitloss, "");
  Pitch::print(Pitch::info, "TOF1->TKU constant: "+std::to_string(fitloss->GetParameter(1)));
  hist->Write(hist->GetName());
  c->SaveAs("cal_ploss.pdf");
  delete c;

  // Compute the estimated TOF momenta in TKU
  std::vector<double> residuals(tofs.size());;
  for (size_t i = 0; i < tofs.size(); i++) {
    tofmoms[i] -= fitloss->Eval(tofmoms[i]);
    residuals[i] = moms[i]-tofmoms[i];
  }

  hist = new TH2F("cal_comp", ";p_{01}^{u}  [MeV/c]; p_{u}  [MeV/c]", 100, 100, 200, 100, 100, 200);
  hist->FillN(tofs.size(), &tofmoms[0], &moms[0], NULL, 1);
  hist->Sumw2();
  hist->Scale(1./hist->GetEntries());
      
  c = new TCanvas("c", "c", 1200, 800);
  hist->Draw("COLZ");
  hist->Write(hist->GetName());
  TLine* line = new TLine(100, 100, 200, 200);
  line->Draw("SAME");
  c->SaveAs("cal_comp.pdf");
  delete c;

  TH1F* hres = new TH1F("cal_res", ";p_{u}-p_{01}^{u} [MeV/c]", 100, -50, 50);
  hres->FillN(residuals.size(), &residuals[0], NULL);
  hres->Sumw2();
  hres->Scale(1./hres->GetEntries());

  gStyle->SetOptStat(1);
  Limits(residuals, llim, ulim);
  c = new TCanvas("c", "c", 1200, 800);
  hres->Draw("");
  hres->Write(hres->GetName());
  c->SaveAs("cal_res.pdf");
  delete c;

  // Compute the weighted average estimated TKU momentum
  std::vector<double> meanmoms(ids.size()), plossesd(ids.size());
  for (size_t i = 0; i < ids.size(); i++) {
    // Compute the upstream estimate according to the energy loss
    tofmom = tofmoms[ids[i]];
    tkumom = moms[ids[i]];

    // Take the weighted average momentum estimate
    tkumom_err = momes[ids[i]];
    tofmom_err = Uncertainty(globals["diffuser_thickness"]);
    meanmom = tkumom/pow(tkumom_err, 2) + tofmom/pow(tofmom_err, 2);
    meanmom /= (1./pow(tkumom_err, 2) + 1./pow(tofmom_err, 2));

    meanmoms[i] = meanmom;
    plossesd[i] = meanmom-dmoms[i];
  }

  // Initialize the fit
  fitloss->FixParameter(0, mass);
  fitloss->SetParameter(1, 10);
  fitloss->SetParLimits(1, 0, 1000);
  gStyle->SetOptStat(0);

  Limits(meanmoms, llim, ulim);
  hist = new TH2F("cal_ploss", ";p_{01}  [MeV/c]; #Deltap  [MeV/c]", 100, llim, ulim, 40, 0, 40);
  hist->FillN(meanmoms.size(), &meanmoms[0], &plossesd[0], NULL, 1);
  hist->Sumw2();
  hist->Scale(1./hist->GetEntries());

  c = new TCanvas("c", "c", 1200, 800);
  hist->Draw("COLZ");
  hist->Fit(fitloss, "");
  Pitch::print(Pitch::info, "TKU->TKD constant: "+std::to_string(fitloss->GetParameter(1)));
  hist->Write(hist->GetName());
  c->SaveAs("cal_ploss_tkd.pdf");
  delete c;

  // Compute the estimated momenta in TKD
  std::vector<double> tkdmoms(ids.size()), dresiduals(ids.size());
  for (size_t i = 0; i < dmoms.size(); i++) {
    tkdmoms[i] = meanmoms[i]-fitloss->Eval(meanmoms[i]);
    dresiduals[i] = dmoms[i]-tkdmoms[i];
  }

  hist = new TH2F("cal_comp_tkd", ";p_{u}^{d}  [MeV/c]; p_{d}  [MeV/c]", 100, 100, 200, 100, 100, 200);
  hist->FillN(dmoms.size(), &tkdmoms[0], &dmoms[0], NULL, 1);
  hist->Sumw2();
  hist->Scale(1./hist->GetEntries());
      
  c = new TCanvas("c", "c", 1200, 800);
  hist->Draw("COLZ");
  hist->Write(hist->GetName());
  line->Draw("SAME");
  c->SaveAs("cal_comp_tkd.pdf");
  delete c;

  hres = new TH1F("cal_res_tkd", ";p_{d}-p_{u}^{d} [MeV/c]", 100, -50, 50);
  hres->FillN(dresiduals.size(), &dresiduals[0], NULL);
  hres->Sumw2();
  hres->Scale(1./hres->GetEntries());

  gStyle->SetOptStat(1);
  Limits(residuals, llim, ulim);
  c = new TCanvas("c", "c", 1200, 800);
  hres->Draw("");
  hres->Write(hres->GetName());
  c->SaveAs("cal_res_tkd.pdf");
  delete c;

  outfile->Close();
}
