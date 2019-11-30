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
  bool upstream = !planeid || planeid == (int)globals["ref_vid"];
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
      cst = 18.7816;
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
     cst = 0.000;
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
      return 10.;
  } else if ( diff == 4 ) { 	// 6 mm beam
      return 15.;
  } else if ( diff == 15 ) { 	// 10 mm beam
      return 20.;
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

TMatrixD TransferMatrix(const TMatrixD& X, const TMatrixD& Y) {

  // Compute the inverse scatter matrix
  size_t d = X.GetNcols();
  size_t n = X.GetNrows();
  TMatrixD S(d, d), Xt(n, d);
  Xt = X;
  Xt.T();
  S = Xt*X;
  S.Invert();

  // Return the least-square solution
  return S*Xt*Y;
}

/** @file  MigrationMatrix.cc
 *
 *  @brief Computes transfer matrices between two beam line positions
 *
 *	   The code deconvolves transfer matrices from one position in the 
 *	   in the beam line to another. It assumes linear transfer.
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

  // If reconstructed data is included, set up particle identification to select muons
  ParticleIdentification particle_id;
  double tofdz = 0.;
  if ( (globals["data"] || globals["recmc"]) && !globals["notofs"] ) {
//    particle_id = ParticleIdentification(globals.GetDataFiles(), 
//						{0, 1}, globals["tof_min"], globals["tof_max"], 1e4);
//    tofdz = particle_id.GetTofDistance();
    tofdz = 7.676337;	// TODO TEMPORARY NO ELECTRONS IN THE SIMULATION...
    Pitch::print(Pitch::info, "Distance between TOFs: "+std::to_string(tofdz)+" m");
  }

  // Initialize cuts
  std::vector<std::string> types;
  std::vector<Cut> ucuts = {None, TOF0SP, TOF1SP, TOF01, Track, Chi2,
				Fiducial, Momentum, EnergyLoss, All};
  std::vector<Cut> dcuts = {None, Track, Chi2, Fiducial, Momentum, All};
  std::map<Cut, bool> upassed, dpassed;
  for (const std::string& type : {"truth", "recmc", "data"}) {
    if ( !globals[type] )
	continue;
    types.push_back(type);
  }

  // Intialize the drawer and the info box
  Beam::Drawer drawer;
  if ( globals["mice"] ) {
    std::string data_type = "[simulation]";
    InfoBox* info = new InfoBox(data_type, globals["maus_version"].AsString(),
  	globals["run_name"].AsString(), globals["user_cycle"].AsString());
    info->SetPosition("tl", .02, .3);
    drawer.SetInfoBox(info);
  }

  // Initialize the vectors
  std::map<std::string, std::vector<std::vector<double>>> pz_t;
  std::vector<std::string> dets = {"tku", "tkd"};
  std::vector<std::string> vars = {"x", "y", "px", "py", "pz"};
  std::vector<std::string> labels = {"x", "y", "p_{x}", "p_{y}", "p_{z}"};
  std::vector<std::string> units = {"mm", "mm", "MeV/c", "MeV/c", "MeV/c"};
//  std::vector<double> lims = {1.999, 1.999, 4.999, 4.999, 9.999};
  size_t d = vars.size();
  std::string name, title;

  std::vector<Vector<double>> samples_ref, samples_tar;
  Vector<double> sample(vars.size());

  // Estimate the resolutions on each of the phase space variables
  std::string maus_version = globals["maus_version"];
  size_t ref_vid, tar_vid;
  ref_vid = 64;				// TKU station 5
  tar_vid = 80;				// Target plane id
  TVector3 pos_ref, mom_ref, pos_tar, mom_tar;
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
    T->SetBranchAddress("Truth", &track_tru);  	// Set the address of data_ptr
    T->SetBranchAddress("Recon", &track_rec);  	// Set the address of data_ptr

    // If the truth or the recmc are missing, cannot proceed
    if ( !track_tru || !track_rec )
	continue;

    // Loop over the tracks
    ProgressBar pbar;
    for (size_t i = 0; i < (size_t)T->GetEntries(); ++i) { // TODO FILE TRUNCATED TODO

      // Display the progress in %
      pbar.GetProgress(i, (size_t)T->GetEntries());

      // Update the track pointed to by *track_tru and *track_rec
      T->GetEntry(i);

      // If there is nothing in the upstream tracker, proceed
      if ( track_tru->mom.find(tar_vid) == track_tru->mom.end() ||
	   track_tru->mom.at(ref_vid).Mag() < 135 ||
	   track_tru->mom.at(ref_vid).Mag() > 145 )
	  continue;

      // Fill the transfer matrix samples
      if ( track_tru->pid.at(ref_vid) == -13 &&
	   track_tru->tk_maxr[0] < 150 ) {

 	pos_ref = track_tru->pos.at(ref_vid);
	mom_ref = track_tru->mom.at(ref_vid);
	sample[0] = pos_ref.x();
	sample[1] = pos_ref.y();
	sample[2] = mom_ref.x();
	sample[3] = mom_ref.y();
	sample[4] = mom_ref.z();
	samples_ref.push_back(sample);

 	pos_tar = track_tru->pos.at(tar_vid);
	mom_tar = track_tru->mom.at(tar_vid);
	sample[0] = pos_tar.x();
	sample[1] = pos_tar.y();
	sample[2] = mom_tar.x();
	sample[3] = mom_tar.y();
	sample[4] = mom_tar.z();
	samples_tar.push_back(sample);
      }
      
    }
  }

  // Initialize the origin and target sample matrices
  size_t n = samples_ref.size();
  TMatrixD X(n, d), Y(n, d);
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < d; j++) {
      X[i][j] = samples_ref[i][j];
      Y[i][j] = samples_tar[i][j];
    }

  // Compute the transfer matrix
  TMatrixD transfer = TransferMatrix(X, Y);
  transfer.Print();

  // Use the matrix to make a radial prediction, compare with the truth
  std::vector<std::vector<double>> res(d);
  std::vector<double> radres(n);
  double rad, rad_prop;
  TMatrixD x(d, 1), y(d, 1);
  for (size_t k = 0; k < n; k++) {
    rad = sqrt(pow(Y[k][0], 2)+pow(Y[k][1], 2));

    for (size_t i = 0; i < d; i++)
	x[i] = X[k][i];
    y = transfer*x;

    for (size_t i = 0; i < d; i++)
	res[i].push_back(y[i][0]-Y[k][i]);

    rad_prop = sqrt(pow(y[0][0], 2)+pow(y[1][0], 2));

    radres[k] = rad-rad_prop;
  }

  gStyle->SetOptStat(1);
  TH1F* hradres = new TH1F("radres", ";R_{diff}-R_{u}^{diff}  [mm]", 100, -100, 100);
  hradres->FillN(radres.size(), &radres[0], NULL);
  TCanvas *c = new TCanvas("c", "c", 1200, 800);
  hradres->Draw();
  c->SaveAs("radres.pdf");
  delete c;

  TFile *outfile = new TFile("transfer_histograms.root", "RECREATE");
  for (size_t i = 0; i < d; i++) {
    TH1F* hres = new TH1F(TString::Format("transfer_res_%s", vars[i].c_str()),
	TString::Format(";%s-#hat{%s} [%s]", labels[i].c_str(),
	labels[i].c_str(), units[i].c_str()),100, -100, 100);
    hres->FillN(res[i].size(), &res[i][0], NULL);
    TCanvas *c = new TCanvas("c", "c", 1200, 800);
    hres->Draw();
    hres->Write(hres->GetName());
    c->SaveAs(TString::Format("%s.pdf", hres->GetName()));
    delete c;
  }

  outfile->Close();
}
