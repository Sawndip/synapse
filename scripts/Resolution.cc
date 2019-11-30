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

/** @file  Resolution.cc
 *
 *  @brief Evaluates the tracker resolution
 *
 *	   Evalutes the resolution of the tracker on the phase space variables (x, y, px, py, pz).
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
    particle_id = ParticleIdentification(globals.GetDataFiles(), 
						{0, 1}, globals["tof_min"], globals["tof_max"], 1e4);
    tofdz = particle_id.GetTofDistance();
//    tofdz = 7.676337;	// TODO TEMPORARY NO ELECTRONS IN THE SIMULATION...
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

  // Initialize the resolution graphs
  std::map<std::string, std::vector<std::vector<double>>> pz_t;
  std::map<std::string, std::map<std::string, TH1F*>> residuals;
  std::map<std::string, TH2F*> pzpt, ptpt, ppt, rec, counts;
  std::vector<std::string> dets = {"tku", "tkd"};
  std::vector<std::string> vars = {"x", "y", "px", "py", "pz"};
  std::vector<std::string> labels = {"x", "y", "p_{x}", "p_{y}", "p_{z}"};
  std::vector<std::string> units = {"mm", "mm", "MeV/c", "MeV/c", "MeV/c"};
  std::vector<double> lims = {1.999, 1.999, 4.999, 4.999, 9.999};
  std::string name, title;
  for (const std::string& det : dets) {

    for (size_t i = 0; i < vars.size(); i++) {
      name = det+"_"+vars[i];
      title = ";<#hat{"+labels[i]+"}-"+labels[i]+"> ["+units[i]+"];Density [";
      labels[i].find("p_") != std::string::npos ? title += "("+units[i]+")^{-1}]"
						: title += units[i]+"^{-1}]";
      residuals[det][vars[i]] = new TH1F(name.c_str(), title.c_str(), 50, -lims[i], lims[i]);
     }

    pzpt[det] = new TH2F(("pzpt_"+det).c_str(), ";p_{T}  [MeV/c]; <#hat{p}_{z}-p_{z}> [MeV/c]", 
		12, 0, 60, 100, -50, 50);
    ptpt[det] = new TH2F(("ptpt_"+det).c_str(),
	";p_{T}  [MeV/c]; <#hat{p}_{T}  -p_{T} > [MeV/c]", 12, 0, 60, 100, -50, 50); 
    ppt[det] = new TH2F(("ppt_"+det).c_str(),
	";p_{T}  [MeV/c]; <#hat{p}-p> [MeV/c]", 12, 0, 60, 100, -25, 25);
    ppt[det+"fix"] = new TH2F(("ppt_"+det+"_fix").c_str(),
	";p_{T}  [MeV/c]; <#hat{p}-p> [MeV/c]", 12, 0, 60, 100, -25, 25);
    rec[det] = new TH2F(("rec_"+det).c_str(), ";R [mm]; p_{T} [MeV/c]", 15, 0, 150, 12, 0, 60);
    counts[det] = new TH2F(("counts_"+det).c_str(), ";R [mm]; p_{T} [MeV/c]", 15, 0, 150, 12, 0, 60);
 }

  ppt["tof"] = new TH2F("ppt_tof", ";p_{T}  [MeV/c]; <#hat{p}-p> [MeV/c]", 12, 0, 60, 100, -25, 25);
  ppt["tktk"] = new TH2F("ppt_tktk", ";p_{T}  [MeV/c]; <#hat{p}-p> [MeV/c]", 12, 0, 60, 100, -25, 25);


  // Estimate the resolutions on each of the phase space variables
  std::string maus_version = globals["maus_version"]; // MUST GET IT FROM FILE TODO
  std::vector<float> array_tru(9), array_rec(16);
  double tofmom(0.), tkumom(0.), tofmom_err(0.), tkumom_err(0.), tolerance(0.), meanmom(0.);
  size_t nplanes, tku_vid, tkd_vid, max_vid;
  tku_vid = globals["tku_vid"];
  tkd_vid = globals["tkd_vid"];
  max_vid = globals["max_vid"];
  std::vector<double> tofs, moms;
  TVector3 pos_rec, mom_rec, pos_tru, mom_tru;
  std::vector<double> pxs, pys, rpxs, rpys, pzs, rpzs, tofpzs;
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
      if ( track_tru->mom.find(tku_vid) == track_tru->mom.end() ||
	   track_tru->mom.at(tku_vid).Mag() < 135 ||
	   track_tru->mom.at(tku_vid).Mag() > 145 )
	  continue;

      // Fill the upstream sample if the particle is a muon at the ref. plane 
      double tkdmom = -1;
      if ( track_tru->pid.at(tku_vid) == -13 &&
	   track_tru->tk_maxr[0] < 150 ) {

        // Fill the main counter
 	pos_tru = track_tru->pos.at(tku_vid);
	mom_tru = track_tru->mom.at(tku_vid);
	counts["tku"]->Fill(pos_tru.Pt(), mom_tru.Pt());

	// Fill the residuals if a reconstructed track is present
        if ( track_rec->tk_nt[0] == 1 ) {

	  pos_rec = track_rec->pos.at(0);
	  mom_rec = track_rec->mom.at(0);

	  residuals["tku"]["x"]->Fill(pos_rec.x()-pos_tru.x());
	  residuals["tku"]["y"]->Fill(pos_rec.y()-pos_tru.y());
	  residuals["tku"]["px"]->Fill(mom_rec.x()-mom_tru.x());
	  residuals["tku"]["py"]->Fill(mom_rec.y()-mom_tru.y());
	  residuals["tku"]["pz"]->Fill(mom_rec.z()-mom_tru.z());

	  pzpt["tku"]->Fill(mom_tru.Pt(), mom_rec.z()-mom_tru.z());
	  ptpt["tku"]->Fill(mom_tru.Pt(), mom_rec.Pt()-mom_tru.Pt());
	  ppt["tku"]->Fill(mom_tru.Pt(), mom_rec.Mag()-mom_tru.Mag());

	  rec["tku"]->Fill(pos_tru.Pt(), mom_tru.Pt());

          // If the time-of-flight is available, measure the residual total momentum
          if ( track_rec->tof_nsp[0] == 1 && track_rec->tof_nsp[1] == 1 ) {
	    // Get the reconstructed TOF momentum and the tracker momentum
	    tofmom = TOFMomentum(track_rec->t[1]-track_rec->t[0], tofdz, globals["diffuser_thickness"]);
	    ppt["tof"]->Fill(mom_tru.Pt(), tofmom-mom_tru.Mag());

	    // Evaluate the weighted average upstream momentum
	    tkumom = mom_rec.Mag();
            tkumom_err = track_rec->mome.at(0).Mag();
            tofmom_err = Uncertainty(globals["diffuser_thickness"]);
	    meanmom = tkumom/pow(tkumom_err, 2) + tofmom/pow(tofmom_err, 2);
            meanmom /= (1./pow(tkumom_err, 2) + 1./pow(tofmom_err, 2));
	    ppt["tkufix"]->Fill(mom_tru.Pt(), meanmom-mom_tru.Mag());

	    // Propagate to the downstream tracker
	    tkdmom = DownstreamMomentum(meanmom, globals["absorber"]);
	  }
	}
      } else {
	continue;
      }

      // Fill the downstream sample if the particle is a muon at the last plane
      if ( track_tru->pos.find(max_vid) != track_tru->pos.end() &&
	   track_tru->pid.at(max_vid) == -13 &&
	   track_tru->tk_maxr[1] < 150 ) {

        // Fill the main counter
 	pos_tru = track_tru->pos.at(tkd_vid);
	mom_tru = track_tru->mom.at(tkd_vid);
	counts["tkd"]->Fill(pos_tru.Pt(), mom_tru.Pt());

	// Fill the residuals if a reconstructed track is present
        if ( track_rec->tk_nt[1] == 1 ) {

	  pos_rec = track_rec->pos.at(5);
	  mom_rec = track_rec->mom.at(5);

	  residuals["tkd"]["x"]->Fill(pos_rec.x()-pos_tru.x());
	  residuals["tkd"]["y"]->Fill(pos_rec.y()-pos_tru.y());
	  residuals["tkd"]["px"]->Fill(mom_rec.x()-mom_tru.x());
	  residuals["tkd"]["py"]->Fill(mom_rec.y()-mom_tru.y());
	  residuals["tkd"]["pz"]->Fill(mom_rec.z()-mom_tru.z());

	  pzpt["tkd"]->Fill(mom_tru.Pt(), mom_rec.z()-mom_tru.z());
	  ptpt["tkd"]->Fill(mom_tru.Pt(), mom_rec.Pt()-mom_tru.Pt());
	  ppt["tkd"]->Fill(mom_tru.Pt(), mom_rec.Mag()-mom_tru.Mag());

	  rec["tkd"]->Fill(pos_tru.Pt(), mom_tru.Pt());

	  if ( tkdmom > 0 ) {
	    ppt["tktk"]->Fill(mom_tru.Pt(), tkdmom-mom_tru.Mag());		

	    // Evaluate the weighted average
	    tkumom = tkdmom;
            tkumom_err = 2.;
	    tkdmom = mom_rec.Mag();
            double tkdmom_err = track_rec->mome.at(5).Mag();
	    meanmom = tkumom/pow(tkumom_err, 2) + tkdmom/pow(tkdmom_err, 2);
            meanmom /= (1./pow(tkumom_err, 2) + 1./pow(tkdmom_err, 2));

	    ppt["tkdfix"]->Fill(mom_tru.Pt(), meanmom-mom_tru.Mag()); 
	  }
	}
      }
      
    }
  }

  // Draw the efficiency distributions
  TFile *outfile = new TFile("resol_histograms.root", "RECREATE");
  for (const std::string& det : dets) {
    TCanvas* c = new TCanvas("c", "c", 1200, 800);
    rec[det]->Divide(counts[det]);
    rec[det]->Draw("COLZ");
    rec[det]->Write(rec[det]->GetName());
    c->SaveAs(("resolution_eff_"+det+".pdf").c_str());
    delete c;
  }

  // Draw the residual distributions
//  gStyle->SetTitleSize(.1, "XY");
  std::map<std::string, TH1F*> stack;
  TCanvas *c = new TCanvas("c", "c", 1400, 2000);
  c->Divide(2, 5, 1e-6, 1e-6);
  double max;
  for (size_t j = 0; j < vars.size(); j++) {
    max = 0.;
    for (const std::string& det : dets) {
      residuals[det][vars[j]]->Scale(1./residuals[det][vars[j]]->GetEntries());
      if ( residuals[det][vars[j]]->GetMaximum() > max )
	  max = residuals[det][vars[j]]->GetMaximum();
    }
    max *= 1.1;

    for (size_t i = 0; i < dets.size(); i++) {
      // Set the Pad style
      c->cd(j*2+i+1);
      gPad->SetGridy();
      gPad->SetBottomMargin(0.15);
      if ( !i )
	  gPad->SetRightMargin(0);
      else
	  gPad->SetLeftMargin(0);
      if ( j )
          gPad->SetTopMargin(0.05);

      // Draw
      drawer.SetStyle(residuals[dets[i]][vars[j]], "recmc");
      residuals[dets[i]][vars[j]]->SetMaximum(max);
      residuals[dets[i]][vars[j]]->Draw("HIST");

      // Set the axes style
      residuals[dets[i]][vars[j]]->GetXaxis()->SetTitleSize(.06);
      residuals[dets[i]][vars[j]]->GetXaxis()->SetLabelSize(.06);
      residuals[dets[i]][vars[j]]->GetYaxis()->SetTitleSize(.06);
      residuals[dets[i]][vars[j]]->GetYaxis()->SetTitleOffset(.85);
      residuals[dets[i]][vars[j]]->GetYaxis()->SetLabelSize(.06);
    }
  }

  // Draw info box
  c->cd(1);
  drawer.GetInfoBox()->Draw();

  // Draw tags
  c->cd();
  size_t id = 0;
  for (const std::string& var : {"TKU", "TKD"}) {
    TPaveText* title = new TPaveText(.225+id*.45, .981, .325+id*.45, .999, "NDC");
    title->SetFillStyle(0);
    title->SetLineColor(kBlack);
    TText* text = title->AddText(var.c_str());
    text->SetTextAlign(22);
    text->SetTextSize(0.03);
    title->Draw("SAME");
    id++;
  }

  c->SaveAs("resolutions.pdf");
  delete c;

  // Draw the residual distributions as a function of pT
  for (const std::string& det : {"tku", "tkd", "tof", "tktk", "tkufix", "tkdfix"}) {

    c = new TCanvas("c", "c", 1200, 800);
    ppt[det]->Draw("COLZ");
    ppt[det]->Write(ppt[det]->GetName());
    c->SaveAs(("resolution_ppt_2d_"+det+".pdf").c_str());
    delete c;

    TH1F* ppt_1dm = drawer.ProjectionMeanX(ppt[det], true);
    TH1F* ppt_1d = drawer.ProjectionRMSX(ppt[det], true);
    ppt_1dm->Write(ppt_1dm->GetName());
    ppt_1d->Write(ppt_1d->GetName());
    for (size_t i = 0; i < ppt_1dm->GetNbinsX(); i++) {
      ppt_1d->SetBinError(i+1, ppt_1d->GetBinContent(i+1));
      ppt_1d->SetBinContent(i+1, ppt_1dm->GetBinContent(i+1));
    }
    c = new TCanvas("c", "c", 1200, 800);

    ppt_1d->SetLineWidth(2);
    ppt_1d->SetLineColor(kBlue);
    ppt_1d->SetFillColorAlpha(kBlue, .25);
    ppt_1d->Draw("LE2");
    ppt_1dm->SetMarkerStyle(21);
    ppt_1dm->SetMarkerSize(2);
    ppt_1dm->SetMarkerColor(4);
    ppt_1dm->SetLineWidth(2);
    ppt_1dm->SetLineColor(kBlue);
    ppt_1dm->SetFillColorAlpha(kBlue, .5);
    ppt_1dm->Draw("LE2SAME");

    c->SaveAs(("resolution_ppt_"+det+".pdf").c_str());
    delete c;

    if ( det != "tku" && det != "tkd" )
	continue;

    TCanvas* c = new TCanvas("c", "c", 1200, 800);
    pzpt[det]->Draw("COLZ");
    pzpt[det]->Write(pzpt[det]->GetName());
    c->SaveAs(("resolution_pzpt_2d_"+det+".pdf").c_str());
    delete c;

    c = new TCanvas("c", "c", 1200, 800);
    ptpt[det]->Draw("COLZ");
    ptpt[det]->Write(ptpt[det]->GetName());
    c->SaveAs(("resolution_ptpt_2d_"+det+".pdf").c_str());
    delete c;

    // Create a profile of the mean as a function of pT
    TH1F* pzpt_1dm = drawer.ProjectionMeanX(pzpt[det], true);
    TH1F* pzpt_1d = drawer.ProjectionRMSX(pzpt[det], true);
    pzpt_1dm->Write(pzpt_1dm->GetName());
    pzpt_1d->Write(pzpt_1d->GetName());
    for (size_t i = 0; i < pzpt_1dm->GetNbinsX(); i++) {
      pzpt_1d->SetBinError(i+1, pzpt_1d->GetBinContent(i+1));
      pzpt_1d->SetBinContent(i+1, pzpt_1dm->GetBinContent(i+1));
    }
    c = new TCanvas("c", "c", 1200, 800);

    pzpt_1d->SetLineWidth(2);
    pzpt_1d->SetLineColor(kBlue);
    pzpt_1d->SetFillColorAlpha(kBlue, .25);
    pzpt_1d->Draw("LE2");
    pzpt_1dm->SetMarkerStyle(21);
    pzpt_1dm->SetMarkerSize(2);
    pzpt_1dm->SetMarkerColor(4);
    pzpt_1dm->SetLineWidth(2);
    pzpt_1dm->SetLineColor(kBlue);
    pzpt_1dm->SetFillColorAlpha(kBlue, .5);
    pzpt_1dm->Draw("LE2SAME");

    c->SaveAs(("resolution_pzpt_"+det+".pdf").c_str());
    delete c;

    TH1F* ptpt_1dm = drawer.ProjectionMeanX(ptpt[det], true);
    TH1F* ptpt_1d = drawer.ProjectionRMSX(ptpt[det], true);
    ptpt_1dm->Write(ptpt_1dm->GetName());
    ptpt_1d->Write(ptpt_1d->GetName());
    for (size_t i = 0; i < ptpt_1dm->GetNbinsX(); i++) {
      ptpt_1d->SetBinError(i+1, ptpt_1d->GetBinContent(i+1));
      ptpt_1d->SetBinContent(i+1, ptpt_1dm->GetBinContent(i+1));
    }
    c = new TCanvas("c", "c", 1200, 800);

    ptpt_1d->SetLineWidth(2);
    ptpt_1d->SetLineColor(kBlue);
    ptpt_1d->SetFillColorAlpha(kBlue, .25);
    ptpt_1d->Draw("LE2");
    ptpt_1dm->SetMarkerStyle(21);
    ptpt_1dm->SetMarkerSize(2);
    ptpt_1dm->SetMarkerColor(4);
    ptpt_1dm->SetLineWidth(2);
    ptpt_1dm->SetLineColor(kBlue);
    ptpt_1dm->SetFillColorAlpha(kBlue, .5);
    ptpt_1dm->Draw("LE2SAME");

    c->SaveAs(("resolution_ptpt_"+det+".pdf").c_str());
    delete c;
  }

  outfile->Close();
}
