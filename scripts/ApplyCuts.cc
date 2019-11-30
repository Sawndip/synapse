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
      cst = 11.5598; 		// 11.5598 +/- 0.00044644 MeV
//     cst = 12.0159;
  } else if ( diff == 4 ) { 	// 6 mm beam
      cst = 17.6816; 		// 17.6816 +/- 0.000368964 MeV
//      cst = 18.1207;
  } else if ( diff == 15 ) { 	// 10 mm beam
      cst = 37.7302; 		// +/- 0.00032925 MeV
//      cst = 39.7302; 		// +/- 0.00032925 MeV
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
  gSystem->Load("libMausCpp.so");
  gSystem->Load("../lib/libMiceTrack.so");

  // If reconstructed data is included, set up particle identification to select muons
  ParticleIdentification particle_id;
  double tofdz = 0.;
  if ( (globals["data"] || globals["recmc"]) && !globals["notofs"] ) {
    particle_id = ParticleIdentification(globals.GetDataFiles(),
						{0, 1}, globals["tof_min"], globals["tof_max"]);
    tofdz = particle_id.GetTofDistance();
//    tofdz = 7.661326;
//    tofdz = 7.662702;
//    tofdz = 7.676337;	// TODO TEMPORARY NO ELECTRONS IN THE SIMULATION...
//    tofdz = 7.660134;
//    tofdz = 7.651550;
    Pitch::print(Pitch::info, "Distance between TOFs: "+std::to_string(tofdz)+" m");
  }

  // Set up and draw the apertures according to the requested model
  Beam::Aperture apertures;
  if ( (double)globals["aperture"] > 0 ) {
    const GeometryHandler& geoh = globals.GetGeometryHandler();
    apertures.Add(geoh["tku"].z()-1e3, geoh["tkd"].z()+1e3, globals["aperture"]);
  } else {
    apertures.SetMICEDefault(globals["geometry_filename"]);
  }
  apertures.Draw();

  // Initialize cuts and counters
  std::vector<std::string> types;
  std::vector<Cut> ucuts = {None, TOF0SP, TOF1SP, TOF01, Track, Chi2,
				Fiducial, Momentum, EnergyLoss, All};
  std::vector<Cut> dcuts = {None, Track, Chi2, Fiducial, Momentum, All};
  std::map<std::string, std::map<Cut, size_t>> ucounters, dcounters;
  std::map<Cut, bool> upassed, dpassed;
  for (const std::string& type : {"truth", "recmc", "data"}) {
    if ( !globals[type] )
	continue;
    types.push_back(type);
    for (const Cut& cut : ucuts)
        ucounters[type][cut] = 0;
    for (const Cut& cut : dcuts)
        dcounters[type][cut] = 0;
  }

  // If diagnostics are requested, initialize histograms
  std::map<Cut, HistDef> uhistdefs, dhistdefs;
  std::map<Cut, std::map<std::string, TH1F*>> uhists, dhists;
  if ( globals["diagnostics"] ) {

    uhistdefs[TOF0SP] = {"tof0_nsp", "Number of TOF0 space points", "", 5, -0.5, 4.5};
    uhistdefs[TOF1SP] = {"tof1_nsp", "Number of TOF1 space points", "", 5, -0.5, 4.5};
    uhistdefs[TOF01] = {"tof01", "t_{01}", "ns", 80, 0., 8.};
    uhistdefs[Track] = {"tku_nt", "Number of tracks in TKU", "", 5, -0.5, 4.5};
    uhistdefs[Chi2] = {"tku_chi2", "#chi^{2}_{u}/ndf", "", 80, 0, 8};
    uhistdefs[Fiducial] = {"tku_fid", "R_{max}^{u}", "mm", 100, 0, 200};
    uhistdefs[Momentum] = {"tku_mom", "p_{u}", "MeV/c", 100, 100, 180};
    uhistdefs[EnergyLoss] = {"eloss", "p_{u}-p_{01}^{u}", "MeV/c", 100, -50, 50};

    dhistdefs[Track] = {"tkd_nt", "Number of tracks in TKD", "", 5, -0.5, 4.5};
    dhistdefs[Chi2] = {"tkd_chi2", "#chi^{2}_{d}/ndf", "", 80, 0, 8};
    dhistdefs[Fiducial] = {"tkd_fid", "R_{max}^{d}", "mm", 100, 0, 200};
    dhistdefs[Momentum] = {"tkd_mom", "p_{d}", "MeV/c", 100, 50, 210};

    for (const std::string& type : types) {

      // Initialize the upstream diagnostics histograms
      for (const std::pair<Cut, HistDef>& el : uhistdefs) {
	const Cut &cut = el.first;
	const HistDef &def = el.second;
	std::string axis_title = !def.unit.size() ? ";"+def.label : ";"+def.label+" ["+def.unit+"]";
	uhists[cut][type] = new TH1F((def.name+"_"+type).c_str(), axis_title.c_str(),
								def.nbins, def.llim, def.ulim);
      }

      // Initialize the downstream diagnostics histograms
      for (const std::pair<Cut, HistDef>& el : dhistdefs) {
	const Cut &cut = el.first;
	const HistDef &def = el.second;
	std::string axis_title = !def.unit.size() ? ";"+def.label : ";"+def.label+" ["+def.unit+"]";
	dhists[cut][type] = new TH1F((def.name+"_"+type).c_str(), axis_title.c_str(),
								def.nbins, def.llim, def.ulim);
      }
    }
  }

  // Initialize the output
  std::string run_name = std::regex_replace(globals.GetDataFiles()[0],
	std::regex("import_|(.*/)|(.root)|_[0-9]{5}"), std::string(""));
  std::string out_file = "import_"+run_name+".root";
  TFile *out = new TFile(out_file.c_str(), "RECREATE");

  TNtuple* truth_samples = new TNtuple("truth_ntuple", "",
	"SpillID:EventID:VirtualPlaneID:x:y:z:px:py:pz");
  TNtuple* uncut_truth_samples = new TNtuple("uncut_truth_ntuple", "",
	"SpillID:EventID:VirtualPlaneID:x:y:z:px:py:pz");
  TNtuple* recmc_samples = new TNtuple("recmc_ntuple", "",
	"SpillID:EventID:TrackerID:StationID:x:y:z:px:py:pz:ex:ey:ez:epx:epy:epz");
  TNtuple* data_samples = new TNtuple("data_ntuple", "",
	"SpillID:EventID:TrackerID:StationID:x:y:z:px:py:pz:ex:ey:ez:epx:epy:epz");

  // Apply the cuts and write the information that passes them to the files
  std::string maus_version = globals["maus_version"];
  std::vector<float> array_tru(9), array_rec(16);
  double tofmom(0.), tkumom(0.), tofmom_err(0.), tkumom_err(0.),
	 tkdmom(0.), tkdmom_err(0.), scale(1.), meanmom(0.);
  size_t nplanes, tkuid, tkdid;
  MiceTrack *track_cut(NULL);
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
    for (size_t i = 0; i < (size_t)T->GetEntries(); ++i) {

      // Display the progress in %
      pbar.GetProgress(i, (size_t)T->GetEntries());

      // Update the track pointed to by *track_tru and *track_rec
      T->GetEntry(i);

      // If there is no reconstruction, apply the cuts on the MC truth
      track_cut = track_rec ? track_rec : track_tru;

      // Check all the upstream cuts
      for (const Cut& cut : ucuts)
	if ( cut != All )
	    upassed[cut] = CheckCut(track_cut, cut, tkuid);

      // If they the energy lost between TOF1 and TKU is not consistent with a muon seeing
      // the expected amount of material given the diffuser setting, cut out
      upassed[EnergyLoss] = false;
      if ( (upassed[TOF01] || globals["notofs"]) && upassed[Track] ) {

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
//	  meanmom = tkumom/pow(tkumom_err, 2) + tofmom/pow(tofmom_err, 2);
//        meanmom /= (1./pow(tkumom_err, 2) + 1./pow(tofmom_err, 2));
	meanmom = tkumom;

	// Scale all the TKU momenta accordingly
	scale = meanmom/tkumom;
	for (size_t planeid = 0; planeid < 5; planeid++) {
//	  track_cut->mom[planeid].SetX(scale*track_cut->mom[planeid].x());
//	  track_cut->mom[planeid].SetY(scale*track_cut->mom[planeid].y());
//	  track_cut->mom[planeid].SetZ(scale*track_cut->mom[planeid].z());

//	  track_cut->mome[planeid].SetX(scale*track_cut->mome[planeid].x());
//	  track_cut->mome[planeid].SetY(scale*track_cut->mome[planeid].y());
//	  track_cut->mome[planeid].SetZ(scale*track_cut->mome[planeid].z());
	}

	// Check if the new momentum is within acceptance
	if ( meanmom > (double)globals["tku_mom_min"] && meanmom < (double)globals["tku_mom_max"] )
	    upassed[Momentum] = true;

	// Check the energy loss cut
        if ( fabs(meanmom-tofmom) < Tolerance(globals["diffuser_thickness"]) )
      	    upassed[EnergyLoss] = true;
      }

      // If the TOFs are absent, pass all the related cuts automatically
      if ( globals["notofs"] )
        for (const Cut& cut : {TOF0SP, TOF1SP, TOF01, EnergyLoss})
	    upassed[cut] = true;

      // Make an diffuser aperture cut in the pure simulation, issue with the outer module (TODO)
/*      if ( track_tru && !globals["notofs"] )
        if ( track_tru->pos[41].Pt() > 130 || track_tru->pos[46].Pt() > 130 ) {
//        if ( track_tru->pos[41].Pt() > 110 || track_tru->pos[46].Pt() > 110 ) {
	  track_cut->tk_nt[0] = 0;
          for (const Cut& cut : {Track, Fiducial, Chi2, Momentum, EnergyLoss})
	      upassed[cut] = false;
        }*/

      // If any of the cuts failed, set "All" cut to false
      upassed[All] = true;
      for (const Cut& cut : ucuts)
	if ( cut != All && !upassed[cut] ) {
	  upassed[All] = false;
	  break;
	}

      // Increment the upstream counters
      for (const Cut& cut : ucuts)
	if ( upassed[cut] )
	    ucounters[type][cut]++;

      // Check the downstream cuts only if particles passes the upstream cuts
      if ( upassed[All] ) {
        for (const Cut& cut : dcuts)
	  if ( cut != All )
	    dpassed[cut] = CheckCut(track_cut, cut, tkdid);

        if ( dpassed[Track] ) {
	  // Correct the momentum by taking the weighted average with the TKU extrapolation
	  tkumom = DownstreamMomentum(track_cut->mom[tkuid].Mag(), globals["absorber"]);
	  tkumom_err = 2.;
	  tkdmom = track_cut->mom[tkdid].Mag();
	  tkdmom_err = track_cut->mome[tkdid].Mag();
//	    meanmom = tkumom/pow(tkumom_err, 2) + tkdmom/pow(tkdmom_err, 2);
//          meanmom /= (1./pow(tkumom_err, 2) + 1./pow(tkdmom_err, 2));
	  meanmom = tkdmom;

	  // Scale all the TKD momenta accordingly
	  scale = meanmom/tkdmom;
	  for (size_t planeid = 5; planeid < 10; planeid++) {
//	    track_cut->mom[planeid].SetX(scale*track_cut->mom[planeid].x());
//	    track_cut->mom[planeid].SetY(scale*track_cut->mom[planeid].y());
//	    track_cut->mom[planeid].SetZ(scale*track_cut->mom[planeid].z());

//	    track_cut->mome[planeid].SetX(scale*track_cut->mome[planeid].x());
//	    track_cut->mome[planeid].SetY(scale*track_cut->mome[planeid].y());
//	    track_cut->mome[planeid].SetZ(scale*track_cut->mome[planeid].z());
	  }

	  // Check if the corrected momentum is within acceptance
	  if ( meanmom > (double)globals["tkd_mom_min"] && meanmom < (double)globals["tkd_mom_max"] )
	      dpassed[Momentum] = true;
	}

        // If any of the cuts failed, set "All" cut to false
        dpassed[All] = true;
        for (const Cut& cut : dcuts)
	  if ( cut != All && !dpassed[cut] ) {
	    dpassed[All] = false;
	    break;
	  }

        // Increment the downstream counters
        for (const Cut& cut : dcuts)
	  if ( dpassed[cut] )
	      dcounters[type][cut]++;
      }

      // If all the upstream cuts are passed, fill the output file
      if ( upassed[All] ) {

        // If the truth is present, fill the correponding tree
	if ( track_tru ) {
	  // If the track does not make it to the upstream reference plane, proceed
	  if ( track_tru->pos.find(globals["tku_vid"]) == track_tru->pos.end() )
	      continue;

	  // Fill the truth information
          array_tru[0] = track_tru->spillid;
  	  array_tru[1] = track_tru->eventid;
 	  for (const std::pair<size_t, TVector3> el : track_tru->pos) {
	    array_tru[2] = el.first;
	    array_tru[3] = track_tru->pos.at(el.first).x();
	    array_tru[4] = track_tru->pos.at(el.first).y();
	    array_tru[5] = track_tru->pos.at(el.first).z();
	    array_tru[6] = track_tru->mom.at(el.first).x();
	    array_tru[7] = track_tru->mom.at(el.first).y();
	    array_tru[8] = track_tru->mom.at(el.first).z();

	    // In the space between the two trackers, apply default apertures
	    if ( el.first > (size_t)globals["tku_vid"] && el.first < (size_t)globals["tkd_vid"] &&
		!apertures.IsIn(array_tru[3], array_tru[4], array_tru[5]) && !dpassed[All] )
		continue;

	    // Fill the truth downstream if there is a valid TKD track only
	    if ( el.first < (size_t)globals["tkd_vid"] || dpassed[All] )
	        truth_samples->Fill(&array_tru[0]);

	    // Fill the uncut truth even if a TKD track is absent, but only if it is in the tracker
	    if ( el.first < (size_t)globals["tkd_vid"] || (track_tru->tk_maxr[1] < 150 &&
							   track_tru->pid.at(el.first) == -13) )
	        uncut_truth_samples->Fill(&array_tru[0]);
	  }
        }

	if ( track_rec ) {
          array_rec[0] = track_rec->spillid;
  	  array_rec[1] = track_rec->eventid;
	  nplanes = dpassed[All] ? 10 : 5;
	  for (size_t sid = 0; sid < nplanes; sid++) {
	    if ( track_rec->pos.find(sid) != track_rec->pos.end() ) {
	      array_rec[2] = (int)(sid/5);
	      array_rec[3] = (int)(sid%5);
	      array_rec[4] = track_rec->pos.at(sid).x();
	      array_rec[5] = track_rec->pos.at(sid).y();
	      array_rec[6] = track_rec->pos.at(sid).z();
	      array_rec[7] = track_rec->mom.at(sid).x();
	      array_rec[8] = track_rec->mom.at(sid).y();
	      array_rec[9] = track_rec->mom.at(sid).z();
	      array_rec[10] = track_rec->pose.at(sid).x();
	      array_rec[11] = track_rec->pose.at(sid).y();
	      array_rec[12] = track_rec->pose.at(sid).z();
	      array_rec[13] = track_rec->mome.at(sid).x();
	      array_rec[14] = track_rec->mome.at(sid).y();
	      array_rec[15] = track_rec->mome.at(sid).z();
	      if ( !track_tru ) {
	        data_samples->Fill(&array_rec[0]);
	      } else {
		if ( sid == 5 && (track_tru->pos.find((size_t)globals["tkd_vid"])) ==
				  track_tru->pos.end() )
		    continue;
	        recmc_samples->Fill(&array_rec[0]);
	      }
	    }
	  }
	}
      }

      // If requested, fill the diagnostic graphs
      if ( globals["diagnostics"] ) {
	if ( upassed[TOF1SP] && upassed[Track] && upassed[Chi2] &&
	     upassed[Fiducial] && upassed[Momentum] )
	    uhists[TOF0SP][type]->Fill(track_cut->tof_nsp[0]);
	if ( upassed[TOF0SP] && upassed[Track] && upassed[Chi2] &&
	     upassed[Fiducial] && upassed[Momentum] )
	    uhists[TOF1SP][type]->Fill(track_cut->tof_nsp[1]);
	if ( upassed[TOF0SP] && upassed[TOF1SP] && upassed[Track] && upassed[Chi2] &&
	     upassed[Fiducial] && upassed[Momentum] && upassed[EnergyLoss] )
	    uhists[TOF01][type]->Fill(track_cut->t[1]-track_cut->t[0]-1e9*tofdz/299792458);
	if ( upassed[TOF01] )
	    uhists[Track][type]->Fill(track_cut->tk_nt[0]);
	if ( upassed[TOF01] && upassed[Track] && upassed[Fiducial]
	     && upassed[Momentum] && upassed[EnergyLoss] )
	    uhists[Chi2][type]->Fill(track_cut->tk_chi2[0]);
	if ( upassed[TOF01] && upassed[Track] && upassed[Chi2]
	     && upassed[Momentum] && upassed[EnergyLoss] )
	    uhists[Fiducial][type]->Fill(track_cut->tk_maxr[0]);
	if ( upassed[TOF01] && upassed[Track] && upassed[Chi2]
	     && upassed[Fiducial] && upassed[EnergyLoss] )
	    uhists[Momentum][type]->Fill(track_cut->mom[tkuid].Mag());
	if ( upassed[TOF01] && upassed[Track] && upassed[Chi2]
	     && upassed[Fiducial] && upassed[Momentum] )
	    uhists[EnergyLoss][type]->Fill(track_cut->mom[tkuid].Mag()-tofmom);

	if ( upassed[All] ) {
	  dhists[Track][type]->Fill(track_cut->tk_nt[1]);
	  if ( dpassed[Track] && dpassed[Fiducial] && dpassed[Momentum] )
	      dhists[Chi2][type]->Fill(track_cut->tk_chi2[1]);
	  if ( dpassed[Track] && dpassed[Chi2] && dpassed[Momentum] )
	      dhists[Fiducial][type]->Fill(track_cut->tk_maxr[1]);
	  if ( dpassed[Track] && dpassed[Chi2] && dpassed[Fiducial] )
	      dhists[Momentum][type]->Fill(track_cut->mom[tkdid].Mag());
	}
      }
    }
  }

  // Fill the output file with the Ntuples and information about the sample
  types.resize(0);
  out->cd();
  if ( truth_samples->GetEntries() )
    truth_samples->Write("Truth");
  delete truth_samples;

  if ( uncut_truth_samples->GetEntries() )
      uncut_truth_samples->Write("UncutTruth");
  delete uncut_truth_samples;

  if ( recmc_samples->GetEntries() ) {
    recmc_samples->Write("RecMC");
    types.push_back("recmc");
  }
  delete recmc_samples;

  if ( data_samples->GetEntries() ) {
    data_samples->Write("Data");
    types.push_back("data");
  }
  delete data_samples;

  TNamed("MausVersion", maus_version.c_str()).Write();
  out->Close();

  // Write the counters to a file
  ofstream fcounters;
  fcounters.open("diag_counters.txt");
  for (const std::string type : types ) {
    fcounters << type << "\n";
    for (const Cut& cut : ucuts)
        fcounters << std::left << std::setw(30)
		  << "US "+CutName.at(cut) << "  " << ucounters[type][cut] << "\n";
    fcounters << "\n";
    for (const Cut& cut : dcuts)
        fcounters << std::left << std::setw(30)
		  << "DS "+CutName.at(cut) << "  " << dcounters[type][cut] << "\n";
    fcounters << "\n";
  }
  fcounters.close();

  // Move the counters to the appropriate directory
  Pitch::print(Pitch::info, "Moving the counters to "+run_name+"/diag");
    std::string sysCmd = "mkdir -p "+run_name+"/diag; "
	  "mv diag_counters.txt "+run_name+"/diag; done";
  if ( std::system(sysCmd.c_str()) )
      Pitch::print(Pitch::error, "Couldn't move the counters");

  // If requested, draw and save the diagnostic graphs
  if ( globals["diagnostics"] ) {

    // Remove the histograms that are empty
    for (const std::string type : {"truth", "recmc", "data"} )
      if ( std::find(types.begin(), types.end(), type) == types.end() ) {
        for (std::pair<Cut, std::map<std::string, TH1F*>> el : uhists)
	    el.second.erase(type);
        for (std::pair<Cut, std::map<std::string, TH1F*>> el : dhists)
  	    el.second.erase(type);
      }

    // Intialize the drawer and the info box
    TFile *outfile = new TFile("diag_histograms.root", "RECREATE");
    Beam::Drawer drawer;
    if ( globals["mice"] ) {
      std::string data_type = globals["data"] ? "Preliminary" : "[simulation]";
      std::string maus_version = globals["maus_version"];
      InfoBox* info = new InfoBox(data_type, maus_version,
    	globals["run_name"].AsString(), globals["user_cycle"].AsString());
      info->SetPosition("tl", .02, .24);
      drawer.SetInfoBox(info);
    }

    // Print the upstream cut diagnostic histograms
    for (std::pair<Cut, std::map<std::string, TH1F*>> el : uhists) {
      // Scale all the histograms to 1
      for (const std::string& type : types) {
	el.second[type]->Write();
	el.second[type]->Sumw2();
        el.second[type]->Scale(1./el.second[type]->GetEntries());
      }

      // Draw
      drawer.SaveStack(el.second, "diag_"+uhistdefs[el.first].name, false, "EMR");
    }

    // Print the downstream cut diagnostic histograms
    for (std::pair<Cut, std::map<std::string, TH1F*>> el : dhists) {
      // Scale all the histograms to 1
      for (const std::string& type : types) {
	el.second[type]->Write();
	el.second[type]->Sumw2();
        el.second[type]->Scale(1./el.second[type]->GetEntries());
      }

      // Draw
      drawer.SaveStack(el.second, "diag_"+dhistdefs[el.first].name, true, "EMR");
    }

    // Close the output file
    outfile->Close();

    // Move the plots to the appropriate directory
    Pitch::print(Pitch::info, "Moving the diagnostics graphs to "+run_name+"/diag");
    std::string sysCmd = "mkdir -p "+run_name+"/diag; "
	  "mv diag*.pdf diag*.root "+run_name+"/diag; done";
    if ( std::system(sysCmd.c_str()) )
        Pitch::print(Pitch::error, "Couldn't move the diagnostics plots");
  }
}
