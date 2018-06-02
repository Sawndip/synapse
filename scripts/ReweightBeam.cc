// Cpp includes
#include <iostream>
#include <algorithm>
#include <regex>
#include <map>

// Root includes
#include "TStyle.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TRandom3.h"

// Additional modules
#include "Globals.hh"
#include "ProgressBar.hh"
#include "InfoBox.hh"
#include "Bunch.hh"

/** @file  ReweightBeam.cc
 *
 *  @brief Reweights the MAUS simulation to reproduce the amplitude distribution observed in data.
 *
 *	   Algorithm that reweights the MAUS simulation to faithfully reproduce the amplitude
 *	   distribution observed in hte corresponding data
 *   	   downstream of the MICE absorber. It uses the requested density estimator to reconstruct
 *	   the minimum local density as a function of the fraction of the beam included.
 **/

/** @brief Reweighting method */
enum Method {
  maxdist, 	///< Fixes everything to the largest discrepancy
  binbybin	///< Fixes bin by bin from largest to smallest
};

/** @brief Reweighting metric */
enum Metric {
  own, 		///< Simulation metric
  data		///< Data metric
};

/** @brief Computes the amplitudes of particles in the metric defined by covmat
 *
 *  @param	beam		Particle beam object
 *  @param	covmat		Custom metric
 *
 *  @return			Vector of amplitudes
 */
std::vector<double> CustomAmplitudes(const Beam::Bunch& beam,
				     const Matrix<double>& covmat) {

  std::vector<double> amps(beam.Size());
  double eps = beam.NormEmittance().GetValue();
  Matrix<double> x(4, 1), xt(1, 4);
  for (size_t i = 0; i < amps.size(); i++) {
    x[0][0] = beam.Samples("px")[i];
    x[1][0] = beam.Samples("py")[i];
    x[2][0] = beam.Samples("x")[i];
    x[3][0] = beam.Samples("y")[i];
    xt = x.Transpose();
    amps[i] = eps*(xt*covmat*x)[0][0];
  }

  return amps;
}

/** @brief	Main function
 *
 *  @param	argc		Number of command line arguments
 *  @param	argv		List of command line aguments
 */
int main(int argc, char ** argv) {

  // Reroot the info messages to the terminal for now
  Pitch::setAnOutput(Pitch::info, std::cerr);

  // Reweighting algorithm global parameters (stored in globals)
  Globals &globals = Globals::GetInstance(argc, argv);

  // Get the run name from the first file name (remove superfluous with regex)
  std::string run_name = std::regex_replace(globals.GetDataFiles()[0],
				std::regex("import_|(.*/)|(.root)"), std::string(""));
  Pitch::print(Pitch::info, "Run name: "+run_name);

  // Set the algorithm and the metric
  Method method = (Method)(int)globals["reweight_alg"];
  Metric metric = (Metric)(int)globals["reweight_metric"];

  // Do not draw the default stat box
  gStyle->SetOptStat(0);
  gErrorIgnoreLevel = kWarning;	// Silence chi^2 test

  // Containers for the reconstructed MC samples and real data, one per tracker station
  std::vector<std::string> vars = {"x", "y", "px", "py", "pz"}; // Variables
  std::map<std::string, std::map<std::string, std::vector<double>>> rsamples;
  std::vector<size_t> inids, rewids; 	// Unique recon MC event ids
  size_t factor = 1e6;			// Factor to apply to the spill id
  std::string maus_version, buff;

  // Fill the sample maps from the imported ROOT file
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
      if ( !buff.size() ) {
	buff == maus_version;
      } else {
	if ( maus_version != buff )
	    Pitch::print(Pitch::warning, "The MAUS versions do not match, should not proceed");
      }
    }

    // Import the reconstructed MC/data if it is provided
    for (const std::string& type_uc : {"RecMC", "Data"}) {
      std::string type = type_uc;
      std::transform(type.begin(), type.end(), type.begin(), ::tolower);
      if ( data_file.GetListOfKeys()->Contains(type_uc.c_str()) ) {
        TNtuple* rec_samples = (TNtuple*)data_file.Get(type_uc.c_str());
        float* ntuple;
	std::string det;
        size_t st, spill, evt;

        Pitch::print(Pitch::info, "Processing the "+type_uc);
        ProgressBar pbar;
        for (size_t i = 0; i < (size_t)rec_samples->GetEntries(); ++i) {

          // Display the progress in %
          pbar.GetProgress(i, (size_t)rec_samples->GetEntries());

	  // Fetch the variables from the Ntuple
	  rec_samples->GetEntry(i);
	  ntuple = rec_samples->GetArgs();
          spill = ntuple[0];
          evt = ntuple[1];
	  det = ntuple[2] ? "tkd" : "tku";
	  st = ntuple[3];

	  // Only keep the reference plane information in TKU
	  if ( det != "tku" || st != 0 )
	      continue;

	  // Store the unique ids of the events
	  if ( type == "recmc" )
	      inids.push_back(factor*spill+evt);

	  // Fill the reconstructed samples
          rsamples[type]["x"].push_back(ntuple[4]);
          rsamples[type]["y"].push_back(ntuple[5]);
          rsamples[type]["px"].push_back(ntuple[7]);
          rsamples[type]["py"].push_back(ntuple[8]);
          rsamples[type]["pz"].push_back(ntuple[9]);
        }
      }
    }

  } // End of the list of files

  // Compute the amplitudes of the data and recmc. Data will not change
  Beam::Bunch beam_data, beam_recmc;
  std::vector<double> amps_data, amps_recmc;
  beam_data = Beam::Bunch(rsamples["data"], 0., "data");
  if ( globals["corrected"] )
      beam_data.SetCorrectedAmplitudes();
  amps_data = beam_data.Amplitudes();
  Matrix<double> covmat = beam_data.CovarianceMatrix().Inverse();

  beam_recmc = Beam::Bunch(rsamples["recmc"], 0., "recmc");
  if ( metric == own ) {
    if ( globals["corrected"] )
        beam_recmc.SetCorrectedAmplitudes();
    amps_recmc = beam_recmc.Amplitudes();
  } else if ( metric == data ) {
    amps_recmc = CustomAmplitudes(beam_recmc, covmat);
  }

  // Extract the ellipses
  TEllipse *ell_data = beam_data.Ellipse("x", "px");
  ell_data->SetLineStyle(2);
  ell_data->SetLineColor(kGreen+2);

  TEllipse *ell_before = beam_recmc.Ellipse("x", "px");

  // Draw the amplitude scatter plot before reweighting 
  TCanvas* cscat = new TCanvas("c", "c", 1200, 800);
  beam_recmc.AmplitudeScatter("x", "px")->Draw();
  cscat->SaveAs("reweight_scat_before.pdf");
  delete cscat;

  // Choose the binning depending on the input emittance
  int epsset;
  double meanamp = Math::Mean(amps_data);
  if ( meanamp < 15 ) {
    epsset = 3;
  } else if ( meanamp < 30 ) {
    epsset = 6;
  } else {
    epsset = 10;
  }
  size_t nbins = 100;
  if ( epsset == 6 )
    nbins = 70;
  if ( epsset == 10 )
    nbins = 50;

  // Initialize and fill data and a recmc histogram, scale them to 1
  double maxamp = 200;
  TH1F* hamps_data = new TH1F("amps_data", ";A_{#perp}  [mm];PDF", nbins, 0, maxamp);
  hamps_data->FillN(amps_data.size(), &amps_data[0], NULL);
  hamps_data->Sumw2();
  hamps_data->Scale(1./hamps_data->GetEntries());

  TH1F* hamps_recmc = new TH1F("amps_recmc", ";A_{#perp}  [mm];PDF", nbins, 0, maxamp);
  hamps_recmc->FillN(amps_recmc.size(), &amps_recmc[0], NULL);
  hamps_recmc->Sumw2();
  hamps_recmc->Scale(1./hamps_recmc->GetEntries());

  // Initialize the drawer and the stack
  Beam::Drawer drawer;
  std::map<std::string, TH1F*> hists;
  hists["data"] = hamps_data;
  hists["recmc"] = hamps_recmc;

  // Initialize the info box
  InfoBox *info = new InfoBox("Preliminary", maus_version, globals["run_name"], globals["user_cycle"]);
  drawer.SetInfoBox(info);

  // Find the ids of the points to remove from the recmc sample in order for it to agree with data
  TH1F* hamps_ratio;
  TRandom3 rdmzer;
  double max(0), maxratio(1), ratio, ran, kstest(0), chi2test(0);
  size_t binmax, binid, ite(0);
  std::vector<size_t> ids;
  Pitch::print(Pitch::info, "Reweighting the simulation");
  while ( chi2test < .9 || kstest < .9 ) {

    Pitch::print(Pitch::debug, "Reweighting iteration "+std::to_string(ite));

    // Draw the input histograms
    if ( !ite )
        drawer.SaveStack(hists, "reweight_amps_before", true);

    // Divide the data sample by the recmc sample to get probability ratios
    hamps_ratio = (TH1F*)hamps_data->Clone();
    hamps_ratio->Divide(hamps_recmc);

    if ( method == binbybin ) { // Scale the from highest bin to lowest

      // Loop over the bins from the overflow to bin 1
      for (size_t i = 0; i < (size_t)hamps_ratio->GetNbinsX()+1; i++) {
        ratio = hamps_ratio->GetBinContent((size_t)hamps_ratio->GetNbinsX()+1-i);

        // Loop over the amplitudes and randomly reject some of them to fit the true PDF (data)
	// if they fall inside the bin currently under consideration
        // Retain the ids of the particles that make it through the cut
        ids.resize(0);
        for (size_t j = 0; j < amps_recmc.size(); j++) {
          // Find the bin in which the amplitude lies 
          binid = hamps_ratio->GetXaxis()->FindBin(amps_recmc[j]);
	  if ( binid == (size_t)hamps_ratio->GetNbinsX()+1-i ) {
    
            // Randomly sample a number between 0 and 1
            ran = rdmzer.Uniform(0, 1);

            // If that point is smaller than the ratio in the bin, let it through, remove otherwise
            if ( ran < hamps_ratio->GetBinContent(binid) )
	        ids.push_back(j);

	  } else {
	    ids.push_back(j);
	  }
        }

        // Remove the particles that correspond to the amplitudes to ditch
        rsamples["rew"].clear();
        rewids.resize(ids.size());
        for (size_t j = 0; j < ids.size(); j++) {
          for (const std::string& var : vars )
	      rsamples["rew"][var].push_back(rsamples["recmc"][var][ids[j]]);

          rewids[j] = inids[ids[j]];
        }
        rsamples["recmc"] = rsamples["rew"];
        inids = rewids;

        // Update the amplitudes
        beam_recmc = Beam::Bunch(rsamples["recmc"]);
        if ( metric == own ) {
          if ( globals["corrected"] )
              beam_recmc.SetCorrectedAmplitudes();
          amps_recmc = beam_recmc.Amplitudes();
        } else if ( metric == data ) {
          amps_recmc = CustomAmplitudes(beam_recmc, covmat);
        }

        // Fill the histograms, scale them to each other
        hamps_recmc->Reset();
        hamps_recmc->FillN(amps_recmc.size(), &amps_recmc[0], NULL);
        hamps_recmc->Scale(1./hamps_recmc->GetEntries());
        
	drawer.SaveStack(hists,
	  "reweight_amps_ite_"+std::to_string((int)ite)+"_bin"+std::to_string((int)i), true);

        // Update the probability ratios
        hamps_ratio = (TH1F*)hamps_data->Clone();
        hamps_ratio->Divide(hamps_recmc);
      }

    } else if ( method == maxdist ) { // Scale everything at once to the most distant bin

      // Find the relatively most distant bin that is higher than the MC (quadratic)
      max = 0;
      for (size_t i = 0; i < (size_t)hamps_ratio->GetNbinsX(); i++) {
        ratio = hamps_ratio->GetBinContent(i+1);
        if ( ratio > 1 )
          if ( pow((ratio-1.)/hamps_ratio->GetBinError(i+1), 2) > max ) {
	    max = pow((ratio-1.)/hamps_ratio->GetBinError(i+1), 2);
  	    maxratio = ratio;
 	    binmax = i+1;
          }
      }
      Pitch::print(Pitch::debug, "Most distant bin: "+std::to_string(binmax));

      // Scale it so that the most distant bin has a ratio of one (keep it as is)
      hamps_ratio->Scale(1./maxratio);

      // Loop over the amplitudes and randomly reject some of them to fit the true PDF (data)
      // Retain the ids of the particles that make it through the cut
      ids.resize(0);
      for (size_t i = 0; i < amps_recmc.size(); i++) {
        // Find the bin in which the amplitude lies 
        binid = hamps_ratio->GetXaxis()->FindBin(amps_recmc[i]);
    
        // Randomly sample a number between 0 and 1
        ran = rdmzer.Uniform(-3, 1);

        // If that point is smaller than the ratio in the bin, let it through, remove otherwise
        if ( ran < hamps_ratio->GetBinContent(binid) )
	    ids.push_back(i);
      }

      // Remove the particles that correspond to the amplitudes to ditch
      rsamples["rew"].clear();
      rewids.resize(ids.size());
      for (size_t i = 0; i < ids.size(); i++) {
        for (const std::string& var : vars )
	    rsamples["rew"][var].push_back(rsamples["recmc"][var][ids[i]]);

        rewids[i] = inids[ids[i]];
      }
      rsamples["recmc"] = rsamples["rew"];
      inids = rewids;
    }

    // Update the amplitudes
    beam_recmc = Beam::Bunch(rsamples["recmc"]);
    if ( metric == own ) {
      if ( globals["corrected"] )
          beam_recmc.SetCorrectedAmplitudes();
      amps_recmc = beam_recmc.Amplitudes();
    } else if ( metric == data ) {
      amps_recmc = CustomAmplitudes(beam_recmc, covmat);
    }

    // Fill the histograms, scale them to each other and draw
    hamps_recmc->Reset();
    hamps_recmc->FillN(amps_recmc.size(), &amps_recmc[0], NULL);
    hamps_recmc->Scale(1./hamps_recmc->GetEntries());

    drawer.SaveStack(hists, "reweight_amps_ite_"+std::to_string((int)ite), true);

    // Run the KS and chi^2 tests
    kstest = hamps_data->KolmogorovTest(hamps_recmc);
    Pitch::print(Pitch::debug, "Kolmogov-Smirnov test: "+std::to_string(kstest));

    chi2test = hamps_data->Chi2Test(hamps_recmc, "UUNORM");
    Pitch::print(Pitch::debug, "Chi2 test: "+std::to_string(chi2test));

    // Increment
    ite++;
  }

  // Draw the exit comparison
  drawer.SaveStack(hists, "reweight_amps_after", true);

  // Draw the amplitude scatter plot after reweighting 
  cscat = new TCanvas("c", "c", 1200, 800);
  beam_recmc.AmplitudeScatter("x", "px")->Draw();
  cscat->SaveAs("reweight_scat_after.pdf");
  delete cscat;

  // Draw the beam RMS ellipses of the data and the RecMC before and after reweighting
  TEllipse* ell_after = beam_recmc.Ellipse("x", "px");
  ell_after->SetLineColor(kBlue);

  TCanvas* cell = new TCanvas("c", "c", 1200, 800);
  TH2F *frame = new TH2F("ell_frame", ";x [mm]; p_{x}  [MeV/c]", 1, -100, 100, 1, -50, 50);
  frame->Draw();
  ell_before->Draw("SAME");
  ell_after->Draw("SAME");
  ell_data->Draw("SAME");
  cell->SaveAs("reweight_ellipses.pdf");
  delete cell;

  // Sort the recon event ids to optimize the speed of the binary search
  std::sort(rewids.begin(), rewids.end());

  // Move the diagnostics histograms
  Pitch::print(Pitch::info, "Moving the reweighting graphs to "+run_name+"_reweighted/rew");
  std::string sysCmd = "mkdir -p "+run_name+"_reweighted/rew; mv reweight*.pdf "
							+run_name+"_reweighted/rew; done";
  if ( std::system(sysCmd.c_str()) )
      Pitch::print(Pitch::error, "Couldn't move reweighting plots");

  // Now update the file with the reweighted sample
  size_t id, lastid(0);
  bool lastfound(false);
  std::string outname = "import_"+run_name+"_reweighted.root";
  TFile *out = new TFile(outname.c_str(), "RECREATE");

  std::map<std::string, TNtuple*> samples;
  samples["truth"] = new TNtuple("truth_ntuple", "",
	"SpillID:EventID:VirtualPlaneID:x:y:z:px:py:pz");
  samples["uncuttruth"] = new TNtuple("uncut_truth_ntuple", "",
	"SpillID:EventID:VirtualPlaneID:x:y:z:px:py:pz");
  samples["recmc"] = new TNtuple("recmc_ntuple", "",
	"SpillID:EventID:TrackerID:StationID:x:y:z:px:py:pz:ex:ey:ez:epx:epy:epz");
  samples["data"] = new TNtuple("data_ntuple", "",
	"SpillID:EventID:TrackerID:StationID:x:y:z:px:py:pz:ex:ey:ez:epx:epy:epz");

  std::vector<size_t> tids;
  Pitch::print(Pitch::info, "Writing the reweighted file");
  for (const std::string& file : globals.GetDataFiles()) {
    // Set up the ROOT file and data pointer
    TFile data_file(file.c_str());		// Load the MAUS output file
    if ( !data_file.IsOpen() ) {
      Pitch::print(Pitch::error, "Data file not found: "+file);
      return 2;
    }
    Pitch::print(Pitch::info, "Processing "+file);

    // Import the reconstructed MC/data if it is provided
    for (const std::string& type_uc : {"Truth", "UncutTruth", "RecMC", "Data"}) {
      std::string type = type_uc;
      std::transform(type.begin(), type.end(), type.begin(), ::tolower);
      std::string dictype = (type == "data") ? "data" : "recmc";
      if ( data_file.GetListOfKeys()->Contains(type_uc.c_str()) ) {
        TNtuple* iSize = (TNtuple*)data_file.Get(type_uc.c_str());
        float* ntuple;
        size_t spill, evt;

        Pitch::print(Pitch::info, "Writing the "+type_uc);
        ProgressBar pbar;
        for (size_t i = 0; i < (size_t)iSize->GetEntries(); ++i) {

          // Display the progress in %
          pbar.GetProgress(i, (size_t)iSize->GetEntries());

	  // Fetch the variables from the Ntuple
	  iSize->GetEntry(i);
	  ntuple = iSize->GetArgs();

	  // If it is part of the accepted sample, write it to the output tree
	  if ( type != "data" ) {
            spill = ntuple[0];
            evt = ntuple[1];
	    id = factor*spill+evt;

	    if ( id == lastid ) { // If unique ID has already been searched for, rely on last search
	      if ( lastfound ) {
	          samples[type]->Fill(ntuple);
	      }

	    } else if ( std::binary_search(rewids.begin(), rewids.end(), id) ) {
	      samples[type]->Fill(ntuple);
	      lastfound = true;

	    } else {
	      lastfound = false;
	    }
	    lastid = id;

	  } else {
	    samples[type]->Fill(ntuple);
	  }

        }
      }
    }

  } // End of the list of files

  // Fill the output file with the Ntuples and information about the sample
  out->cd();
  if ( samples["truth"]->GetEntries() )
      samples["truth"]->Write("Truth");
  delete samples["truth"];

  if ( samples["uncuttruth"]->GetEntries() )
      samples["uncuttruth"]->Write("UncutTruth");
  delete samples["uncuttruth"];

  if ( samples["recmc"]->GetEntries() )
      samples["recmc"]->Write("RecMC");
  delete samples["recmc"];

  if ( samples["data"]->GetEntries() )
      samples["data"]->Write("Data");
  delete samples["data"];

  TNamed("MausVersion", maus_version.c_str()).Write();
  out->Close();
}
