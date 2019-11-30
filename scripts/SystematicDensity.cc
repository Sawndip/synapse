// Cpp includes
#include <iostream>
#include <vector>
#include <map>

// Additional modules
#include "Globals.hh"
#include "DPoisson.hh"
#include "Extractor.hh"

/** @file  SystematicErrors.cc
 *
 *  @brief Evaluates systematics errors in the MICE data.
 *
 *	   Algorithm that evaluates systematic errors to the MICE data provided with
 *	   with an array of MAUS simulation of the beam.
 **/


/** @brief	Returns the level at a given alpha by using linear interpolation
 *
 *  @param	levs		Sorted array of levels
 *  @param	trans		True transmission
 *  @param	alpha		Sorted array of levels
 */
double Level(const std::vector<double>& levs,
	     const double& trans,
	     const double alpha) {

  // Compute the ordering from alpha
  double ord = (levs.size()+1.)*alpha/trans - 1.;

  // If the ordering is below 1, return the maximum density
  if ( ord < 1. )
      return levs.front();

  // If the ordering is above the transmission, return 0
  if ( ord > levs.size()-1 )
      return 0.;

  // In any oder case, return the interpolated value
  size_t id = (int)ord;
  double xi = ord-id;
  return levs[id]*xi+levs[id]*(1.-xi);
}

/** @brief	Returns the correction graph evaluated every 0.01\,% (10k+1 points)
 *
 *  @param	tlevs		Sorted array of true levels
 *  @param	rlevs		Sorted array of reconstructed levels
 *  @param	ttrans		True transmission
 *  @param	rtrans		Reconstructed transmission
 */
TGraphErrors* CorrectionGraph(const std::vector<double>& tlevs,
			      const std::vector<double>& rlevs,
			      const double& ttrans,
			      const double& rtrans) {

  // Initialize the graph
  TGraphErrors* graph = new TGraphErrors();

  // Loop over alpha ranging from 0 to 1 by steps of 0.01\,\%
  // Interpolate the density of rec and the truth, compute ratio
  double alpha, tlev, rlev, ratio;
  for (size_t i = 0; i < 1e4+1; i++) {
    alpha = i*1e-4;
    tlev = Level(tlevs, ttrans, alpha);
    rlev = Level(rlevs, rtrans, alpha);
    ratio = 0.;
    if ( rlev > 0. )
        ratio = tlev/rlev;
    graph->SetPoint(i, alpha, ratio);
  }
  
  return graph;
}

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

  // Get the run name from the first file name (remove superfluous with regex)
  std::string run_name = std::regex_replace(globals.GetDataFiles()[0],
				std::regex("import_|(.*/201)|(/tk.*)"), std::string(""));
  run_name.insert(0, "201");
  Pitch::print(Pitch::info, "Run name: "+run_name);

  // Set up the requested plane IDs
  std::map<std::string, size_t> vids = {{"tku", globals["tku_vid"]}, {"tkd", globals["tkd_vid"]}};
  std::map<std::string, size_t> stids = {{"tku", 0}, {"tkd", 5}};
  std::vector<std::string> types = {"utruth", "recmc"};
  std::vector<std::string> default_dets = {"tku", "tkd"};
  std::vector<std::string> dets;
  std::map<std::string, std::vector<size_t>> plane_ids;
  std::map<std::string, std::vector<size_t>> default_plane_ids = 
	{{"utruth",{vids["tku"], vids["tkd"]}}, {"recmc",{0, 5}}};

  // Extract the types of systematics
  size_t nsets = globals.GetDataFiles().size();
  std::vector<std::string> set_names(nsets);
  std::string set_name;
  for (size_t i = 0; i < nsets; i++) {
    set_name = globals.GetDataFiles()[i];
    set_name = std::regex_replace(set_name, std::regex("(.*/tk)"), std::string(""));
    set_name = std::regex_replace(set_name, std::regex("(/.*)"), std::string(""));
    set_names[i] = "tk"+set_name;
  }

  // Loop over the simulation files, evaluate and store the corrections for each
  // amplitude bin upstream and downstream of the absorber
  TFile *outfile = new TFile(TString::Format("%s_density_syst.root", run_name.c_str()), "RECREATE");
  std::vector<std::map<std::string, TGraphErrors*>> corrections(nsets);
  std::map<std::string, std::map<std::string, std::vector<double>>> levels;
  std::map<std::string, std::map<std::string, double>> trans;
  std::string set_det;
  for (size_t iFile = 0; iFile < nsets; iFile++) {
    // Set up the extractor
    Beam::Extractor ext({globals.GetDataFiles()[iFile]});

    // Get the setting, determine which stations are required
    set_name = set_names[iFile];
    set_det = set_name.substr(0, 3);
    set_name = set_name.substr(4);
    Pitch::print(Pitch::info, "Importing the "+set_name+" systematics set for "+set_det);

    plane_ids = default_plane_ids;
    if ( set_name != "base" ) {
      dets = {set_det};
    } else {
      set_det = "";
      dets = default_dets;
    }

    // Extract the upstream and downstream stations only
    std::map<std::string, Beam::Stream> streams = ext.GetStreams(types, plane_ids, globals["maxn"]);

    // Set the density levels
    for (const std::string& type : types) {
      for (size_t i = 0; i < plane_ids[type].size(); i++)
          trans[type][default_dets[i]] = streams[type].Transmission(plane_ids[type][i]);

      for (size_t i = 0; i < dets.size(); i++) {
        streams[type][plane_ids[type][i]].SetDensityLevels();
        levels[type][dets[i]] = streams[type][plane_ids[type][i]].DensityLevels();
        std::sort(levels[type][dets[i]].rbegin(), levels[type][dets[i]].rend());
      }
    }

    // Compute the correction vectors
    for (const std::string det : dets) {
      // Fill correction graph, save it to file
      corrections[iFile][det] = CorrectionGraph(levels["utruth"][det], levels["recmc"][det],
						trans["utruth"][det], trans["recmc"][det]);
      for (size_t i = 5000; i < 5010; i++)
	  std::cerr << corrections[iFile][det]->GetY()[i] << std::endl;
      corrections[iFile][det]->SetName(TString::Format(
					"syst_corr_%s_%s", det.c_str(), set_name.c_str()));

      // If the run corresponds to the baseline, draw the histograms
      if ( set_name == "base" ) {
	//Correction histogram
        TCanvas* c = new TCanvas("c", "c", 1200, 800);
        corrections[iFile][det]->Draw("AL");
        corrections[iFile][det]->SetMinimum(0.8);
        corrections[iFile][det]->SetMaximum(1.2);
	corrections[iFile][det]->GetXaxis()->SetRangeUser(0, 1);

        TLine *line = new TLine(0, 1, 1, 1);
        line->SetLineColorAlpha(kBlack, .5);
        line->Draw("SAME");

        c->SaveAs(TString::Format("syst_density_corr_%s.pdf", det.c_str()));
        delete c;
      }

      // Save graphs to file
      outfile->cd();
      corrections[iFile][det]->Write(corrections[iFile][det]->GetName());
    }
  } // End of the list of files

  outfile->Close();

  Pitch::print(Pitch::info, "Moving the systematics graphs to "+run_name+"/syst");
  std::string sysCmd = "mkdir -p "+run_name+"/syst; mv syst*.pdf *_syst.root "+run_name+"/syst; done";
  if ( std::system(sysCmd.c_str()) )
      Pitch::print(Pitch::error, "Couldn't move systematics plots");
}
