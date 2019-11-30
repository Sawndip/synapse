// Cpp includes
#include <iostream>
#include <vector>
#include <map>

// Additional modules
#include "Globals.hh"
#include "Extractor.hh"

#include "Generator.hh"
#include "Scattering.hh"
#include "EnergyLoss.hh"
#include "Transport.hh"

#include "TPaletteAxis.h"

/** @file  DensityProfile.cc
 *
 *  @brief Produces non-parametric density profiles of the beam
 *
 *  	   Algorithm that produces non-parametric density profiles of the beam upstream and
 *   	   downstream of the MICE absorber. It uses the requested density estimator to reconstruct
 *	   the minimum local density as a function of the fraction of the beam included.
 **/


/** @brief	Main function
 *
 *  @param	levels		Array of density levels
 *  @param	alpha		Fraction of the beam contained in the contour
 *  @param	trans		Transmission
 */
double FindLevel(const std::vector<double>& levels, const double alpha, const double trans) {

  // Return if trivial
  double minalpha = trans/levels.size();
  double maxalpha = trans*levels.size()/(levels.size()+1);
  if ( alpha <= minalpha )
      return levels.front();
  if ( alpha >= maxalpha )
      return 0.;

  // Find the point just before alpha, interpolate
  double pos = (levels.size()+1.)*alpha/trans-1.;
  size_t id = (size_t)pos;
  double xi = pos - id;
  return xi*levels[id+1]+(1.-xi)*levels[id];
}

/** @brief	Main function
 *
 *  @param	argc		Number of command line arguments
 *  @param	argv		List of command line aguments
 */
int main(int argc, char ** argv) {

  // Amplitude algorithm global parameters (stored in globals)
  Globals &globals = Globals::GetInstance(argc, argv);

  // Types and plane ids of the data to extract. Empty means import all (recmc and data)
  size_t inid((size_t)globals["tku_vid"]), outid((size_t)globals["tkd_vid"]);
  std::vector<std::string> tvars = {"x", "y", "px", "py"};
  std::vector<std::string> types;
  std::map<std::string, std::vector<size_t>> plane_ids;
  for (const std::string type : {"truth", "recmc", "data"})
    if ( globals[type] ) {
      types.push_back(type);
      plane_ids[type] = (type == "truth") ? std::vector<size_t>({inid, outid}) :
					    std::vector<size_t>({0, 5});
    }

  // Get the requested streams from the beam extractor
  Beam::Extractor ext(globals.GetDataFiles());
  std::string run_name = ext.GetRunName();
  std::map<std::string, Beam::Stream> streams = ext.GetStreams(types, plane_ids, globals["maxn"]);

  // Evaluate the density profiles in the upstream and downstream reference planes
  TFile *outfile = new TFile("profiles_de.root", "RECREATE");
  std::map<std::string, TGraphErrors*> profiles;
  TGraphErrors* gratio;
  std::vector<double> levels;
  std::map<std::string, int> colors = {{"in", kBlack}, {"out", kBlue}};
  double alpha, err, ratio, trans, level;
  size_t planeid, size;
  for (const std::string& type : types) {

    // Initialize the multigraph
    TMultiGraph* mg = new TMultiGraph();
    mg->SetTitle(";Fraction #alpha;#rho_{#alpha}  [mm^{-2}(MeV/c)^{-2}]");

    // Loop over the input and output beams
    for (const std::string loc : {"in", "out"}) {  

      // Get the transmission
      Pitch::print(Pitch::info, "Profiling "+loc+" DE for "+type);
      planeid = (loc == "in") ? plane_ids[type][0] : plane_ids[type][1];
      streams[type][planeid].SetName(loc+"_"+type);
      size = streams[type][planeid].Size();
      trans = (double)size/streams[type].Front().Size();

      // Evaluate the density in all of the training points (necessary to get ordering)
      streams[type][planeid].SetDensityLevels(trans);
      levels = streams[type][planeid].DensityLevels();

      // Sort them by decreasing order of density
      std::sort(levels.rbegin(), levels.rend());

      // Fill the graph with the levels as a function of the sample fraction (1e4 points)
      profiles[loc] = new TGraphErrors();
      profiles[loc]->SetName(TString::Format("de_profile_%s_%s", type.c_str(), loc.c_str()));
      for (size_t i = 0; i < 1e4; i++) {
	alpha = (i+1.)*1e-4;
	level = FindLevel(levels, alpha, trans);
	err = level ? pow(4/(alpha*(1-alpha)), 1./4)/sqrt(size) : 0.;

	profiles[loc]->SetPoint(i, alpha, level);
	profiles[loc]->SetPointError(i, 0., err*level);
      }

      // Add them to the multigraph to draw, write them to the file to be saved
      mg->Add(profiles[loc], "LE3");
      profiles[loc]->SetLineWidth(2);
      profiles[loc]->SetFillStyle(1000);
      profiles[loc]->SetLineColor(colors[loc]);
      profiles[loc]->SetFillColorAlpha(colors[loc], .25);
      profiles[loc]->Write(profiles[loc]->GetName());
    }

    // Draw
    TCanvas *c = new TCanvas("c", "c", 1200, 800);
    mg->Draw("A");
    mg->SetMinimum(0);
    mg->GetXaxis()->SetLimits(0., 1.);
    c->SaveAs(TString::Format("de_profile_%s.pdf", type.c_str()));
    delete c;

    // Draw the ratio of the output over the input
    gratio = new TGraphErrors();
    gratio->SetLineColor(kBlue);
    gratio->SetLineWidth(2);
    gratio->SetFillStyle(1000);
    gratio->SetFillColorAlpha(kBlue+2, .25);
    gratio->SetName(TString::Format("de_ratio_%s", type.c_str()));
    gratio->SetTitle(";Fraction #alpha;#rho_{#alpha}^{u}/#rho_{#alpha}^{d}");
    for (size_t i = 0; i < (size_t)profiles["in"]->GetN(); i++) {
      if ( !profiles["out"]->GetY()[i] )
	  continue;

      alpha = profiles["in"]->GetX()[i];
      ratio = profiles["out"]->GetY()[i]/profiles["in"]->GetY()[i];
      gratio->SetPoint(i, alpha, ratio);
//      err = ratio*sqrt(pow(profiles["in"]->GetEY()[i]/profiles["in"]->GetY()[i], 2) +
//			pow(profiles["out"]->GetEY()[i]/profiles["out"]->GetY()[i], 2));
      err = ratio*profiles["out"]->GetEY()[i]/profiles["out"]->GetY()[i];
      gratio->SetPointError(i, 0., err);
    }

    c = new TCanvas("c", "c", 1200, 800);
    gratio->Draw("ALE3");
    gratio->Write(gratio->GetName());
    gratio->SetMinimum(0);
    gratio->SetMaximum(1.5);
    gratio->GetXaxis()->SetLimits(0., 1.);
    TLine *baseline = new TLine(0, 1., 1., 1.);
    baseline->SetLineColorAlpha(kBlack, .6);
    baseline->SetLineWidth(2);
    baseline->Draw("SAME");
    c->SaveAs(TString::Format("de_ratio_%s.pdf", type.c_str()));
    delete c;
  }
  outfile->Close();

  // Produce Poincaré sections of the estimator, if requested
  if ( globals["poincare"]  ) {
    Pitch::print(Pitch::info, "Saving the Poincaré sections");
    TFile* outfile = new TFile("poincare_de.root", "RECREATE");
    for (const std::string& type : types) {
      for (const std::string loc : {"in", "out"}) {
        // Get the sections and save them to a file, move them to a directory
        planeid = (loc == "in") ? plane_ids[type][0] : plane_ids[type][1];
	for (size_t i = 0; i < 4; i++) {
	  for (size_t j = i+1; j < 4; j++) {
            Pitch::print(Pitch::info, "Saving "+loc+" "+tvars[i]+tvars[j]+" for "+type);
	    TH2F* graph = streams[type][planeid].DensityEstimate().Graph2D(
						-100, 100, -100, 100, i, j, {0, 0, 0, 0});
	    graph->Write(graph->GetName());
	  }
	}
      }  
    }
    outfile->Close();
  }

  // Move the profiles to the appropriate directory
  Pitch::print(Pitch::info, "Moving the density graphs to "+run_name+"/de");
  std::string sysCmd = "mkdir -p "+run_name+"/de; mv de_*.pdf *_de.root "+run_name+"/de; done";
  if ( std::system(sysCmd.c_str()) )
      Pitch::print(Pitch::error, "Couldn't move density plots");
}
