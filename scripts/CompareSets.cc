// C++ includes
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <algorithm>

// Root includes
#include "TSystem.h"
#include "TFile.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TLine.h"
#include "TCanvas.h"

// Other includes
#include "InfoBox.hh"
#include "Pitch.hh"
#include "Scattering.hh"
#include "EnergyLoss.hh"
#include "Globals.hh"

/** @file  CompareSets.cc
 *
 *  @brief Compares the phase space density evolution of different magnetic channel settings.
 *
 *	   Algorithm that imports the phase space density evolution from several settings
 *	   stored in individual ROOT histogram files and combines them to produce summary
 *	   graphs compared to the theoretical expected emittance change.
 */

/** @brief	Main function
 *
 *  @param	argc		Number of command line arguments
 *  @param	argv		List of command line aguments
 */
int main(int argc, char ** argv) {

  // Reroot the info messages to the terminal for now
  Pitch::setAnOutput(Pitch::info, std::cerr);

  // Compare algorithm global parameters (stored in globals)
  Globals &globals = Globals::GetInstance(argc, argv);
  if ( !globals.GetDataFiles().size() ) {
    Pitch::print(Pitch::error, "No data file specified, please execute in the following fashion:");
    Pitch::print(Pitch::error, "./compare_sets [options] setting0.root [... settingN.root]");
    return 2;
  } else {
    for (const std::string& file : globals.GetDataFiles())
      if ( std::string(file).find(".root") > std::string(file).size() ) {
        Pitch::print(Pitch::error, "Invalid data file type: "+std::string(file));
	return 2;
      }
  }

  // Load libtree to get rid of the warnings...
  gSystem->Load("libTree");

  // Style settings
  std::vector<std::string> settings;
//  std::vector<std::string> vars = {"acut", "subeps", "vol", "de"};
  std::vector<std::string> vars = {"subeps", "vol"};
  std::map<std::string, int> markers = {{"acut",24}, {"subeps",25}, {"vol",26}, {"de",27}};
  std::vector<int> colors = {kBlack, kMagenta+1, kBlue+1, kCyan+1, kGreen+1, kOrange+1, kRed+1};

  // Initialize the graphs
  std::map<std::string, std::map<std::string, std::map<std::string, TGraphErrors*>>> graphs;
  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("Reconstructed emittance change;#epsilon_{in}  [mm];#Delta#epsilon/#epsilon_{in}  [%]");

  // Loop over the input ROOT files and fill the graphs
  double z, in, out, inerr, outerr, val, err;
  int id(-1);
  TLegend *setting_leg = new TLegend(.5, .625, .725, .875);
  setting_leg->SetLineColorAlpha(0, 0);
  setting_leg->SetFillStyle(0);
  for (const std::string& file : globals.GetDataFiles()) {
    // Set up the ROOT file and data pointer
    TFile *data_file = new TFile(file.c_str());		// Load the MAUS output file
    if ( !data_file->IsOpen() ) {
      Pitch::print(Pitch::error, "Histogram file not found: "+file);
      return 2;
    }

    // Get the input emittance and the setting name
    std::string base(file);
    base = base.substr(base.find("/")+1);

    std::string eps = base.substr(base.find("mm")-2, 2);
    if ( eps.find("_") != std::string::npos )
        eps = eps.substr(1);

    std::string setting = base;
    setting.erase(setting.end()-14-eps.size(), setting.end());

//    std::cerr << "\n" << eps << std::endl;
//    std::cerr << setting << std::endl;

    if ( std::find(settings.begin(), settings.end(), setting) == settings.end() )
	++id;

    // Fill graphs
    for (std::string& var : vars) {
      // Get the histogram that contains the information on this variable
      TGraphErrors* graph_ptr = (TGraphErrors*)data_file->Get(("data_"+var).c_str());      

      if ( !graph_ptr )
	  continue;

      // Initialize the graphs, set their style
      graphs[setting][eps][var] = new TGraphErrors();
      graphs[setting][eps][var]->SetMarkerStyle(markers[var]);
      graphs[setting][eps][var]->SetMarkerColor(colors[id+1]);
      graphs[setting][eps][var]->SetLineColor(colors[id+1]);
      graphs[setting][eps][var]->SetMarkerSize(1.5);
      mg->Add(graphs[setting][eps][var], "P");

      // Get the input and output variables and their errors
      graph_ptr->GetPoint(0, z, in);
      inerr = graph_ptr->GetErrorY(0);
      graph_ptr->GetPoint(5, z, out);
      outerr = graph_ptr->GetErrorY(5);     

      // Compute the change and the error, fill the graph
      val = 100*(out-in)/in;
      err = (val+100)*sqrt(pow(inerr/in, 2)+pow(outerr/out, 2));
      if ( var == "vol" ) {
//	val = 100*(1-sqrt(1-val/100));
	val /= 2;
	err /= 2.;
      }
//	std::cerr << var << "  " << val << "  " << err << std::endl;

//      std::cerr << settings.size() << std::endl;
      if ( setting == "2016-05_2" )
        graphs[setting][eps][var]->SetPoint(0, std::atof(eps.c_str())+0.1*((int)settings.size()-1), val);
      else
        graphs[setting][eps][var]->SetPoint(0, std::atof(eps.c_str())+0.1*((int)settings.size()-2), val);
      graphs[setting][eps][var]->SetPointError(0, 0, err);
    }

    // Close the graph file
    if ( std::find(settings.begin(), settings.end(), setting) == settings.end() ) {
      settings.push_back(setting);
      setting_leg->AddEntry(graphs[setting][eps]["subeps"], setting.c_str(), "l");
    }
    data_file->Close();
  }
  

  // Produce predictions for different settings
  double depth = 6.5;		// [cm], thickness of the absorber
  double density = 0.69;	// [g/cm^3], measured for the whole cylinder
  double Zeff = pow(.75*pow(3, 2.94)+.25*pow(1, 2.94), 1./2.94); // Effective atomic number of LiH
  double Aeff = Zeff/0.50321;					 // Effective atomitc mass
  double X_0 = 716.4*Aeff/(Zeff*(Zeff+1)*log(287/sqrt(Zeff)));	 // Radiation length of LiH [g/cm^2]
  X_0 /= density;

  double scatfunc = 13.6*sqrt(depth/X_0); // simp
  double bg = 140./105.66;
  double beta = bg/sqrt(1.+bg*bg);

  Material mat;
  mat.name = "LiH";
  mat.rho = density;
  mat.Z = Zeff;
  mat.A = Zeff/0.50321; // from PDG
  EnergyLoss eloss(mat, 10*depth); 

  
  std::vector<double> betas = {450, 600, 800};
  std::map<double, TF1*> feps;
  id = 0;
  colors = {kGreen+1, kCyan+1, kMagenta+1};
  for (const double& betaperp : betas) {
    feps[betaperp] = new TF1("feps", "100*[0]*sqrt(1+[1]*[2]/x)-100", 1, 20);
    feps[betaperp]->SetParameter(0, 1-eloss.GetMomentumLoss(140, 105.66)/140);
    feps[betaperp]->SetParameter(1, scatfunc*scatfunc/(105.66*beta*beta*140));
    feps[betaperp]->SetParameter(2, betaperp);
    feps[betaperp]->SetLineColor(colors[id]);
    ++id;
  }
//  std::cerr << feps->GetParameter(0) << (feps

  // Draw the multigraph
  TCanvas *c = new TCanvas("c", "c", 1200, 800);
  mg->Draw("A");
  mg->GetXaxis()->SetLimits(4, 12);

  TLine *zeroline = new TLine(mg->GetXaxis()->GetXmin(), 0, mg->GetXaxis()->GetXmax(), 0);
  zeroline->SetLineColor(kGray);
  zeroline->Draw("SAME");

  InfoBox info("Preliminary", "2.9.1", "sol. + flip", "2016/04-05");
  info.Draw();

  for (const double& betaperp : betas)
    feps[betaperp]->Draw("SAME");

  setting_leg->Draw("SAME");
  c->SaveAs("comparison.pdf");
  delete c;
}
