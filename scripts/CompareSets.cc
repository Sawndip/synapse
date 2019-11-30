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

double ThFunction(double* x, double* par) {

  if ( x[0] > 5.5 && x[0] < 6.5 )
      return 100*par[0]*sqrt(1+par[1]*par[2]/6)-100;

  if ( x[0] > 9.5 && x[0] < 10.5 )
      return 100*par[0]*sqrt(1+par[1]*par[2]/10)-100;

  return 100*par[0]*sqrt(1+par[1]*par[2]/x[0])-100;
}

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
      if ( std::string(file).find(".txt") > std::string(file).size() ) {
        Pitch::print(Pitch::error, "Invalid data file type: "+std::string(file));
	return 2;
      }
  }

  // Style settings
  std::vector<std::string> settings;
//  std::vector<std::string> vars = {"amp", "subeps", "vol", "de"};
//  std::vector<std::string> vars = {"amp", "subeps", "vol"};
//  std::vector<std::string> vars = {"subeps"};
  std::vector<std::string> vars = {"vol", "vol_de"};
  std::map<std::string, int> markers = {{"amp",27}, {"subeps",24}, {"vol",25}, {"vol_de",26}};
  std::map<std::string, double> colors =
	{{"2017/02-7",kBlue+1}, {"2016/05-1",kRed+1}, {"2016/04-1.2",kBlack}};

  // Initialize the graphs
  std::map<std::string, std::map<int, std::map<std::string, TGraphErrors*>>> graphs, fullgraphs;
  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle(";#epsilon_{#perp}^{u}  [mm];#Delta#epsilon_{#perp}  /#epsilon_{#perp}^{u}  [%]");

  // Loop over the input ROOT files and fill the graphs
  std::map<std::string, std::map<int, std::map<std::string, double>>> changes, errors, syst_errors;
  double change, error, syst_error;
  std::string var, abs;
  int id(-1), eps;
  for (const std::string& file : globals.GetDataFiles()) {
    // Set up the text file
    /*TFile *data_file = new TFile(file.c_str());		// Load the MAUS output file
    if ( !data_file->IsOpen() ) {
      Pitch::print(Pitch::error, "Histogram file not found: "+file);
      return 2;
    }*/
    std::ifstream data_file;
    data_file.open(file);

    // Get the setting name, add a new category if it does not exist
    std::string setting(file);
    setting = setting.substr(0, setting.size()-4);
    bool isde(false);
    if ( setting.find("de") != std::string::npos ) {
      setting = setting.substr(11);
      isde = true;
    } else {
      setting = setting.substr(8);
    }
    setting.replace(4, 1, "/");
    setting.replace(7, 1, "-");
    if ( setting.size() > 9 )
        setting.replace(9, 1, ".");

    if ( std::find(settings.begin(), settings.end(), setting) == settings.end() )
	++id;
    if ( std::find(settings.begin(), settings.end(), setting) == settings.end() )
      settings.push_back(setting);

    // Extract all the information from the text file
    changes.clear();
    errors.clear();
    syst_errors.clear();
    while ( data_file >> var >> eps >> abs >> change >> error >> syst_error ) {
      if ( var == "vol" ) {
//	if ( isde )
//	    continue;
	change = sqrt(1.+change)-1.;
	error /= 2;
      } else if (var == "den") {
	std::cerr << eps << "  " << change << std::endl;
	change = 1./sqrt(1.+change)-1.;
	std::cerr << change << std::endl;
	error /= 2;
      }
      if ( isde )
	var += "_de";
      changes[var][eps][abs] = 1e2*change;
      errors[var][eps][abs] = 1e2*error;
      syst_errors[var][eps][abs] = 1e2*syst_error;
    }

    // Fill graphs
    for (const std::string& var : vars) {
      for (const int& eps : {6, 10}) {

	// Skip if it wasn't found in the file
	if ( changes.find(var) == changes.end() )
	    continue;

        // Initialize the graphs once, set their style
	if ( !graphs[setting][eps][var] ) {
          graphs[setting][eps][var] = new TGraphErrors();
          graphs[setting][eps][var]->SetMarkerStyle(markers[var]);
          graphs[setting][eps][var]->SetLineColor(colors[setting]);
          graphs[setting][eps][var]->SetLineWidth(3);
          graphs[setting][eps][var]->SetMarkerSize(2);
          fullgraphs[setting][eps][var] = new TGraphErrors(*graphs[setting][eps][var]);
          fullgraphs[setting][eps][var]->SetLineWidth(1.5);
          mg->Add(graphs[setting][eps][var], "P");
          mg->Add(fullgraphs[setting][eps][var], "P");
        }

	double sign = - (var == "vol") + (var == "den");
        double x = eps+0.05/*+0.3*((int)settings.size()-2)*/+sign*0.2;
	double full_error = sqrt(pow(errors[var][eps]["LiH"], 2)+pow(syst_errors[var][eps]["LiH"], 2));
        graphs[setting][eps][var]->SetPoint(0, x, changes[var][eps]["LiH"]);
        graphs[setting][eps][var]->SetPointError(0, 0, errors[var][eps]["LiH"]);
        fullgraphs[setting][eps][var]->SetPoint(0, x, changes[var][eps]["LiH"]);
        fullgraphs[setting][eps][var]->SetPointError(0, 0, full_error);
      }
    }

    // Close the graph file
    data_file.close();
  }

  // Initialize the legend
  TLegend *setting_leg = new TLegend(.5, .88-.05*settings.size(), .725, .89);
  setting_leg->SetLineColorAlpha(0, 0);
  setting_leg->SetFillStyle(0);
  for (const std::string& setting : settings)
      setting_leg->AddEntry(graphs[setting][6][vars[0]], setting.c_str(), "l");

  TLegend *type_leg = new TLegend(.125, .125, .35, .4);
  type_leg->SetLineColorAlpha(0, 0);
  type_leg->SetFillStyle(0);
  std::vector<std::string> names = {"Amplitude", "kNN DE"};
  size_t ii = 0;
  for (const std::string& var : vars) {
      type_leg->AddEntry(graphs[settings[0]][6][var], names[ii].c_str(), "p");
      ii++;
  }
 
  TGraphErrors *gdummythick = new TGraphErrors();
  gdummythick->SetLineWidth(3);
  type_leg->AddEntry(gdummythick, "Statistical", "e");

  TGraphErrors *gdummy = new TGraphErrors();
  gdummy->SetLineWidth(1.5);
  type_leg->AddEntry(gdummy, "Total", "e");

  TF1 *fdummy = new TF1("fdummy", "0.", 0, 1);
  fdummy->SetLineWidth(1);
  fdummy->SetLineColor(kBlack);
  type_leg->AddEntry(fdummy, "Model", "l");

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

  EnergyLoss eloss(AbsLiH, 10*depth); 

  std::map<std::string, double> betas = {{"2017/02-7",600}, {"2016/05-1",650}, {"2016/04-1.2",900}};
  std::map<double, TF1*> feps;
//  colors = {kGreen+1, kCyan+1, kMagenta+1};
  double betaperp;
  for (const std::string& setting : settings) {
    betaperp = betas[setting];
    feps[betaperp] = new TF1("feps", ThFunction, 1, 20, 3);
    feps[betaperp]->SetParameter(0, 1-eloss.GetMomentumLoss(140, 105.66)/140);
    feps[betaperp]->SetParameter(1, scatfunc*scatfunc/(105.66*beta*beta*140));
    feps[betaperp]->SetParameter(2, betaperp);
    feps[betaperp]->SetLineWidth(1);
    feps[betaperp]->SetLineColor(colors[setting]);
    feps[betaperp]->SetNpx(1e3);
  }

  // Draw the multigraph
  TCanvas *c = new TCanvas("c", "c", 1200, 800);
  mg->Draw("A");
  mg->GetYaxis()->SetTitleOffset(.8);
  mg->GetXaxis()->SetLimits(2, 12);
  mg->SetMaximum(3);
  mg->SetMinimum(-12);

  // Draw a band for each emittance setting
  for (const double eps : {/*3, */6, 10}) {
    TBox *box = new TBox(eps-.5, -12, eps+.5, 3);
    box->SetFillColorAlpha(kBlack, .15);
    box->Draw("FSAME");
  }
      

  // Draw a line for 0 emittance change
  TLine *zeroline = new TLine(mg->GetXaxis()->GetXmin(), 0, mg->GetXaxis()->GetXmax(), 0);
  zeroline->SetLineColorAlpha(kBlack, .5);
  zeroline->Draw("SAME");

  // Draw the info box
  InfoBox info("Internal", "3.2.0", "sol. + flip", "2016/2017");
  info.Draw();

  // Draw the functions
  for (const std::string& setting : settings)
      feps[betas[setting]]->Draw("SAME");

  setting_leg->Draw("SAME");
  type_leg->Draw("SAME");
  c->SaveAs("cool_summary.pdf");
  delete c;
}
