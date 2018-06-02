// Cpp includes
#include <ctime>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <regex>
#include <map>
#include <chrono>

// Root includes
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TObject.h"
#include "TString.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TF1.h"
#include "TVector3.h"
#include "TMath.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TRandom3.h"
#include "TCanvas.h"

// Additional modules
#include "Pitch.hh"
#include "ProgressBar.hh"
#include "InfoBox.hh"
#include "Bunch.hh"
#include "Scattering.hh"
#include "EnergyLoss.hh"
#include "Transport.hh"
#include "ToyTools.hh"
#include "DChiSquared.hh"
#include "Globals.hh"

/** @file  ToySimulation.cc
 *
 *  @brief Runs a toy simulation of the beam line.
 *
 *	   Measures the evolution of phase space parameters across a toy absorber.
 *
 *	   Supports the addition of drift spaces, solenoids and apertures.
 **/

double AmplitudeDistribution(double *x, double *par) {

  // Return the theroretical distribution of the amplitude
  DChiSquared fchi2(4);
  return par[0]*fchi2(x[0]/par[1])/par[1];
}

double AmplitudeChangeDistribution(double *x, double *par) {

  // Return the theroretical distribution of the amplitude
  double delta = (par[2]-par[1])/par[1];
  return par[0]*x[0]*(exp(-x[0]/(2*par[2]))/pow(1+delta, 2)-exp(-x[0]/(2*par[1])))/(4*pow(par[1], 2));
}

double RelativeAmplitudeChangeDistribution(double *x, double *par) {

  // Return the theroretical distribution of the amplitude
  DChiSquared fchi2(4);
  return AmplitudeChangeDistribution(x, par)/(par[0]*fchi2(x[0]/par[1])/par[1]);
}

double ChiSquaredDist(double x, double d) {

  return  exp(-x/2)*pow(x, d/2-1)/pow(2, d/2)/tgamma(d/2);
}

/** @brief	Main function
 *
 *	       ./toy_sim [options] [external_input.root]
 *
 *  @param	argc		Number of command line arguments
 *  @param	argv		List of command line aguments
 */
int main(int argc, char ** argv) {

  // Reroot the info messages to the terminal for now
  Pitch::setAnOutput(Pitch::info, std::cerr);

  // Toy MC algorithm global parameters (stored in globals)
  Globals &globals = Globals::GetInstance(argc, argv);
  if ( globals.GetDataFiles().size() ) {
    for (const std::string& file : globals.GetDataFiles())
      if ( std::string(file).find(".root") > std::string(file).size() ) {
        Pitch::print(Pitch::error, "Invalid data file type: "+std::string(file));
	return 2;
      }
  }

  // Set the Stat box to be along the top right corner of the TPad in the TCanvas
  gSystem->Load("libTree");
  gStyle->SetTitleSize(.04, "XY");
  gStyle->SetLabelSize(.04, "XY");
  gStyle->SetStatX(.9);
  gStyle->SetStatY(.9);
  gStyle->SetOptStat(2210);
  gStyle->SetOptFit(1);

  // Toy Monte Carlo variables
  std::vector<std::string> vars = {"x", "y", "px", "py", "pz"};
  
  // Labels of the cooling variables
  double frac = globals["frac"];
  std::string frac_str = std::to_string((int)(100*frac));
  std::vector<std::string> pars = {"meanp", "eps", "betax", "betay", "beta", "trans",
	"subeps", "acut", "vol", "mode", "quan"};
  std::vector<std::string> fracpars = {"acut", "subeps", "vol"};
  std::map<std::string, std::string> par_labels = {{"meanp","#bar{p}"}, {"eps","#epsilon_{n}"}, 
	{"betax","#beta_{x}"}, {"betay","#beta_{y}"}, {"beta","#beta_{#perp}"}, {"trans","T"},
	{"subeps","e"}, {"acut","A"}, {"vol","#epsilon"}, {"mode","A_{#perp}^{*}"}, {"quan","P"}};
  std::map<std::string, std::string> par_units = {{"meanp","MeV/c"}, {"eps","mm"},
	{"betax","mm"}, {"betay","mm"}, {"beta","mm"}, {"trans","%"}, {"subeps","mm"},
	{"acut","mm"}, {"vol","mm^{2}MeV^{2}/c^{2}"}, {"mode","mm"}, {"quan","mm"}};
  std::map<std::string, std::string> par_titles =
	{{"meanp","Mean total momentum"},
	 {"eps","Normalised transverse RMS emittance"},
	 {"betax","Optical beta function in x"},
	 {"betay","Optical beta function in y"},
	 {"beta","Transverse beta function"},
	 {"trans","Tansmission"},
	 {"subeps","Subsample emittance"},
	 {"acut","Amplitude cut"},
	 {"vol","Fractional emittance"},
	 {"mode","Most probable amplitude"},
	 {"quan","Amplitude quantile"}};
  std::map<std::string, std::string> par_names =
	{{"eps","Emittance"}, {"subeps","Subemittance"}, {"vol","Fractional"}};

  // Calculate the LiH parameters in MICE
  double depth = 6.5;		// [cm], thickness of the absorber
  double density = 0.693;	// [g/cm^3], measured for the whole cylinder
  double Z = 4;			// 3+1	
  double A = 7.05274;		// [g/mol], 95.52% of Li6 (7), 4.48% of Li7 (7), natural H (1.00794)
  double X_0 = 102;	 	// [cm], radiation length
  double I = 36.5;		// [eV], mean excitation potential
  Scattering scat(depth/X_0);

  Material mat;
  mat.name = "LiH";
  mat.rho = density;
  mat.Z = Z; 		
  mat.A = A; 	
  mat.I = I;
  EnergyLoss eloss(mat, 10*depth);
  std::cerr << eloss.GetMomentumLoss(140, 105.66) << std::endl;
  std::cerr << eloss.GetMostProbableMomentumLoss(140, 105.66) << std::endl;

  // Print the changes as a function of beta for 3 different input emittance
  std::vector<size_t> colors = {kBlack, kRed+2, kGreen+2, kBlue+2, kMagenta+2, kCyan+2, kYellow+2};
  /*std::vector<double> vbeta = {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};
  std::map<size_t, std::vector<double>> vdeps;
  vdeps[3] = {-10.04,-8.49,-6.84,-5.19,-3.57,-1.95,-0.37,1.19,2.72,4.23};
  vdeps[6] = {-11.25,-10.5,-9.5,-8.81,-7.92,-7.05,-6.19,-5.34,-4.50,-3.67};
  vdeps[10] = {-12,-12,-11.16,-10.46,-9.83,-9.22,-8.62,-8.01,-7.41,-6.82};
  std::map<size_t, TGraphErrors*> gdeps;
  gdeps[3] = new TGraphErrors(10);
  gdeps[6] = new TGraphErrors(10);
  gdeps[10] = new TGraphErrors(10);
  for (size_t i = 0; i < 10; i++) {
    for (size_t eps : {3,6,10} ) {
      gdeps[eps]->SetPoint(i, vbeta[i], vdeps[eps][i]);
      gdeps[eps]->SetPointError(i, 0, 100*sqrt((pow(1+vdeps[eps][i]/100, 2)-0.88*0.88)*2./1e5));
    }
  }
  size_t id(0);
  for (auto pair : gdeps) {
    pair.second->SetMarkerStyle(20+id);
    pair.second->SetMarkerSize(1.5);
    pair.second->SetLineWidth(2);
    pair.second->SetLineColor(colors[id]);
    id++;
  }*/


  std::map<size_t, TF1*> fbeta;
  std::vector<size_t> epss = {3, 6, 10};
  size_t id = 0;
  TCanvas* canv0 = new TCanvas("c", "c", 1200, 800);
  gPad->SetGridx();
  gPad->SetGridy();
  TLegend* leg0 = new TLegend(.15, .7, .35, .85);
  leg0->SetLineColor(0);
  for (const double& eps : epss) {
    fbeta[eps] = new TF1(TString::Format("fbeta_%d", (int)eps), "100*[0]*sqrt(1+[1]*x/[2])-100", 0, 1025);
    fbeta[eps]->SetTitle(";#beta_{#perp}  [mm];#Delta#epsilon_{#perp}  /#epsilon_{#perp}^{i}  [%]");
    fbeta[eps]->SetParameter(0, 1-eloss.GetMomentumLoss(140, 105.66)/140);
    fbeta[eps]->SetParameter(1, 1e-3);
    fbeta[eps]->SetParameter(2, eps);
    fbeta[eps]->SetLineColor(colors[id]);

    if ( eps == 3 )
      fbeta[eps]->Draw("");
    fbeta[eps]->Draw("SAME");

    TFile ifile(TString::Format("beta_eps_%dmm.root", (int)eps), "READ");
    TGraphErrors* graph = (TGraphErrors*)ifile.Get("");
    graph->SetLineColor(colors[id]);
    graph->SetMarkerStyle(24+id);
    graph->Draw("PE SAME");
    leg0->AddEntry(graph, TString::Format("#epsilon_{#perp}^{i}  = %d mm", (int)eps), "lep"); 

    ++id;
  }
  leg0->Draw("SAME");
  canv0->SaveAs("deps_beta_toy.pdf");
  delete canv0;

  // Print the bias on emittance change as a function of transmission
  std::vector<size_t> lengths = {0, 174, 372, 639, 1137, 5986};
  std::map<size_t, double> corr = {{0,0}, {86,0.1}, {174,0.2}, {268,0.3}, {372,0.4}, {492,0.5},
	{639,0.6},{836,0.7},{1137,0.8},{1761,0.9},{5986,0.99}};
  id = 0;
  TCanvas* canv1 = new TCanvas("c", "c", 1200, 800);
  gPad->SetGridx();
  gPad->SetGridy();
  TLegend* leg1 = new TLegend(.6, .15, .85, .4);
  leg1->SetLineColor(0);
  TFile ifile1("trans_amp.root", "READ");
  TGraphErrors* graph1 = (TGraphErrors*)ifile1.Get("");
  graph1->SetTitle(";Transmission [%];#Delta#epsilon_{#perp}  /#epsilon_{#perp}^{i}  [%]");
  graph1->SetLineColor(kBlack);
  graph1->SetMarkerStyle(20);
  graph1->SetMarkerSize(1.5);
  leg1->AddEntry(graph1, "Amplitude", "ep"); 
  graph1->Draw("APE");
  for (const double& len : lengths) {

    TFile ifile(TString::Format("trans_rad_%d.root", (int)len), "READ");
    TGraphErrors* graph = (TGraphErrors*)ifile.Get("");
    graph->SetLineColor(colors[id+1]);
    graph->SetMarkerStyle(21+id);
    graph->SetMarkerColor(colors[id+1]);
    graph->SetMarkerSize(1.5);
    graph->Draw("PE SAME");
    std::cerr << len << " " << 100*corr[len] << std::endl;
    std::cerr << len << " " << 100*corr[len] << std::endl;
    leg1->AddEntry(graph, TString::Format("Radial, #rho_{xp}  = %d %%", (int)(100*corr[len])), "ep"); 

    ++id;
  }
  leg1->Draw("SAME");
  canv1->SaveAs("deps_trans_toy.pdf");
  delete canv1;

  // Initialize the transporter
  Transport transport(4);
  if ( globals["toy_drift"] )
      transport.AddDriftSpace(globals["toy_drift"], globals["toy_mom"]);
  if ( globals["toy_sol"] )
      transport.AddSolenoid((double)globals["toy_sol"]/100, .5);

  // Define the input beam. Get the input sample from an external file if specified
  size_t Size = (double)globals["toy_n"];
  std::map<std::string, std::vector<double>> iSize;
  size_t ref_id = globals["ref_id"];
  if ( globals.GetDataFiles().size() ) {
    Pitch::print(Pitch::info, "Importing"+globals.GetDataFiles()[0]);
    TFile data_file(globals.GetDataFiles()[0].c_str());	// Load the MAUS input file
    if ( !data_file.IsOpen() ) {
      Pitch::print(Pitch::error, "Data file not found: "+globals.GetDataFiles()[0]);
      return 2;
    }

    TNtuple* samples = (TNtuple*)data_file.Get("Truth");
    float* ntuple;
    size_t vid;
    TVector3 mom, pos;
    ProgressBar pbar;
    for (size_t i = 0; i < (size_t)samples->GetEntries(); ++i) {

      // Display the progress in %
      pbar.GetProgress(i, (size_t)samples->GetEntries());

      // Fetch the variables from the Ntuple
      samples->GetEntry(i);
      ntuple = samples->GetArgs();
      vid = ntuple[2];
      if ( vid != ref_id )
	  continue;
      pos = TVector3(ntuple[3], ntuple[4], ntuple[5]);
      mom = TVector3(ntuple[6], ntuple[7], ntuple[8]);
	
      // Fill the external samples
      iSize["x"].push_back(pos.x());
      iSize["y"].push_back(pos.y());
      iSize["px"].push_back(mom.x());
      iSize["py"].push_back(mom.y());
      iSize["pz"].push_back(mom.z());
    }
  } else {
    iSize = GaussianBunch(globals["toy_mass"], globals["toy_mom"], Size,
			      globals["toy_eps"], globals["toy_beta"], globals["toy_alpha"]);
  }

  // Draw the input beam
  Beam::Bunch inbeam(iSize, 0, "in");
  TH2F* hinxpx = inbeam.Histogram("x", "px", -400, 400, -100, 100);
  TEllipse* inxpxell = inbeam.Ellipse("x", "px");

  TCanvas* canv = new TCanvas("c", "c", 1200, 800);
  hinxpx->Draw("COLZ");
  inxpxell->Draw("SAME");
  canv->SaveAs("xpx_in.pdf");
  delete canv;

  // Define and draw the density estimation of the input beam
  if ( globals["de"] ) {
    std::vector<std::vector<double>> samples(iSize["x"].size());
    for (size_t i = 0; i < samples.size(); i++) {
      samples[i].push_back(iSize["x"][i]);
      samples[i].push_back(iSize["px"][i]);
      samples[i].push_back(iSize["y"][i]);
      samples[i].push_back(iSize["py"][i]);
    }

    DensityEstimator de(samples, "knn", false);  
    TH2F* hxpx = new TH2F("xpx_in_de", "in xp_{x} density estimation;x [mm];p_{x} [MeV/c]",
			  100, -400, 400, 100, -100, 100);
    for (size_t i = 0; i < (size_t)hxpx->GetNbinsX(); i++)
      for (size_t j = 0; j < (size_t)hxpx->GetNbinsY(); j++)
	  hxpx->SetBinContent(i+1, j+1, de({hxpx->GetXaxis()->GetBinCenter(i+1),
					    hxpx->GetYaxis()->GetBinCenter(j+1), 0, 0}));

    gStyle->SetOptStat(0);
    TCanvas *cbeam = new TCanvas("c", "c", 1200, 800);
    hxpx->Draw("CONTZ");
    cbeam->SaveAs("xpx_in_de.pdf");
    delete cbeam;
  }

  // Diffuse and scatter the input beam to create the output beam
  std::map<std::string, std::vector<double>> outsamples = iSize;
  if ( globals["toy_eloss"] )
      eloss.IonizeBunch(outsamples, globals["toy_mom"], globals["toy_mass"]);
  if ( globals["toy_scat"] )
      scat.ScatterBunch(outsamples, globals["toy_mom"], globals["toy_mass"]);
  if ( globals["toy_drift"] || globals["toy_sol"] )
      transport.TransportBunch(outsamples);

  // Define the output beam
  double outpos = 10*depth+(int)globals["toy_drift"]+(int)globals["toy_sol"];
  Beam::Bunch outbeam(outsamples, outpos, "out");

  // Draw the output beam
  gStyle->SetOptStat(1);
  TH2F* houtxpx = outbeam.Histogram("x", "px", -400, 400, -100, 100);
  TEllipse* outxpxell = outbeam.Ellipse("x", "px");
  Matrix<double> covmat = outbeam.CovarianceMatrix();
  covmat.Print();
  std::cerr << "corr: " << covmat[0][2]/sqrt(covmat[0][0]*covmat[2][2]) << std::endl;
  std::cerr << sqrt(140*140*covmat[2][2]/covmat[0][0]) << std::endl;

  canv = new TCanvas("c", "c", 1200, 800);
  houtxpx->Draw("COLZ");
  outxpxell->Draw("SAME");
  canv->SaveAs("xpx_out.pdf");
  delete canv;

  // Define and draw the density estimation of the input beam
  if ( globals["de"] ) {
    std::vector<std::vector<double>> samples(outsamples["x"].size());
    for (size_t i = 0; i < samples.size(); i++) {
      samples[i].push_back(outsamples["x"][i]);
      samples[i].push_back(outsamples["px"][i]);
      samples[i].push_back(outsamples["y"][i]);
      samples[i].push_back(outsamples["py"][i]);
    }

    DensityEstimator de(samples, "knn", false);  
    TH2F* hxpx = new TH2F("xpx_out_de", "out xp_{x} density estimation;x [mm];p_{x} [MeV/c]",
			  100, -400, 400, 100, -100, 100);
    for (size_t i = 0; i < (size_t)hxpx->GetNbinsX(); i++)
      for (size_t j = 0; j < (size_t)hxpx->GetNbinsY(); j++)
	  hxpx->SetBinContent(i+1, j+1, de({hxpx->GetXaxis()->GetBinCenter(i+1),
					    hxpx->GetYaxis()->GetBinCenter(j+1), 0, 0}));

    gStyle->SetOptStat(0);
    TCanvas *cbeam = new TCanvas("c", "c", 1200, 800);
    hxpx->Draw("CONTZ");
    cbeam->SaveAs("xpx_out_de.pdf");
    delete cbeam;
  }

  // Plot the energy loss
  gStyle->SetOptStat(2210);
  TH1F* heloss = new TH1F("eloss", "Energy loss going through the absorber", 100, 0, 30);
  Vector<double> vin, vout;
  for (size_t i = 0; i < Size; i++) {
    vin = Vector<double>({iSize["px"][i], iSize["py"][i], iSize["pz"][i]});
    vout = Vector<double>({outsamples["px"][i], outsamples["py"][i], outsamples["pz"][i]});
    heloss->Fill(vin.mag()-vout.mag());
  }

  canv = new TCanvas("c", "c", 1200, 800);
  heloss->Draw("E");
  canv->SaveAs("eloss_toy.pdf");
  delete canv;

  // Draw the in and out ellipses to compare them on a single canvas
  TH2F* hempty = new TH2F("empty", "1 #sigma ellipses;x [mm]; p_{x} [MeV/c]",
			  100, -100, 100, 100, -25, 25);
  gStyle->SetOptStat(0);
  canv = new TCanvas("c", "c", 1200, 800);
  hempty->Draw("");
  inxpxell->Draw("SAME");
  inxpxell->SetLineColor(4);
  outxpxell->Draw("SAME");
  canv->SaveAs("xpx_ellipses.pdf");
  delete canv;

  // Compute and print the emittance change
  double ieps = inbeam.NormEmittance().GetValue();
  double oeps = outbeam.NormEmittance().GetValue();
  Pitch::print(Pitch::info, "Emittances: "+std::to_string(ieps)+" mm (in),  "
					  +std::to_string(oeps)+" mm (out)");
  Pitch::print(Pitch::info, "Emittance change: "+std::to_string(100*(oeps-ieps)/ieps)+" %");
  if ( globals["de"] ) {
    double iepsest = inbeam.NormEmittanceEstimate(Beam::vol).GetValue();
    double oepsest = outbeam.NormEmittanceEstimate(Beam::vol).GetValue();
    Pitch::print(Pitch::info, "Estimated emittances: "+std::to_string(iepsest)+" mm (in),  "
					    	      +std::to_string(oepsest)+" mm (out)");
    Pitch::print(Pitch::info, "Estimated emittance change: "
	+std::to_string(100*(oepsest-iepsest)/iepsest)+" %");
  }

  // Get the arrays of amplitudes in and out
  std::vector<double> inamps, outamps;
  inamps = inbeam.Amplitudes();
  outamps = outbeam.Amplitudes();

  // Set up the theoretical amplitude distributions
  TF1* finamps = new TF1("finamps", AmplitudeDistribution, 0, 100, 2);
  finamps->SetParameters(5*Size, ieps);
  finamps->SetLineColor(4);
  finamps->SetLineWidth(1);

  TF1* foutamps = new TF1("foutamps", AmplitudeDistribution, 0, 100, 2);
  foutamps->SetParameters(5*Size, oeps);
  foutamps->SetLineWidth(1);

  // Draw the two amplitude distributions
  TH1F* hinamps = new TH1F("inamps", "Transverse amplitude;A_{#perp} [mm]", 20, 0, 100);
  hinamps->SetLineWidth(2);
  hinamps->SetMarkerStyle(24);
  hinamps->SetMarkerSize(1.5);
  hinamps->FillN(inamps.size(), &(inamps[0]), NULL);
  TH1F* houtamps = new TH1F("outamps", "Transverse amplitude;A_{#perp} [mm]", 20, 0, 100);
  houtamps->SetLineWidth(2);
  houtamps->SetLineColor(2);
  houtamps->SetMarkerStyle(25);
  houtamps->SetMarkerSize(1.5);
  houtamps->FillN(outamps.size(), &(outamps[0]), NULL);

  THStack *hamps_stack = new THStack("amp_stack", "Transverse amplitude;A_{#perp}  [mm]");
  hamps_stack->Add(hinamps);
  hamps_stack->Add(houtamps);

  TLegend* leg = new TLegend(.7, .55, .89, .69);
  leg->SetLineColor(0);
  leg->SetFillColor(0);
  leg->AddEntry(hinamps, "up", "lpe");
  leg->AddEntry(houtamps, "down", "lpe");

  InfoBox info("toy", "2.9.1");

  canv = new TCanvas("c", "c", 1200, 800);
  hamps_stack->Draw("PE NOSTACK");
  info.Draw();
  leg->Draw("SAME");
  canv->SaveAs("amp_toy.pdf");

  finamps->Draw("SAME");
  foutamps->Draw("SAME");
  canv->SaveAs("amp_toy_fitted.pdf");
  delete canv;

  // Draw the amplitude change distribution
  Vector<double> vamps = Vector<double>(outamps)-Vector<double>(inamps);
  TH2F* hvamp = new TH2F("vamp",
	"Transverse amplitude change;A_{#perp}^{U} [mm];A_{#perp}^{D}-A_{#perp}^{U} [mm]",
	20, 0, 100, 25, -25, 25);
  hvamp->FillN(inamps.size(), &inamps[0], &vamps[0], NULL, 1);

  canv = new TCanvas("c", "c", 1200, 800);
  hvamp->Draw("COLZ");
  Beam::Drawer drawer;
  TH1F* vamp_x = drawer.ProjectionMeanX(hvamp);
  vamp_x->SetLineColor(2);
  vamp_x->SetLineWidth(3);
  vamp_x->SetMarkerStyle(21);
  vamp_x->Draw("P SAME");

  info.Draw();
  canv->SaveAs(std::string("vamp_toy.pdf").c_str());
  delete vamp_x;
  delete canv;

  // Draw the relative change in density in each bin
  TF1* frelampdiff = new TF1("frelampdiff", RelativeAmplitudeChangeDistribution, 1e-2, 100, 3);
  frelampdiff->SetParameters(5*Size, ieps, oeps);
  frelampdiff->SetLineWidth(1);
  std::cerr << frelampdiff->Eval(1e-2) << std::endl;

  canv = new TCanvas("c", "c", 1200, 800);
  hinamps->Sumw2();
  houtamps->Sumw2();
  houtamps->SetLineColor(9);
  houtamps->Add(hinamps, -1);
  houtamps->Divide(hinamps);
  houtamps->SetTitle("Fractional change in density");
  houtamps->GetYaxis()->SetTitle("(N_{D}-N_{U})/N_{U} [%]");
  houtamps->Draw("E");
  frelampdiff->Draw("SAME");

  info.Draw();
  canv->SaveAs("famp_toy.pdf");
  delete canv;

  // Compare the Voronoi cell volumes (TODO TODO TODO)
  if ( globals["voronoi"] ) {
    inbeam.SetVoronoiVolumes();
    Pitch::print(Pitch::info, "Computed upstream Voronoi volumes");
    std::vector<double> invvols = inbeam.VoronoiVolumes();
    outbeam.SetVoronoiVolumes();
    Pitch::print(Pitch::info, "Computed downstream Voronoi volumes");
    std::vector<double> outvvols = outbeam.VoronoiVolumes();

    double minin = *std::min_element(invvols.begin(), invvols.end());
    double minout = *std::min_element(outvvols.begin(), outvvols.end());
    std::cerr << minin << "  " << minout << std::endl;
    std::cerr << "means: " << Math::Mean(invvols) << "  " << Math::Mean(outvvols) << std::endl;
    for (size_t i = 0; i < invvols.size(); i++) {
//	invvols[i] = pow(2*invvols[i], .5)/105.66/M_PI;
//	outvvols[i] = pow(2*outvvols[i], .5)/105.66/M_PI;
	invvols[i] = log(invvols[i]);
	outvvols[i] = log(outvvols[i]);
    }

    TCanvas *ccomp = new TCanvas("c", "c", 1200, 800);
    TH2F* hincomp = new TH2F("incomp", ";A_{#perp}  [mm];v_{#perp}  [mm^{2}MeV^{2}/c^{2}]",
	20, 0, 100, 20, 0, 20);
    hincomp->FillN(inamps.size(), &inamps[0], &invvols[0], NULL, 1);
    hincomp->Draw("COLZ");
    ccomp->SaveAs("incomp.pdf");
    delete ccomp;

    ccomp = new TCanvas("c", "c", 1200, 800);
    TH2F* houtcomp = new TH2F("outcomp", ";A_{#perp}  [mm];v_{#perp}  [mm^{2}MeV^{2}/c^{2}]",
	20, 0, 100, 20, 0, 20);
    houtcomp->FillN(outamps.size(), &outamps[0], &outvvols[0], NULL, 1);
    houtcomp->Draw("COLZ");
    ccomp->SaveAs("outcomp.pdf");
    delete ccomp;

    TH1F* hinvvols = new TH1F("invvols", ";v_{#perp}  [mm^{2}MeV^{2}/c^{2}]", 20, 0, 20);
    hinvvols->SetLineWidth(2);
    hinvvols->SetMarkerStyle(24);
    hinvvols->SetMarkerSize(1.5);
    hinvvols->FillN(invvols.size(), &(invvols[0]), NULL);
    TH1F* houtvvols = new TH1F("outvvols", "v_{#perp}  [mm^{2}MeV^{2}/c^{2}]", 20, 0, 20);
    houtvvols->SetLineWidth(2);
    houtvvols->SetLineColor(2);
    houtvvols->SetMarkerStyle(25);
    houtvvols->SetMarkerSize(1.5);
    houtvvols->FillN(outvvols.size(), &(outvvols[0]), NULL);

    THStack *hvvols_stack = new THStack("amp_stack", ";v_{#perp}  [mm^{2}MeV^{2}/c^{2}]");
    hvvols_stack->Add(hinvvols);
    hvvols_stack->Add(houtvvols);

    leg = new TLegend(.7, .55, .89, .69);
    leg->SetLineColor(0);
    leg->SetFillColor(0);
    leg->AddEntry(hinvvols, "up", "lpe");
    leg->AddEntry(houtvvols, "down", "lpe");

    canv = new TCanvas("c", "c", 1200, 800);
    hvvols_stack->Draw("PE NOSTACK");
    info.Draw();
    leg->Draw("SAME");
    canv->SaveAs("vvol_toy.pdf");

  //  finvvols->Draw("SAME");
  //  foutvvols->Draw("SAME");
    canv->SaveAs("vvol_toy_fitted.pdf");
    delete canv;
  }

  // Plot the fractional quantities as a function of the fraction
  Pitch::print(Pitch::info, "Filling the fractional graphs");
  Beam::Bunch tempin, tempout;
  std::map<std::string, TGraphErrors*> fracgraphs;
  for (const std::string& fracpar : fracpars)
      fracgraphs[fracpar] = new TGraphErrors();

  double maxfrac(1.), minfrac(.1), fracstep(.05);
  double tempfrac, in, out, transmission, ine, oute, value, error;
  ProgressBar pbar;
  for (tempfrac = minfrac; tempfrac <= maxfrac+0.01; tempfrac += fracstep) {

    transmission = (double)outbeam.Size()/inbeam.Size();
    if ( transmission < tempfrac )
	continue;

    inbeam.SetCoreFraction(tempfrac);
    outbeam.SetCoreFraction(tempfrac/transmission);    
    for (const std::string& fracpar : fracpars )  {
      if ( fracpar == "subeps" ) {
        in = tempin.NormEmittance().GetValue();
        out = tempout.NormEmittance().GetValue();
//        in = tempin.Mean("amp").GetValue();
//        out = tempout.Mean("amp").GetValue();
        ine = tempin.NormEmittance().GetError();
        oute = tempout.NormEmittance().GetError();
      } else if ( fracpar == "acut" ) {
        in = Math::Max(tempin.Amplitudes());
        ine = in/sqrt(tempin.Size());
        out = Math::Max(tempout.Amplitudes());
        oute = out/sqrt(tempin.Size());
      } else if ( fracpar == "vol" ) {
	if ( !globals["de"] ) {
          in = tempin.FracEmittance().GetValue();
          out = tempout.FracEmittance().GetValue();
        } else {
	  in = inbeam.FracEmittance().GetValue();
	  out = outbeam.FracEmittance().GetValue();
        }
        ine = in/sqrt(tempin.Size());
        oute = out/sqrt(tempin.Size());
      } else if ( fracpar == "quan" ) {
        in = Math::Quantile(inbeam.Amplitudes(), tempfrac);
        ine = in/sqrt(tempin.Size());
        out = Math::Quantile(outbeam.Amplitudes(), tempfrac);
        oute = out/sqrt(tempin.Size());
      } else {
	Pitch::print(Pitch::warning, "Fractional quantity not recognized, abort");
	break;
      }

      if ( fracpar != "vol" ) {
        value = 100*(out-in)/in;
        error = fabs(value)*sqrt(pow(ine/in, 2)+pow(oute/out, 2));
      } else {
        value = 100*(sqrt(1.+(out-in)/in)-1.);
        error = fabs(value)*sqrt(pow(ine/in, 2)+pow(oute/out, 2));
      }
      fracgraphs[fracpar]->SetPoint((int)((tempfrac+0.01-minfrac)/fracstep), tempfrac, value);
      fracgraphs[fracpar]->SetPointError((int)((tempfrac+0.01-minfrac)/fracstep), 0, error);
    }

    // Display the progress in %
    pbar.GetProgress((tempfrac+0.01-minfrac)/fracstep-1, (maxfrac+0.01-minfrac)/fracstep);
  }

  id = 0;
  TMultiGraph* mgfrac = new TMultiGraph();
  fracgraphs["acut"]->SetTitle("Amplitude, A_{#alpha}"); 
  fracgraphs["subeps"]->SetTitle("Subemittance, e_{#alpha}"); 
  fracgraphs["vol"]->SetTitle("Fractional emittance, #epsilon_{#alpha}"); 
  TLegend *legfrac = new TLegend(.5, .65, .89, .89);
  legfrac->SetLineColorAlpha(0, 0);
  for (const std::string& fracpar : fracpars ) {
    // Set the graph parameters
    fracgraphs[fracpar]->GetXaxis()->SetTitle("Fraction #alpha");
    std::string y_label =
	"#Delta"+par_labels[fracpar]+"_{#alpha}/"+par_labels[fracpar]+"_{#alpha}^{i} [%]";
    fracgraphs[fracpar]->GetYaxis()->SetTitle(y_label.c_str());
    fracgraphs[fracpar]->SetMarkerSize(2);
    fracgraphs[fracpar]->SetMarkerStyle(24+id);
    fracgraphs[fracpar]->SetLineColor(colors[id]);
    fracgraphs[fracpar]->SetLineWidth(3);
    legfrac->AddEntry(fracgraphs[fracpar], fracgraphs[fracpar]->GetTitle(), "PE");
    id++;

    mgfrac->Add(fracgraphs[fracpar], "PE");
  }

  TCanvas *canvfrac = new TCanvas("c", "c", 1200, 800);
  mgfrac->Draw("A");
  mgfrac->SetMinimum(-8);
  mgfrac->SetMaximum(-3);
  mgfrac->SetTitle(";Fraction #alpha;#Delta#epsilon_{#perp}  /#epsilon_{#perp}^{i}  [%]");
  double truechange = fbeta[6]->Eval(globals["toy_beta"]);
  TLine *ltrue = new TLine(0.06, truechange, 0.99, truechange);
  ltrue->Draw("SAME");
  legfrac->Draw("SAME");
  canvfrac->SaveAs("super_frac_toy.pdf");
  delete mgfrac;
  delete canvfrac;

  /////////////////////////////////////////////
  ////////// VARIABLE TRANSMISSION ////////////
  /////////////////////////////////////////////
  Pitch::print(Pitch::info, "Filling the variable transmission graphs");
  std::map<std::string, TGraphErrors*> transgraphs;
  std::vector<std::string> transpars = {"eps"};  
  for (const std::string par : transpars)
      transgraphs[par] = new TGraphErrors();

  // First get the array of amplitudes from the output sample and sort them
  // Should be implemented in the transport object. Select an aperture (TODO TODO TODO)
  std::vector<std::pair<size_t, double>> outamps_id;
//  for (size_t i = 0; i < Size; i++)
//	outamps_id.push_back(std::pair<size_t, double>(i, outbeam.Amplitude(i)));
  for (size_t i = 0; i < Size; i++)
	outamps_id.push_back(std::pair<size_t, double>(i, outbeam.Radius(i)));

  std::sort(outamps_id.begin(), outamps_id.end(),
	    [] (const std::pair<size_t, double>& a, const std::pair<size_t, double>& b) {
		return a.second < b.second;
	    });   

  std::map<std::string, std::vector<double>> tempiSize, tempoutsamples;
  size_t tempn(0);
  double maxtrans(1.), mintrans(.1), transstep(.05);
  double temptrans;
  pbar = ProgressBar();
  for (temptrans = mintrans; temptrans < maxtrans+0.001; temptrans += transstep) {
    // Truncate the in and out samples according to out order
    for (size_t j = tempn; j < (size_t)(Size*temptrans); j++)
      for (const std::string& var : vars) {
	tempiSize[var].push_back(iSize[var][outamps_id[j].first]);
	tempoutsamples[var].push_back(outsamples[var][outamps_id[j].first]);
      }
    tempn = Size*temptrans;

    // Initialize the observed beams
    tempin = Beam::Bunch(tempiSize);
    tempout = Beam::Bunch(tempoutsamples);   

    // Compare and fill the beta plots
    for (const std::string& par : transpars )  {
      if ( par == "eps" ) {
        in = tempin.NormEmittance().GetValue();
	ine = tempin.NormEmittance().GetError();
        out = tempout.NormEmittance().GetValue();
	oute = tempout.NormEmittance().GetError();
      } else {
	std::cerr << "Quantity " << par << " not recognized, abort" << std::endl;
	break;
      }

      value = 100*(out-in)/in;
//      error = fabs(value)*sqrt(pow(ine/in, 2)+pow(oute/out, 2));
      error = 100*sqrt((pow(1+value/100, 2)-0.88*0.88)/(2*temptrans*1e5));
      transgraphs[par]->SetPoint((int)((temptrans+0.001-mintrans)/transstep), 100*temptrans, value);
      transgraphs[par]->SetPointError((int)((temptrans+0.001-mintrans)/transstep), 0., error);
    }

    // Display the progress in %
    pbar.GetProgress((temptrans+1e-3-mintrans)/transstep-1, (maxtrans+1e-3-mintrans)/transstep);
  }

  // Set a legend for the plots as a function of transmission
  TLegend *trans_leg = new TLegend(.65, .11, .89, .3);
  trans_leg->SetLineColorAlpha(0, 0);
  trans_leg->SetFillColor(0);
  trans_leg->AddEntry(transgraphs[transpars[0]],
		      TString::Format("p_{i}: %d MeV/c", (int)globals["toy_mom"]), "");
  trans_leg->AddEntry(transgraphs[transpars[0]],
		      TString::Format("#epsilon_{i}: %d mm", (int)globals["toy_eps"]), "");
  trans_leg->AddEntry(transgraphs[transpars[0]],
		      TString::Format("#beta_{#perp}: %d mm", (int)globals["toy_beta"]), "");

  for (const std::string& par : transpars ) {
    // Set the graph parameters
//    transgraphs[par]->SetTitle(std::string(par_titles[par]+" change").c_str());
    transgraphs[par]->SetTitle("");
    transgraphs[par]->GetXaxis()->SetTitle("T [%]");
    std::string y_label =
	"#Delta"+par_labels[par]+"_{#alpha}/"+par_labels[par]+"_{#alpha}^{in} [%]";
    transgraphs[par]->GetYaxis()->SetTitle(y_label.c_str());
    transgraphs[par]->SetMarkerSize(1.5);
    transgraphs[par]->SetMarkerStyle(21);
    transgraphs[par]->SetLineColor(2);
    transgraphs[par]->SetLineWidth(3);

    TFile ofile(TString::Format("trans_rad_sol_%d.root", (int)globals["toy_drift"]), "RECREATE");
    transgraphs[par]->Write();
    ofile.Close();

    TCanvas *canv = new TCanvas("c", "c", 1200, 800);
    gPad->SetGridx();
    gPad->SetGridy();
    transgraphs[par]->Draw("APE");
//    trans_leg->Draw("SAME");
    canv->SaveAs(std::string(par+"_trans_toy.pdf").c_str());
    delete transgraphs[par];
    delete canv;
  }
  delete trans_leg;

  /////////////////////////////////////////////
  ///////// VARIABLE BETA FUNCTION ////////////
  /////////////////////////////////////////////
  Pitch::print(Pitch::info, "Filling the variable beta graphs");
  // Change the beta function and check the change as a function of it
  std::map<std::string, TGraphErrors*> betagraphs;
  std::vector<std::string> betapars = {"eps"/*, "subeps", "vol"*/};
  for (const std::string par : betapars)
      betagraphs[par] = new TGraphErrors();
  double maxbeta(1000), minbeta(100), betastep(100);
  double tempbeta;
  pbar = ProgressBar();
  for (tempbeta = minbeta; tempbeta < maxbeta+1; tempbeta += betastep) {

    // Create a new input beam
    iSize.clear();
    iSize = GaussianBunch(globals["toy_mass"], globals["toy_mom"], Size,
		     	     globals["toy_eps"], tempbeta, globals["toy_alpha"]);
    Beam::Bunch ibeam(iSize, 0., "in"); 

    // Diffuse and scatter to create the output beam
    outsamples.clear();
    outsamples = iSize;
    if ( globals["toy_scat"] )
        scat.ScatterBunch(outsamples, globals["toy_mom"], globals["toy_mass"]);
    if ( globals["toy_eloss"] )
        eloss.IonizeBunch(outsamples, globals["toy_mom"], globals["toy_mass"]);
    Beam::Bunch obeam(outsamples, 10*depth, "out");

    transmission = (double)obeam.Size()/ibeam.Size();
    ibeam.SetCoreFraction(frac);
    obeam.SetCoreFraction(frac/transmission);   

    // Compare and fill the beta plots
    for (const std::string& par : betapars )  {
      if ( par == "eps" ) {
        in = ibeam.NormEmittance().GetValue();
        out = obeam.NormEmittance().GetValue();
      } else if ( par == "subeps" ) {
        in = tempin.NormEmittance().GetValue();
        out = tempout.NormEmittance().GetValue();
//        in = tempin.Mean("amp").GetValue();
//        out = tempout.Mean("amp").GetValue();
      } else if ( par == "acut" ) {
        in = Math::Max(tempin.Amplitudes());
        out = Math::Max(tempout.Amplitudes());
      } else if ( par == "vol" ) {
	if ( !globals["de"] ) {
          in = tempin.FracEmittance().GetValue();
          out = tempout.FracEmittance().GetValue();
        } else {
	  in = inbeam.FracEmittance().GetValue();
	  out = outbeam.FracEmittance().GetValue();
        }
        ine = in/sqrt(tempin.Size());
        oute = out/sqrt(tempin.Size());
      } else if ( par == "quan" ) {
        in = Math::Quantile(ibeam.Amplitudes(), tempfrac);
        out = Math::Quantile(obeam.Amplitudes(), tempfrac);
      } else {
	Pitch::print(Pitch::warning, "Fractional quantity not recognized, abort");
	break;
      }
      betagraphs[par]->SetPoint((int)((tempbeta+1-minbeta)/betastep), tempbeta, 100*(out-in)/in);
      betagraphs[par]->SetPointError((int)((tempbeta+1-minbeta)/betastep), 0,
					100*sqrt((pow(1+(out-in)/in, 2)-0.88*0.88)*2./1e5));
    }

    // Display the progress in %
    pbar.GetProgress((tempbeta+1e-3-minbeta)/betastep-1, (maxbeta+1e-3-minbeta)/betastep);
  }


  // Set a legend for the plots as a function of beta perp
  TLegend *beta_leg = new TLegend(.65, .11, .89, .3);
  beta_leg->SetLineColorAlpha(0, 0);
  beta_leg->SetFillColor(0);
  beta_leg->AddEntry(betagraphs[betapars[0]],
		     TString::Format("p_{i}: %d MeV/c", (int)globals["toy_mom"]), "");
  beta_leg->AddEntry(betagraphs[betapars[0]],
		     TString::Format("#epsilon_{i}: %d mm", (int)globals["toy_eps"]), "");

  for (const std::string& par : betapars ) {
    // Set the graph parameters
//    betagraphs[par]->SetTitle(std::string(par_titles[par]+" change").c_str());
    betagraphs[par]->SetTitle("");
    betagraphs[par]->GetXaxis()->SetTitle("#beta_{#perp}  [mm]");
    std::string y_label =
	"#Delta"+par_labels[par]+"_{#alpha}/"+par_labels[par]+"_{#alpha}^{in} [%]";
    betagraphs[par]->GetYaxis()->SetTitle(y_label.c_str());
    betagraphs[par]->SetMarkerStyle(24);
    betagraphs[par]->SetMarkerSize(1.5);
    betagraphs[par]->SetLineColor(2);
    betagraphs[par]->SetLineWidth(3);

    TFile ofile(TString::Format("beta_%s_%dmm.root", par.c_str(), (int)globals["toy_eps"]), "RECREATE");
    betagraphs[par]->Write();
    ofile.Close();

    TCanvas *canv = new TCanvas("c", "c", 1200, 800);
    gPad->SetGridx();
    gPad->SetGridy();
    betagraphs[par]->Draw("APE");
    canv->SaveAs(std::string(par+"_beta_toy.pdf").c_str());
    delete canv;
  }
  delete beta_leg;

  // Make a superimposed plot of all the beta plots
  TCanvas* c = new TCanvas ("c", "c", 1200, 800);
  gPad->SetGridx();
  gPad->SetGridy();
  TLegend* lbeta = new TLegend(.7, .11, .89, .3);
  lbeta->SetLineColorAlpha(0, 0);
  lbeta->SetFillColor(0);
  size_t gid(1);
  TMultiGraph* mgbeta = new TMultiGraph();
  mgbeta->SetTitle(";#beta_{#perp} [mm];#Delta x/x [%]");
  for (const std::string& par : betapars ) {
    // Set the graph parameters
    betagraphs[par]->SetLineColor(gid);
    betagraphs[par]->SetMarkerStyle(19+gid);
    betagraphs[par]->GetYaxis()->SetTitle("#Delta x/x [%]");
    mgbeta->Add(betagraphs[par], "PE");
    lbeta->AddEntry(betagraphs[par], par_names[par].c_str(), "lp");
    gid++;
    if ( gid == 5 ) gid++;
  }
  mgbeta->Draw("A");
  lbeta->Draw("SAME");
  c->SaveAs("beta_super_toy.pdf");
  delete c;
   

  /////////////////////////////////////////////
  ///////// VARIABLE INPUT EMITTANCE //////////
  /////////////////////////////////////////////
  Pitch::print(Pitch::info, "Filling the variable emittance graphs");
  // Change the beta function and check the change as a function of it
  std::map<std::string, TGraph*> epsgraphs;
  std::vector<std::string> epspars = {"eps", "subeps", "vol"};
  for (const std::string par : epspars)
      epsgraphs[par] = new TGraph();
  double maxeps(10), mineps(2), epsstep(1);
  double tempeps;
  pbar = ProgressBar();
  for (tempeps = mineps; tempeps < maxeps+1; tempeps += epsstep) {

    // Create a new input beam
    iSize.clear();
    iSize = GaussianBunch(globals["toy_mass"], globals["toy_mom"], Size,
			     tempeps, globals["toy_beta"], globals["toy_alpha"]);
    Beam::Bunch ibeam(iSize, 0., "in");

    // Diffuse and scatter to create the output beam
    outsamples.clear();
    outsamples = iSize;
    if ( globals["toy_scat"] )
        scat.ScatterBunch(outsamples, globals["toy_mom"], globals["toy_mass"]);
    if ( globals["toy_eloss"] )
        eloss.IonizeBunch(outsamples, globals["toy_mom"], globals["toy_mass"]);
    Beam::Bunch obeam(outsamples, 10*depth, "out");

    transmission = (double)obeam.Size()/ibeam.Size();
    ibeam.SetCoreFraction(frac);
    obeam.SetCoreFraction(frac/transmission);   

    // Compare and fill the eps plots
    for (const std::string& par : epspars )  {
      if ( par == "eps" ) {
        in = ibeam.NormEmittance().GetValue();
        out = obeam.NormEmittance().GetValue();
      } else if ( par == "subeps" ) {
        in = tempin.NormEmittance().GetValue();
        out = tempout.NormEmittance().GetValue();
//        in = tempin.Mean("amp").GetValue();
//        out = tempout.Mean("amp").GetValue();
      } else if ( par == "acut" ) {
        in = Math::Max(tempin.Amplitudes());
        out = Math::Max(tempout.Amplitudes());
      } else if ( par == "vol" ) {
	if ( !globals["de"] ) {
          in = tempin.FracEmittance().GetValue();
          out = tempout.FracEmittance().GetValue();
        } else {
	  in = inbeam.FracEmittance().GetValue();
	  out = outbeam.FracEmittance().GetValue();
        }
        ine = in/sqrt(tempin.Size());
        oute = out/sqrt(tempin.Size());
      } else if ( par == "quan" ) {
        in = Math::Quantile(ibeam.Amplitudes(), tempfrac);
        out = Math::Quantile(obeam.Amplitudes(), tempfrac);
      } else {
	Pitch::print(Pitch::warning, "Fractional quantity not recognized, abort");
	break;
      }
      epsgraphs[par]->SetPoint((int)((tempeps+0.01-mineps)/epsstep), tempeps, 100*(out-in)/in);
    }

    // Display the progress in %
    pbar.GetProgress((tempeps+1e-3-mineps)/epsstep-1, (maxeps+1e-3-mineps)/epsstep);
  }

  // Set a legend for the plots as a function of ineps
  TLegend *eps_leg = new TLegend(.65, .7, .89, .89);
  eps_leg->SetLineColorAlpha(0, 0);
  eps_leg->SetFillColor(0);
  eps_leg->AddEntry(epsgraphs[epspars[0]],
		     TString::Format("p_{i}: %d MeV/c", (int)globals["toy_mom"]), "");
  eps_leg->AddEntry(epsgraphs[epspars[0]],
		     TString::Format("#beta_{#perp}: %d mm", (int)globals["toy_beta"]), "");

  for (const std::string& par : epspars ) {
    // Set the graph parameters
    epsgraphs[par]->SetTitle(std::string(par_titles[par]+" change").c_str());
    epsgraphs[par]->GetXaxis()->SetTitle("#epsilon_{n}^{in} [mm]");
    std::string y_label =
	"#Delta"+par_labels[par]+"_{#alpha}/"+par_labels[par]+"_{#alpha}^{in} [%]";
    epsgraphs[par]->GetYaxis()->SetTitle(y_label.c_str());
    epsgraphs[par]->SetMarkerStyle(21);
    epsgraphs[par]->SetLineColor(2);
    epsgraphs[par]->SetLineWidth(3);

    TCanvas *canv = new TCanvas("c", "c", 1200, 800);
    gPad->SetGridx();
    gPad->SetGridy();
    epsgraphs[par]->Draw("APL");
    eps_leg->Draw("SAME");
    canv->SaveAs(std::string(par+"_eps_toy.pdf").c_str());
    delete canv;
  }
  delete eps_leg;

  // Make a superimposed plot of all the eps plots
  c = new TCanvas ("c", "c", 1200, 800);
  gPad->SetGridx();
  gPad->SetGridy();
  TLegend* leps = new TLegend(.7, .7, .89, .89);
  leps->SetLineColorAlpha(0, 0);
  leps->SetFillColor(0);
  gid = 1;
  TMultiGraph* mgeps = new TMultiGraph();
  mgeps->SetTitle(";#epsilon_{n}^{in} [mm];#Delta x/x [%]");
  for (const std::string& par : epspars ) {
    // Set the graph parameters
    epsgraphs[par]->SetLineColor(gid);
    epsgraphs[par]->SetMarkerStyle(19+gid);
    epsgraphs[par]->GetYaxis()->SetTitle("#Delta x/x [%]");
    mgeps->Add(epsgraphs[par], "pl");
    leps->AddEntry(epsgraphs[par], par_names[par].c_str(), "lp");
    gid++;
    if ( gid == 5 ) gid++;
  }
  mgeps->Draw("A");
  leps->Draw("SAME");
  c->SaveAs("eps_super_toy.pdf");
  delete c;

}
