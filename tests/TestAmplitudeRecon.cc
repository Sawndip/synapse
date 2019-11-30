// Cpp includes
#include <ctime>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <regex>
#include <map>

// Root includes
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH1.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLatex.h"

// Additional modules
#include "Globals.hh"
#include "Generator.hh"
#include "Bunch.hh"
#include "ScatterGraph.hh"

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
  Matrix<double> invcovmat = covmat.Inverse();
  eps = pow(covmat.Determinant(), .25)/105.66;
  Matrix<double> x(4, 1), xt(1, 4);
  for (size_t i = 0; i < amps.size(); i++) {
    x[0][0] = beam.Samples("px")[i];
    x[1][0] = beam.Samples("py")[i];
    x[2][0] = beam.Samples("x")[i];
    x[3][0] = beam.Samples("y")[i];
    xt = x.Transpose();
    amps[i] = eps*(xt*invcovmat*x)[0][0];
  }

  return amps;
}

std::map<std::string, ScatterGraph*> ScatterSections(const std::string& type,
						     Beam::Bunch& beam) {

  std::vector<std::string> vars = {"x", "y", "px", "py", "pz"};
  std::map<std::string, ScatterGraph*> splots;
  std::string name;
  size_t i, j;
  for (i = 0; i < vars.size(); i++) {
    for (j = i+1; j < vars.size(); j++) {
      name = vars[i]+vars[j];
      beam.SetName(type);
      splots[name] = beam.AmplitudeScatter(vars[i], vars[j]);
    }
  }

  return splots;
}

int main(int argc, char ** argv) {

  ////////////////////////////////////////////
  //////////// GLOBAL PARAMETERS /////////////
  ////////////////////////////////////////////
  // Reroot the info messages to the terminal for now
  Pitch::setAnOutput(Pitch::info, std::cerr);

  // Toy MC algorithm global parameters (stored in globals)
  Globals &globals = Globals::GetInstance(argc, argv);

  // Set the Stat box to be along the top right corner of the TPad in the TCanvas
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(.05, "XYZ");
  gStyle->SetLabelSize(.05, "XY");
  gStyle->SetTitleOffset(.9, "X");
  gStyle->SetTitleOffset(1, "Y");
  gStyle->SetTitleOffset(.7, "Z");
//  gSystem->Load("libTree");
//  gStyle->SetStatX(.9);
//  gStyle->SetStatY(.9);
//  gStyle->SetOptStat(2210);
//  gStyle->SetOptFit(1);

  // Characteristics of the histograms
  size_t nsamples = (double)globals["toy_n"];
  size_t nbins = 20;
  double minamp(0), maxamp(100);

  ////////////////////////////////////////////
  ////////// INPUT DISTRIBUTIONS /////////////
  ////////////////////////////////////////////
  // Realistic gaussian muon beam
  Generator gen(4, globals["toy_mass"]);
  gen.SetMatrixParametrisation(140, 3, 150, 0);
  gen.GetCovarianceMatrix().Print();

  gen.SetMatrixParametrisation(globals["toy_mom"], globals["toy_eps"],
				globals["toy_beta"], globals["toy_alpha"]);
//  Beam::Bunch inbeam = gen.GaussianBunch(globals["toy_n"]);
  Beam::Bunch inbeam = gen.SpiralBunch(globals["toy_n"], 20e-6);
  
  std::map<std::string, std::vector<double>> cutsamples;

  // Propagate through an aberrant lens
  for (size_t i = 0; i < nsamples; i++) {
/*      if ( pow(samples["x"][i], 2) + pow(samples["y"][i], 2) > 100*100 ) {
        for (const std::string& var : {"x", "y", "px", "py", "pz"} )
	     samples[var].erase(samples[var].begin()+i);
	continue;
      }
*/

//      samples["px"][i] -= samples["x"][i]*(1-1.5e-4*pow(samples["x"][i], 2))/4.;
//      samples["py"][i] -= samples["y"][i]*(1-1.5e-4*pow(samples["y"][i], 2))/4.;


/*      if ( pow(samples["x"][i], 2) + pow(samples["y"][i], 2) > 150*150 ) {
        for (const std::string& var : {"x", "y", "px", "py", "pz"} )
	     samples[var].erase(samples[var].begin()+i);
	continue;
      }
*/
  }

  // Apply an aperture
//  Beam::Bunch inbeam(samples, 0, "intest");
  for (size_t i = 0; i < nsamples; i++) {
//      if ( inbeam.Radius(i) < 150  ) {
        for (const std::string& var : {"x", "y", "px", "py", "pz"} )
	     cutsamples[var].push_back(inbeam.Samples(var)[i]);
//      }
  }


  // Initialize the beam object
  Beam::Bunch beam(cutsamples, 0, "test");
  beam.SetCorrectedAmplitudes();
  beam.SetCoreFraction(0.25);
//  Bunch beam90 = beam.Fraction(.9);
//  Bunch beam80 = beam.Fraction(.8);

  double eps = beam.NormEmittance().GetValue();
  TLatex* teps = new TLatex(140, 150, TString::Format("#varepsilon_{#perp}*  = %0.2f mm", eps));


  ////////////////////////////////////////////
  ///// DISTRIBUTION PROFILE AND ELLIPSES ////
  ////////////////////////////////////////////
  // Draw a TH2F and the relevant matrices on top of it
  TH2F* hxpx = beam.Histogram("x", "px", -200, 200, -99.99, 99.99);
  hxpx->GetXaxis()->SetTitleOffset(0.9);
  hxpx->GetYaxis()->SetTitleOffset(0.9);
  TCanvas* cell = new TCanvas("c", "c", 1200, 800);
  hxpx->SetTitle("");
  hxpx->Draw("");

  beam.Ellipse("x", "px")->Draw("SAME");

//  beam90.Ellipse("x", "px")->Draw("SAME");
//  beam80.Ellipse("x", "px")->Draw("SAME");

  Matrix<double> truemat(4, 4);
  std::map<size_t, size_t> mapping = {{0,1},{1,3},{2,0},{3,2}};
  for (size_t i = 0; i < 4; i++)
    for (size_t j = 0; j < 4; j++)
        truemat[i][j] = gen.GetCovarianceMatrix()[mapping[i]][mapping[j]];
  DGaus true_dist({0, 0, 0, 0}, gen.GetCovarianceMatrix());
  TEllipse *true_ellipse = (TEllipse*)true_dist.Contour2D(0.09)[0];
  true_ellipse->SetLineColor(kGreen+2);
  true_ellipse->Draw("SAME");

  TEllipse *robust_ellipse = beam.RobustEllipse("x", "py");
  robust_ellipse->SetLineColor(kRed);
  robust_ellipse->Draw("SAME");

  Matrix<double> coremat = beam.CoreCovarianceMatrix();
  Matrix<double> corematxpx(2, 2);
  corematxpx[0][0] = coremat[2][2];
  corematxpx[0][1] = coremat[2][0];
  corematxpx[1][0] = coremat[0][2];
  corematxpx[1][1] = coremat[0][0];
  DGaus core_dist({0, 0}, corematxpx);
  TEllipse *core_ellipse = (TEllipse*)core_dist.Contour2D(0.09)[0];
  core_ellipse->SetLineColor(kBlue+2);
  core_ellipse->Draw("SAME");

  teps->Draw("SAME");

  cell->SaveAs("amp_xpx_ellipses.pdf");
  

  ////////////////////////////////////////////
  //////// AMPLITUDE DISTRIBUTIONS ///////////
  ////////////////////////////////////////////
  // Check which methods are requested to be compared, compute the amplitude arrays
  std::map<std::string, std::vector<double>> amps;
  std::map<std::string, std::map<std::string, ScatterGraph*>> scat;
  std::vector<std::string> types;
//  amps["regular"] = beam.Amplitudes();
  amps["regular"] = CustomAmplitudes(beam, truemat);
  scat["regular"] = ScatterSections("regular", beam);
  types.push_back("regular");
  if ( globals["corrected"] ) {
    beam.SetCorrectedAmplitudes();
    amps["corrected"] = beam.Amplitudes();
    scat["corrected"] = ScatterSections("corrected", beam);
    types.push_back("corrected");
  }
  if ( globals["mcd"] ) {
    beam.SetMCDAmplitudes();
//    amps["mcd"] = beam.Amplitudes();
    amps["mcd"] = CustomAmplitudes(beam, beam.CoreCovarianceMatrix());
    scat["mcd"] = ScatterSections("mcd", beam);
    types.push_back("mcd");
  }

  /*TH1F* histr = new TH1F("ampr", "", nbins, minamp, maxamp);
  histr->FillN(amps["regular"].size(), &amps["regular"][0], NULL);
  TH1F* histg = new TH1F("ampg", "", nbins, minamp, maxamp);
  histg->FillN(amps["corrected"].size(), &amps["corrected"][0], NULL);
  double kolmo = histr->KolmogorovTest(histg);
  std::cerr << kolmo << std::endl;*/

  // Draw all the distributions on a same canvas
  TF1* fchi2 = new TF1("fchi2", "[0]*pow(2*[1], -2)*x*exp(-x/2/[1])", 0, 100);
  fchi2->SetLineColor(1);
  fchi2->SetParameters(/*nsamples**/(maxamp-minamp)/nbins, globals["toy_eps"]);
  THStack* hstack = new THStack("stack", ";A_{#perp}  [mm]; Density [mm^{-1}]");
//  hstack->GetYaxis()->SetTitleOffset(1.2);
  TLegend* leg = new TLegend(.65, .89-.075*amps.size(), .89, .89);
  leg->SetLineColor(0);
  leg->SetFillStyle(0);
  size_t id(0);
  std::vector<size_t> colors = {kBlack, kGreen+2, kRed+2, kBlue+2, kMagenta+2, kCyan+2, kYellow+2};
  for (const std::pair<std::string, std::vector<double>>& amp : amps) {
    TH1F* hist;
    if ( !globals["rebin"] ) {
      hist = new TH1F(TString::Format("amp_%s", amp.first.c_str()), "", nbins, minamp, maxamp);
    } else {
      // Get the edges of bins of equal size in R^4 (A_perp is R^2)
      std::vector<double> edges(nbins+1);
      for (size_t i = 0; i < nbins+1; i++)
          edges[i] = minamp + pow((double)i/nbins, .5)*(maxamp-minamp);
      hist = new TH1F(TString::Format("amp_%s", amp.first.c_str()), "", nbins, &edges[0]);        
    }
    hist->SetLineColor(colors[id+1]);
    hist->SetLineWidth(2);
    hist->SetMarkerStyle(24+id);
    hist->SetMarkerSize(1.5);
    hist->FillN(amp.second.size(), &amp.second[0], NULL);
    hist->Sumw2();
    hist->Scale(1./nsamples);
    leg->AddEntry(hist, amp.first.c_str(), "LPE");
    hstack->Add(hist, "LPE");
    ++id;
  }

  TCanvas *c = new TCanvas("c", "c", 1200, 800);
  hstack->Draw("NOSTACK");
  fchi2->Draw("SAME");
  leg->Draw("SAME");
  c->SaveAs("amp_comp.pdf");
  delete c;

  ////////////////////////////////////////////
  //////// AMPLITUDE SCATTER PLOTS ///////////
  ////////////////////////////////////////////
  for (std::pair<std::string, ScatterGraph*> it : scat["regular"]) {
//    TCanvas* cscat = new TCanvas("c", "c", 1800, 1200);
//    cscat->Divide(2, 2);
    TCanvas* cscat = new TCanvas("c", "c", 1200, 800);
//    cscat->Divide(1, 1);
    size_t id(0);
    for (const std::string& type : types) {
      cscat->cd(1);
      gPad->SetLogz();
      gPad->SetRightMargin(.15);
      scat[type][it.first]->Draw();
//      scat[type][it.first]->SetMinimum(.1);
      ++id;
    }

    cscat->SaveAs(("amp_"+it.first+".pdf").c_str());
    delete cscat;
  }
}
