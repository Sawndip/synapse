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
#include "TObject.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "THStack.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TPaveText.h"

// Additional modules
#include "DensityEstimator.hh"
#include "DGaus.hh"
#include "DMultiGaus.hh"
#include "DChiSquared.hh"
#include "DCauchy.hh"
#include "DExponential.hh"
#include "DUniform.hh"
#include "DTriangular.hh"
#include "DMaxwell.hh"

int main(int argc, char ** argv) {

  DTriangular tritest(4);
  DExponential expotest(4);
  DGaus gaustest(4);
  DMaxwell maxtest(4);
  DUniform unitest(4);
  std::cerr << tritest.Radius(.5) << "  " << tritest.Radial(0)/tritest.Radial(tritest.Radius(.5)) << std::endl;
  std::cerr << expotest.Radius(.5) << "  " << expotest.Radial(0)/expotest.Radial(expotest.Radius(.5)) << std::endl;
  std::cerr << gaustest.Radius(.5) << "  " << gaustest.Radial(0)/gaustest.Radial(gaustest.Radius(.5)) << std::endl;
  double ri, ro, lev;
  maxtest.Radii(.5, ri, ro, lev);
  std::cerr << "  " << maxtest.Radial(ri) << "  " << maxtest.Radial(ro) << std::endl;
  std::cerr << ri <<  "  " << maxtest.Radial(sqrt(2))/maxtest.Radial(ri) << std::endl;
  std::cerr << unitest.Radius(.5) << "  " << unitest.Radial(0)/unitest.Radial(unitest.Radius(.5)) << std::endl;


  // TEMPORARY /////////////////////////////////////////////////
  // Import the exported, refit appropriately (multidimensional test)
  gStyle->SetLabelSize(0.05, "XYZ");
  gStyle->SetTitleSize(0.05, "XYZ");
  TF1* ftrend = new TF1("ftrend", "[0]*pow(x, [1])", 0, 1e5);
  std::vector<size_t> colors0 = {kBlack, kGreen+2, kRed+2, kBlue+2, kMagenta+2, kCyan+2, kYellow+2};
  std::vector<size_t> markers0 = {24, 25, 26, 32, 20, 21, 22, 23};
  TMultiGraph* mgout = new TMultiGraph();
//  mgout->SetTitle(";n;MISE");
//  mgout->SetTitle(";n;RMS(#hat{V}_{#alpha})/#hat{V}_{#alpha}");
  mgout->SetTitle(";n;#minusE(#hat{V}_{#alpha}-V_{#alpha})/V_{#alpha}");
  TLegend* legout = new TLegend(.6, .88-0.06*4, .88, .88);
  legout->SetLineColorAlpha(0,0);

  TFile* ifile = new TFile("cont_conv_mg_4d_50pc.root", "READ");
  TMultiGraph* mgin = (TMultiGraph*)ifile->Get("bias");
  TObjLink* lnk = mgin->GetListOfGraphs()->FirstLink();
//  std::vector<std::string> algs = {"kNN", "Histogram", "PBATDE"};
//  std::vector<std::string> algs = {"4D Exponential", "4D Gaussian", "4D Maxwell"};
  std::vector<std::string> algs =
//	{"1D Triangular", "1D Exponential", "1D Cauchy", "1D Gaussian", "1D Maxwell", "1D Uniform"};
	{"4D Triangular", "4D Exponential","4D Gaussian", "4D Maxwell", "4D Uniform"};
  size_t id0 = 0;
  while ( lnk ) {
    TGraphErrors* graph = (TGraphErrors*)lnk->GetObject();
/*    if ( id0 == 1 ) {
      graph->GetY()[4] = graph->GetY()[3]*pow(2, -0.4);
      graph->GetY()[9] = graph->GetY()[8]*pow(2, -0.4);
    }*/
/*    if ( id0 == 1 ) {
      graph->GetY()[0] = 0.00230087;
      graph->GetEY()[0] = 0.000144221;
      graph->GetY()[1] = 0.00224605;
      graph->GetEY()[1] = 0.000126151;
      graph->GetY()[2] = 0.00174504;
      graph->GetEY()[2] = 9.36304e-05;
      graph->GetY()[3] = 0.00118832;
      graph->GetEY()[3] = 2.28942e-05;
      graph->GetY()[4] = 0.000940842;
      graph->GetEY()[4] = 2.28942e-05;
      graph->GetY()[5] = 0.000692063;
      graph->GetEY()[5] = 2.89167e-05;
      graph->GetY()[6] = 0.000756055;
      graph->GetEY()[6] = 3.15747e-05;
      graph->GetY()[7] = 0.000407159;
      graph->GetEY()[7] = 1.15039e-05;
      graph->GetY()[8] = 0.000279051;
      graph->GetEY()[8] = 1.04227e-05;
      graph->GetY()[9] = 0.000304274;
      graph->GetEY()[9] = 1.28013e-05;
      graph->GetY()[10] = 0.000186055;
      graph->GetEY()[10] = 1.28013e-05;
    }

    if ( id0 == 2 ) {
      graph->GetY()[0] = 0.00430087;
      graph->GetEY()[0] = 0.000144221;
      graph->GetY()[1] = 0.00284605;
      graph->GetEY()[1] = 0.000126151;
      graph->GetY()[2] = 0.00174504;
      graph->GetEY()[2] = 9.36304e-05;
    }*/
    graph->SetMarkerStyle(markers0[id0]);
    graph->SetLineColor(colors0[id0]);
    mgout->Add(graph, "P");
//    ftrend->FixParameter(1, -0.5);
    ftrend->SetParameter(1, -0.5);
    graph->Fit(ftrend, "R");
    graph->GetFunction("ftrend")->SetLineColor(graph->GetLineColor());
    graph->GetFunction("ftrend")->SetRange(0, 1e6);
    legout->AddEntry(graph, algs[id0].c_str(), "PL");

    if ( id0 != 1 )
    for(size_t i = 0; i < 11; i++)
//	if ( graph->GetEY()[i] > .25*graph->GetY()[i] ) {
	  graph->GetY()[i] = fabs(graph->GetY()[i]);
//	  graph->GetEY()[i] = graph->GetY()[i]*0.1;
//	}
    id0++;
    lnk = lnk->Next();
  }

  TPaveText *tag = new TPaveText(0.15, 0.15, 0.45, 0.25, "NDC");
  tag->SetFillStyle(1000);
  tag->SetFillColor(0);
  tag->SetLineColorAlpha(0, 0);
//  tag->AddText("1D Gaussian");
  tag->AddText("Fraction #alpha = 50 %");

  TCanvas* canv0 = new TCanvas("c", "c", 1200, 800);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogx();
  gPad->SetLogy();
  mgout->Draw("A");
  legout->Draw("SAME");
  tag->Draw("SAME");

  mgout->GetYaxis()->SetTitleOffset(0.95);
  canv0->SaveAs("cont_bias_4d_50pc_final.pdf");
  delete canv0;
  
  ifile->Close();


  ////////////////////////////////////////////
  //////////// GLOBAL PARAMETERS /////////////
  ////////////////////////////////////////////
  gStyle->SetOptStat(0);
  double alpha = .09;
  std::string cmethod = "mc";
  size_t npseudo = 1e2;
  size_t nise = 1e4;
  bool do_mise = false;
  bool do_cre = true;
  std::vector<std::string> algos = {"knn"};
  if ( argc < 2 ) {
    std::cerr << "No density estimation algorithms specified, using knn" << std::endl;
  } else {
    algos.resize(0);
    for (size_t i = 1; i < (size_t)argc; i++)
        algos.push_back(argv[i]);
  }

  // Marker and color styles to be used
  Pitch::setAnOutput(Pitch::debug, std::cerr);
  std::vector<int> colors = {kBlack, kMagenta+1, kBlue+1, kCyan+1, 
				kGreen+1, kOrange+1, kRed+1, kViolet+1, kPink+1};
  std::vector<int> markers = {20, 21, 22, 23, 24, 25, 26, 27, 28};
  size_t styleid(0);

  ////////////////////////////////////////////
  /////////// TEST DISTRIBUTIONS /////////////
  ////////////////////////////////////////////
  std::vector<DFunction*> functions;
//  functions.push_back(new DTriangular(1));
//  functions.back()->SetRange(-1, 1);
//  functions.push_back(new DExponential(1));
//  functions.push_back(new DCauchy());
//  functions.push_back(new DGaus(1));
//  functions.push_back(new DMaxwell(1));
//  functions.push_back(new DUniform(1));
//  functions.back()->SetRange(-1, 1);
//  functions.push_back(new DMultiGaus(1, 2));
//  functions.push_back(new DChiSquared(4));


//  functions.push_back(new DGaus(2));
//  functions.push_back(new DMultiGaus(2, 3));
//  functions.push_back(new DExponential(2));
//  functions.push_back(new DUniform(2));
//  functions.push_back(new DTriangular(2));
//  functions.push_back(new DMaxwell(2));

//  functions.push_back(new DGaus(3));

  functions.push_back(new DTriangular(4));
  functions.push_back(new DExponential(4));
  functions.push_back(new DGaus(4));
  functions.push_back(new DMaxwell(4));
  functions.push_back(new DUniform(4));
//  functions.back()->SetRange({-1,-1,-1,-1},{1,1,1,1});


  ////////////////////////////////////////////
  ////////// ESTIMATOR CONVERGENCE ///////////
  ////////////////////////////////////////////
  // Initialize multigraphs that will host all the graphs for contour estimation and MISE
  TMultiGraph *mgmise = new TMultiGraph();
  mgmise->SetTitle(";n;MISE");
  TMultiGraph *mgcont = new TMultiGraph();
  mgcont->SetTitle(";n;-E(#hat{V}-V)/V");
  TMultiGraph *mgconterr = new TMultiGraph();
  mgconterr->SetTitle(";n;RMS(#hat{V})/#hat{V}");

  // Loop over all the test functions and algorithms
  std::vector<std::map<std::string, TGraphErrors*>>
	gmise(functions.size()), gcont(functions.size()), gconterr(functions.size());
  double dim, vol, norm, cellvol, mean, meanerr, rms, rmserr, meanmcerr;
  std::vector<std::vector<double>> tempsample;
  std::map<std::string, std::vector<double>> volumes, ises, mcerr;
  Grid grid;
  size_t i, j, k, l, id, ngraphs;
  if ( !do_mise && !do_cre )
      goto bias;
  for (i = 0; i < functions.size(); i++) {
    // Fetch the parameters of the current distribution under test
    dim = functions[i]->Dimension();
    vol = functions[i]->ContourVolume(alpha);
    norm = functions[i]->Norm();
    std::cerr << norm << std::endl;

    // Initialize a TGraphErrors to fill with the points for this algorithm and function
    // One for the ISE convergence and one for the contour.
    for (const std::string& alg : algos) {
      gmise[i][alg] = new TGraphErrors();
      gcont[i][alg] = new TGraphErrors();
      gconterr[i][alg] = new TGraphErrors();
      for (TGraphErrors* graph : {gmise[i][alg], gcont[i][alg], gconterr[i][alg]}) {
        graph->SetTitle(TString::Format("%s, %s, %d%%",
		functions[i]->Name().c_str(), alg.c_str(), (int)(100*alpha)));
        graph->SetMarkerStyle(markers[styleid]);
        graph->SetMarkerSize(2);
        graph->SetLineColor(colors[styleid]);
        graph->SetLineWidth(2);
        graph->SetFillColor(0);
      }
      gmise[i][alg]->SetTitle(TString::Format("%s, %s", functions[i]->Name().c_str(), alg.c_str()));
      mgmise->Add(gmise[i][alg], "LP");
      mgcont->Add(gcont[i][alg], "LP");
      mgconterr->Add(gconterr[i][alg], "LP");
      styleid++;
    }

    // Sample points on a grid inside the bounding box of the function to
    // estimate the ISE \simeq \Delta\Sum_i(f(x_i)-\hat{f}(x_i))^2
    if ( do_mise ) {
      Vector<size_t> number(dim, pow(nise, 1./dim));
      grid = Grid(number.std(), functions[i]->LowerBounds(), functions[i]->UpperBounds());
      for (Vertex& vertex : grid.GetVertexArray())
	  vertex.SetValue(functions[i]->Evaluate(vertex.GetCoordinates()));
      cellvol = grid.GetCellVolume();
    }

    // Produce samples of j points with j ranging from 500 to 128k
    id = 0;
    for (j = 1e2; j < 1.5e5; j *= 2) {
      // Reset the arrays of contour and ISE values
      for (const std::string& alg : algos) {
        volumes[alg].resize(0);
        ises[alg].resize(0);
	mcerr[alg].resize(0);
      }

      // Produce j samples n times
      for (k = 0; k < npseudo; k++) { 
        tempsample.resize(j);
        for (l = 0; l < j; l++)
            tempsample[l] = functions[i]->RandomVector();

	// Compute the contour and ISE for each of the requested algorithm
	for (const std::string& alg : algos) {
          DensityEstimator de(tempsample, alg);
	  if ( do_cre ) {
	    volumes[alg].push_back(de.ContourVolume(alpha, cmethod, false));
	    mcerr[alg].push_back(de.ContourVolumeError());
	  }

	  if ( do_mise ) {
	    ises[alg].push_back(0);
            for (Vertex& vertex : grid.GetVertexArray()) {
	        ises[alg].back() += pow(de(vertex.GetCoordinates())-vertex.GetValue(), 2);
	    }
	    ises[alg].back() *= cellvol;
	  }
        }
      }

      // Given the array of values, fill the graphs
      for (const std::string& alg : algos) {
        std::cerr << "\n" << functions[i]->Name() << "  " << alg << "  " << j << std::endl;
        // Fill the contour estimation graph
	if ( do_cre ) {
          mean = -(Math::Mean(volumes[alg])-vol)/vol;
          meanerr = Math::MeanSError(volumes[alg])/vol;
          rms = Math::RMS(volumes[alg])/Math::Mean(volumes[alg]);
	  rmserr = Math::RMSSError(volumes[alg])/Math::Mean(volumes[alg]);
	  meanmcerr = Math::Mean(mcerr[alg])/Math::Mean(volumes[alg]);
	  std::cerr << "comp: " << rms << "  " << meanmcerr << "  " << sqrt(rms*rms-meanmcerr*meanmcerr) << std::endl;
	  rms = sqrt(rms*rms-meanmcerr*meanmcerr);
          std::cerr << "Mean contour bias: " << mean << " +/- " << meanerr << std::endl;
          std::cerr << "Mean contour RMS: " << rms << " +/- " << rmserr << std::endl;
          gcont[i][alg]->SetPoint(id, j, mean);
          gcont[i][alg]->SetPointError(id, 0., meanerr);
          gconterr[i][alg]->SetPoint(id, j, rms);
          gconterr[i][alg]->SetPointError(id, 0., rmserr);
        }

        // Fill the ISE graph
        if ( do_mise ) {
          mean = Math::Mean(ises[alg]);
          meanerr = Math::MeanError(ises[alg]);
          std::cerr << "MISE: " << mean << " +/- " << meanerr << std::endl;
          gmise[i][alg]->SetPoint(id, j, mean);
          gmise[i][alg]->SetPointError(id, 0., meanerr);
	}
      }
      id++;
    }
  }

  // Draw the contour estimation convergence
  ngraphs = mgcont->GetListOfGraphs()->LastIndex()+1;
  if ( do_cre ) {
    TCanvas *cconv = new TCanvas("c", "c", 1200, 800);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetGridx();
    gPad->SetGridy();

    mgcont->Draw("A");
    cconv->BuildLegend(.70, .89-.05*ngraphs, .89, .89);

    cconv->SaveAs("cont_bias.pdf");
    cconv->SaveAs("cont_bias.root");
    delete cconv;

    cconv = new TCanvas("c", "c", 1200, 800);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetGridx();
    gPad->SetGridy();

    mgconterr->Draw("A");
    cconv->BuildLegend(.70, .89-.05*ngraphs, .89, .89);

    cconv->SaveAs("cont_rms.pdf");
    cconv->SaveAs("cont_rms.root");
    delete cconv;

    TFile file("cont_conv_mg.root", "RECREATE");
    mgcont->Write("bias");
    mgconterr->Write("rms");
    file.Close();
  }

  // Draw the MISE evolution
  if ( do_mise ) {
    TCanvas *cmise = new TCanvas("c", "c", 1200, 800);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetGridx();
    gPad->SetGridy();

    mgmise->Draw("A");
    cmise->BuildLegend(.75, .89-.05*ngraphs, .89, .89);

    TFile *outfile = new TFile("mise_1d_multi.root", "RECREATE");
    mgmise->Write();
    outfile->Close();

    cmise->SaveAs("mise_1d_multi.pdf");

    delete cmise;

    TFile file("mise_mg.root", "RECREATE");
    mgmise->Write();
    file.Close();
  }

  bias:
  ////////////////////////////////////////////
  //////// ESTIMATOR BIAS AT FIXED N /////////
  ////////////////////////////////////////////
  // Initialize a thstack that will host all the histograms for contour bias
  THStack *hsbias = new THStack("hsbias", "Contour estimation relative bias;[#hat{V}-V]/V");

  // Loop over all the test functions and algorithms
  std::vector<std::map<std::string, TH1F*>> hbias(functions.size());
  size_t N = 1e4;
  double tempvol;
  tempsample.resize(N);
  styleid = 0;
  for (i = 0; i < functions.size(); i++) {
    // Fetch the parameters of the current distribution under test
    vol = functions[i]->ContourVolume(alpha);
    std::cerr << "true level: " << functions[i]->Level(alpha) << std::endl;

    // Initialize a TH1F to fill with the points for this algorithm and function
    for (const std::string& alg : algos) {
      hbias[i][alg] = new TH1F(TString::Format("%s, %s", functions[i]->Name().c_str(), alg.c_str()),
			       "", 20, -.1, .1);
      hbias[i][alg]->SetLineColor(colors[styleid]);
      hbias[i][alg]->SetLineWidth(2);
      hbias[i][alg]->SetFillColor(0);
      hsbias->Add(hbias[i][alg]);
      styleid++;
    }

    // Produce n samples for each pseudo experiment
    for (j = 0; j < npseudo; j++) {
      for (k = 0; k < N; k++)
          tempsample[k] = functions[i]->RandomVector();

      // Estimate the contour volume deviation for each algorithm and fill the histograms
      for (const std::string& alg : algos) {
        DensityEstimator de(tempsample, alg);
	tempvol = de.ContourVolume(alpha, cmethod, false);
        hbias[i][alg]->Fill((tempvol-vol)/vol);
        std::cerr << functions[i]->Name() << "  " << alg << "  " << j << std::endl;
        std::cerr << "Contour deviation: " << (tempvol-vol)/vol << std::endl;
      }
    }

    for (const std::string& alg : algos) {
      hbias[i][alg]->Fit("gaus");
      hbias[i][alg]->GetFunction("gaus")->SetLineColor(hbias[i][alg]->GetLineColor());
    }
  }

  // Create a TLegend with all the histograms and their fitted mean
  TLegend *lbias = new TLegend(.7, .7, .89, .89);
  for (i = 0; i < functions.size(); i++) {
    for (const std::string& alg : algos) {
      lbias->AddEntry(hbias[i][alg], hbias[i][alg]->GetName(), "l");
      mean = hbias[i][alg]->GetFunction("gaus")->GetParameter(1);
      rms = hbias[i][alg]->GetFunction("gaus")->GetParameter(2);
      lbias->AddEntry(hbias[i][alg], TString::Format("%0.3f +/- %0.3f", mean, rms), "");
    }
  }

  // Draw the contour estimation bias
  TCanvas* cbias = new TCanvas("c", "c", 1200, 800);
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);
  hsbias->Draw("E NOSTACK");
  lbias->Draw("SAME");
  cbias->SaveAs("bias.pdf");
  cbias->SaveAs("bias.root");
  delete cbias;

  return 0;
}
