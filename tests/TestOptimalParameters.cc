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
#include "TH2.h"
#include "TF1.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TMarker.h"
#include "TGaxis.h"

// Additional modules
#include "Statistics.hh"
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

  gStyle->SetLabelSize(0.04, "XYZ");
  gStyle->SetTitleSize(0.04, "XYZ");

  std::string algo = "knn";
  if ( argc < 2 ) {
    std::cerr << "No density estimation algorithm specified, using knn" << std::endl;
  } else {
    algo = argv[1];
  }
  std::string var = (algo == "knn") ? "k" : "J";


  std::vector<size_t> colors = {kBlack, kGreen+2, kRed+2, kBlue+2, kMagenta+2, kCyan+2, kYellow+2};
  std::vector<size_t> markers = {24, 25, 26, 32, 20, 21, 22, 23};



  // TEMPORARY ////////////////////////////////////////////////:
/*  TF1* ftrend = new TF1("ftrend", "[0]*pow(x, [1])", 0, 1e5);
  ftrend->SetParameters(1, .5);
  std::vector<std::string> names =
	{"1D Triangular", "1D Exponential", "Cauchy", "1D Gaussian", "1D Maxwell", "1D Uniform"};
  std::vector<double> expos = {.6, .75, .75, .6, .6, .4};
//  std::vector<double> expos = {1./3, 1./3, .2, 1./3, 1./3, .5};

  // Import the exported, refit appropriately
  TFile *ifile = new TFile("pbatde_opt_1D_multi.root", "READ");
  TMultiGraph* mgin = (TMultiGraph*)ifile->Get("");
  TMultiGraph* mgout = new TMultiGraph();
  mgout->SetTitle(";n;J^{*}");
  TLegend* legout = new TLegend(.12, .88-0.06*5, .4, .88);
  legout->SetLineColorAlpha(0, 0);

  TObjLink* lnk = mgin->GetListOfGraphs()->FirstLink();
  size_t id0(0);
  while ( lnk ) {
    if ( id0 != 4 ) {
      TGraphErrors* graph = (TGraphErrors*)lnk->GetObject();
//      if ( id0 == 1 || id0 == 2 ) {
//	graph->GetY()[7] = graph->GetY()[6]*pow(2, .75);
//      }
      
      mgout->Add(graph, "P");
//      ftrend->FixParameter(1, expos[id0]);
      graph->Fit(ftrend, "R");
      graph->GetFunction("ftrend")->SetLineColor(graph->GetLineColor());
      graph->GetFunction("ftrend")->SetRange(0, 1e6);
      legout->AddEntry(graph, names[id0].c_str(), "PL");

      for(size_t i = 0; i < 8; i++)
	if ( graph->GetEY()[i] > .25*graph->GetY()[i] ) {
	  graph->GetY()[i] = ftrend->Eval(graph->GetX()[i]);
	  graph->GetEY()[i] = graph->GetY()[i]*0.1;
	}
    }
    id0++;
    lnk = lnk->Next();
  }

  ifile = new TFile("pbatde_opt_1D_max.root", "READ");
  mgin = (TMultiGraph*)ifile->Get("");
  lnk = mgin->GetListOfGraphs()->FirstLink();
  id0 = 4;
  TGraphErrors* graph = (TGraphErrors*)lnk->GetObject();
  graph->SetMarkerStyle(20);
  graph->SetLineColor(kMagenta+2);
  mgout->Add(graph, "P");
//  ftrend->FixParameter(1, expos[id0]);
  graph->Fit(ftrend, "R");
  graph->GetFunction("ftrend")->SetLineColor(graph->GetLineColor());
  graph->GetFunction("ftrend")->SetRange(0, 1e6);
      for(size_t i = 0; i < 8; i++)
	if ( graph->GetEY()[i] > .25*graph->GetY()[i] ) {
	  graph->GetY()[i] = ftrend->Eval(graph->GetX()[i]);
	  graph->GetEY()[i] = graph->GetY()[i]*0.1;
	}
  legout->AddEntry(graph, names[id0].c_str(), "PL");
  
  TCanvas *canv0 = new TCanvas("c", "c", 1200, 800);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogx();
  gPad->SetLogy();
  mgout->Draw("A");
  legout->Draw("SAME");
  canv0->SaveAs("pbatde_opt_final.pdf");
  delete canv0;

  ifile->Close();

  // Import the exported, refit appropriately (multidimensional test)
  mgout = new TMultiGraph();
  mgout->SetTitle(TString::Format(";n;%s^{*}", var.c_str()));
  legout = new TLegend(.12, .88-0.06*4, .4, .88);
  legout->SetLineColorAlpha(0,0);

  expos = {0.6, 0.4, 0.5, 0.5};
  for (int i = 0; i < 4; i++) {
    ifile = new TFile(TString::Format("%s_opt_%dD.root", algo.c_str(), i+1), "READ");
    mgin = (TMultiGraph*)ifile->Get("");
    lnk = mgin->GetListOfGraphs()->FirstLink();
    graph = (TGraphErrors*)lnk->GetObject();
    graph->SetMarkerStyle(markers[i]);
    graph->SetLineColor(colors[i]);
//    ftrend->FixParameter(1, expos[i]);
    graph->Fit(ftrend, "R");
    graph->GetFunction("ftrend")->SetLineColor(colors[i]);
    graph->GetFunction("ftrend")->SetRange(0, 1e6);
    legout->AddEntry(graph, TString::Format("%dD Gaussian", i+1), "LP");
    mgout->Add(graph, "P");
  }

  canv0 = new TCanvas("c", "c", 1200, 800);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogx();
  gPad->SetLogy();
  mgout->Draw("A");
  legout->Draw("SAME");
  canv0->SaveAs(TString::Format("%s_opt_final_nD.pdf", algo.c_str()));
  delete canv0;
  
  ifile->Close();*/

/*  DExponential exp2(2);
  exp2.Print("CONTZ");    

  std::vector<double> bla1({1, 2, 3}), bla2({2, 4, 9});
  std::cerr << Math::Covariance(bla1, bla2) << std::endl;
//  DGaus(bla1, bla2);


  DExponential expo(4);
  std::cerr << "vol: " << expo.ContourVolume(0.09) << std::endl;

  DMultiGaus multig(1, 2, 3);
  multig.SetRange(1.5, 5);
  multig.PrintRadial();
  std::cerr << multig.Norm() << std::endl;
  std::vector<double> lower, upper;
  multig.Range(lower, upper);
  std::cerr << lower[0] << "  " << upper[0] << std::endl;


  DGaus gauso(1);
  gauso.SetRange(-3, 3);
  std::cerr << gauso.Norm() << std::endl;

  Matrix<double> tmat(1, 2);
  tmat.CofactorMatrix();

  Matrix<double> matt(4, 4, 1);
  matt[0][1] = 5;
  matt.Print();
  std::cerr << matt.IsSymmetric() << std::endl;
  matt.Identity();
  matt.Print();
  std::cerr << matt.IsSymmetric() << std::endl;

  Matrix<double> gaus4_cov({{1, 0, 1, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}});
  std::vector<double> gaus4_means({0, 0, 0, 0});
  DGaus *fgaus4 = new DGaus(gaus4_means, gaus4_cov);
  fgaus4->Norm();*/


  ////////////////////////////////////////////
  //////////// GLOBAL PARAMETERS /////////////
  ////////////////////////////////////////////
  gStyle->SetOptStat(0);
  size_t nsamples = 1e4;
  size_t nise = 1e4;
  size_t npseudo = 10;
  size_t start(20), end(200), step(5);
  size_t npoints = (end+1e-9-start)/step;

  ////////////////////////////////////////////
  /////////// TEST DISTRIBUTIONS /////////////
  ////////////////////////////////////////////
  std::vector<DFunction*> functions;
//  functions.push_back(new DTriangular(1));
//  functions.push_back(new DExponential(1));
//  functions.push_back(new DCauchy());
//  functions.push_back(new DGaus(1));
//  functions.push_back(new DMaxwell(1));
//  functions.push_back(new DUniform(1));

//  functions.push_back(new DMultiGaus(1, 2));
//  functions.push_back(new DChiSquared(4));
//  functions.push_back(new DExponential(1));

  // Limit the range to a certain fraction of the support
/*  double alpha = 0.5;
  std::cerr << "Fraction: " << alpha << std::endl;

  DTriangular tri1(1);
  std::cerr << "Triangular radius: " << tri1.Radius(alpha) << std::endl;
  functions[0]->SetRange(-tri1.Radius(alpha), tri1.Radius(alpha));
  DExponential exp1(1);
  std::cerr << "Exponential radius: " << exp1.Radius(alpha) << std::endl;
  functions[1]->SetRange(-exp1.Radius(alpha), exp1.Radius(alpha));
  DCauchy cauchy;
  std::cerr << "Cauchy radius: " << cauchy.ContourVolume(alpha)/2 << std::endl;
  functions[2]->SetRange(-cauchy.ContourVolume(alpha)/2, cauchy.ContourVolume(alpha)/2);
  DExponential gaus1(1);
  std::cerr << "Gaussian radius: " << gaus1.Radius(alpha) << std::endl;
  functions[1]->SetRange(-gaus1.Radius(alpha), gaus1.Radius(alpha));
  DUniform uni1();
  std::cerr << "Uniform radius: " << alpha << std::endl;
  functions[4]->SetRange(-alpha, alpha);*/


//  functions.push_back(new DGaus(2));
//  functions.push_back(new DMultiGaus(2, 3));
//  functions.push_back(new DExponential(2));
//  functions.push_back(new DUniform(2));
//  functions.push_back(new DTriangular(2));
//  functions.push_back(new DMaxwell(2));

//  functions.push_back(new DGaus(3));

//  functions.push_back(new DTriangular(4));
//  functions.push_back(new DGaus(4));
  functions.push_back(new DUniform(4));



  ////////////////////////////////////////////
  ////////// OPTIMIZE FOR GIVEN N ////////////
  ////////////////////////////////////////////
  // Initialize multigraphs that will host all the optimization graph
  TMultiGraph *mgaic = new TMultiGraph();
  mgaic->SetTitle(TString::Format(";%s;AIC", var.c_str()));
  mgaic->SetName(TString::Format("%s_aic", algo.c_str()));
  TMultiGraph *mgbic = new TMultiGraph();
  mgbic->SetTitle(TString::Format(";%s;BIC", var.c_str()));
  mgbic->SetName(TString::Format("%s_bic", algo.c_str()));
  TMultiGraph *mgise = new TMultiGraph();
  mgise->SetTitle(TString::Format(";%s;ISE", var.c_str()));
  mgise->SetName(TString::Format("%s_ise", algo.c_str()));

  // Test kNN or the PBATDE for different values of k or J
  std::vector<TGraph*> gaic(functions.size()), gbic(functions.size()), gise(functions.size());
  std::vector<TGraph*> maic(functions.size()), mbic(functions.size()), mise(functions.size());
  std::vector<double> aic(npoints), bic(npoints), ise(npoints);
  std::vector<std::vector<double>> tempsample;
  std::string tempalg;
  Grid grid;
  double cellvol, lnL, pen;
  size_t id, minid;
  size_t i, j, k, l;
  for (i = 0; i < functions.size(); i++) {
    // Intialize graphs for the current distribution
    gaic[i] = new TGraph();
    gbic[i] = new TGraph();
    gise[i] = new TGraph();
    for (TGraph* graph : {gaic[i], gbic[i], gise[i]}) {
      graph->SetTitle(TString::Format("%s, %s", functions[i]->Name().c_str(), algo.c_str()));
      graph->SetMarkerStyle(20+i);
      graph->SetLineColor(1+i);
      graph->SetLineWidth(2);
      graph->SetFillColor(0);
    }
  
    mgaic->Add(gaic[i], "PL");
    mgbic->Add(gbic[i], "PL");
    mgise->Add(gise[i], "PL");

    // Sample the requested amount of samples from the current distribution
    tempsample.resize(0);
    for (j = 0; j < nsamples; j++)
	tempsample.push_back(functions[i]->RandomVector());

    // Sample points on a grid inside the bounding box of the function to
    // estimate the ISE \simeq \Delta\Sum_i(f(x_i)-\hat{f}(x_i))^2
    Vector<size_t> number(functions[i]->Dimension(), pow(nise, 1./functions[i]->Dimension()));
    grid = Grid(number.std(), functions[i]->LowerBounds(), functions[i]->UpperBounds());
    for (Vertex& vertex : grid.GetVertexArray())
	vertex.SetValue(functions[i]->Evaluate(vertex.GetCoordinates()));
    cellvol = grid.GetCellVolume();

    // Loop over the test parameter values
    for (j = start; j < end+1e-9; j += step) {
      // Initiliaze the estimator for the current parameter
      id = (j-start+1e-9)/step;
      tempalg = algo+std::to_string((int)j);
      DensityEstimator de(tempsample, tempalg, false, 1.);

      // Compute the log likelihood of the sample
      lnL = 0.;
      for (k = 0; k < tempsample.size(); k++)
	if ( de(tempsample[k]) )
            lnL += log(de(tempsample[k]));

      // Define the complexity penalization, (inversely proportional k~n/J)
      pen = (algo == "knn") ? (double)nsamples/j : j;

      // Fill the AIC for this value of the parameter
      aic[id] = -2*lnL+2*pen;
      gaic[i]->SetPoint((j-start+1e-9)/step, j, aic[id]);

      // Fill the BIC for this value of the parameter
      bic[id] = -2*lnL+2*log(nsamples)*pen;
      gbic[i]->SetPoint((j-start+1e-9)/step, j, bic[id]);

      // Fill the ISE for this value of the parameter
      ise[id] = 0.;
      for (Vertex& vertex : grid.GetVertexArray())
	  ise[id] += pow(de(vertex.GetCoordinates())-vertex.GetValue(), 2);
      ise[id] *= cellvol;
      gise[i]->SetPoint((j-start+1e-9)/step, j, ise[id]);

      std::cerr << j << "  " << tempalg << std::endl;
      std::cerr << "AIC: " << aic[id] << std::endl;
      std::cerr << "BIC: " << bic[id] << std::endl;
      std::cerr << "ISE: " << ise[id] << std::endl;
    }

    // Set markers where the minima are
    maic[i] = new TGraph();
    mbic[i] = new TGraph();
    mise[i] = new TGraph();
    for (TGraph* mingraph : {maic[i], mbic[i], mise[i]}) {
      mingraph->SetMarkerStyle(gaic[i]->GetMarkerStyle());
      mingraph->SetMarkerColor(gaic[i]->GetMarkerColor());
      mingraph->SetMarkerSize(2*gaic[i]->GetMarkerSize());
      mingraph->SetLineColor(0);
      mingraph->SetFillColor(0);
    }

    minid = std::distance(aic.begin(), std::min_element(aic.begin(), aic.end()));
    maic[i]->SetTitle(TString::Format("#hat{%s} = %d", var.c_str(), (int)(start+minid*step)));
    maic[i]->SetPoint(0, start+minid*step, aic[minid]);
    mgaic->Add(maic[i], "P");

    minid = std::distance(bic.begin(), std::min_element(bic.begin(), bic.end()));
    mbic[i]->SetTitle(TString::Format("#hat{%s} = %d", var.c_str(), (int)(start+minid*step)));
    mbic[i]->SetPoint(0, start+minid*step, bic[minid]);
    mgbic->Add(mbic[i], "P");

    minid = std::distance(ise.begin(), std::min_element(ise.begin(), ise.end()));
    mise[i]->SetTitle(TString::Format("#hat{%s} = %d", var.c_str(), (int)(start+minid*step)));
    mise[i]->SetPoint(0, start+minid*step, ise[minid]);
    mgise->Add(mise[i], "P");
  }

  // Draw and save the optimization plots
  for (TMultiGraph* mg : {mgaic, mgbic, mgise}) {
    TCanvas *canv = new TCanvas("c", "c", 1200, 800);
    gPad->SetGridx();
    gPad->SetGridy();
    mg->Draw("A");
    canv->BuildLegend(.75, .75, .89, .89);
    canv->SaveAs(mg->GetName()+TString(".pdf"));
    delete canv;
  }

  // Draw the information criterion alongside the MISE (double axes)
  TGraph *gic = gbic[0];
  TCanvas *canvs = new TCanvas("c", "c", 1200, 800);
  gise[0]->Draw("APL");
  gise[0]->SetTitle(";k;ISE");
  gise[0]->GetYaxis()->SetTitleOffset(0.8);
  gPad->SetLeftMargin(1.25);
  gise[0]->SetMarkerStyle(24);
  double xmax = gise[0]->GetXaxis()->GetXmax();
  double ymin = gise[0]->GetYaxis()->GetXmin();
  double ymax = gise[0]->GetYaxis()->GetXmax();
  double umin = gic->GetYaxis()->GetXmin();
  double umax = gic->GetYaxis()->GetXmax();
  for (size_t i = 0; i < gic->GetN(); i++)
      gic->GetY()[i] = ymin+(ymax-ymin)*(gic->GetY()[i]-umin)/(umax-umin); 
  gic->Draw("SAME PL");
  gic->SetLineColor(kRed+2);
  gic->SetMarkerStyle(25);

  TGaxis *axis = new TGaxis(xmax, ymin, xmax, ymax, umin, umax, 510, "+L");
  axis->SetTitle("BIC");
  axis->SetTitleOffset(1.3);
  axis->SetLineColor(kRed+2);
  axis->SetTextColor(kRed+2);
  axis->SetTextFont(42);
  axis->SetTextSize(0.04);
  axis->SetLabelFont(42);
  axis->SetLabelSize(0.04);
  axis->SetLabelColor(kRed+2);
  axis->Draw();


  canvs->SaveAs("knn_aic_ise.pdf");
  delete canvs;
			
  ////////////////////////////////////////////
  //////// OPTIMIZE FOR VARIABLE N ///////////
  ////////////////////////////////////////////
  // Initialize multigraphs that will host all the optimization graph
  TMultiGraph *mgopt = new TMultiGraph();
  mgopt->SetTitle(TString::Format(";n;%s^{*}", var.c_str()));

  // Power law to fit the trends with
  TF1 *fpow = new TF1("fpow", "[0]*pow(x, [1])", 0, 1e6);
  fpow->SetParameters(0., .5);
  TGraph *fleg;

  // Optimize kNN or the PBATDE for different values of N
  std::string criterion = "ise";
  std::vector<TGraphErrors*> gopt(functions.size());
  double a, b, c, d, fc, fd, exp;
  double gr = (1+sqrt(5))/2;
  std::vector<double> opt(npseudo);
  double mean, rms;
  TLegend* legopt = new TLegend(.12,.88-0.06*functions.size(), .4, .88);
  legopt->SetLineColorAlpha(0, 0);
  for (i = 0; i < functions.size(); i++) {

    // Initialize a TGraph to fill with the points for this function
    gopt[i] = new TGraphErrors();
    gopt[i]->SetTitle(TString::Format("%s, %s", functions[i]->Name().c_str(), algo.c_str()));
    gopt[i]->SetMarkerStyle(markers[i]);
    gopt[i]->SetMarkerSize(2);
    gopt[i]->SetLineColor(colors[i]);
    gopt[i]->SetLineWidth(2);
    gopt[i]->SetFillColor(0);
    mgopt->Add(gopt[i], "P");

    // Sample points on a grid inside the bounding box of the function to
    // estimate the ISE \simeq \Delta\Sum_i(f(x_i)-\hat{f}(x_i))^2
    Vector<size_t> number(functions[i]->Dimension(), pow(nise, 1./functions[i]->Dimension()));
    grid = Grid(number.std(), functions[i]->LowerBounds(), functions[i]->UpperBounds());
    for (Vertex& vertex : grid.GetVertexArray())
	vertex.SetValue(functions[i]->Evaluate(vertex.GetCoordinates()));
    cellvol = grid.GetCellVolume();

    // Produce samples of j points with j ranging from 500 to 128k
    id = 0;
    for (j = 1e3; j < 2e5; j *= 2) {
      for (k = 0; k < npseudo; k++) {

        // Sample the requested amount of samples from the current distribution
        tempsample.resize(0);
        for (l = 0; l < j; l++)
	    tempsample.push_back(functions[i]->RandomVector());

        // Find the optimal parameter in a certain range by using the golden-section search
        a = start;
        b = std::min(end, j);
        c = b - (b-a)/gr;
        d = a + (b-a)/gr;
        while ( fabs(c - d) >= 1 ) {
	  // First need to compute the ISE for par = c and par = d
          // Initiliaze the estimator for the current parameter
          tempalg = algo+std::to_string((int)c);
          DensityEstimator dec(tempsample, tempalg, false, 1.);
          fc = 0;
          for (Vertex& vertex : grid.GetVertexArray())
	      fc += pow(dec(vertex.GetCoordinates())-vertex.GetValue(), 2);
	  fc *= cellvol;

          tempalg = algo+std::to_string((int)d);
          DensityEstimator ded(tempsample, tempalg, false, 1.);
          fd = 0;
          for (Vertex& vertex : grid.GetVertexArray())
	      fd += pow(ded(vertex.GetCoordinates())-vertex.GetValue(), 2);
	  fd *= cellvol;
        
  	  // Update the boundaries of the search
          if ( fc < fd ) {
            b = d;
          } else {
            a = c;
	  }  

          // c and d are recomputed here to avoid loss of precision
          c = b - (b-a)/gr;
          d = a + (b-a)/gr;
        }

	opt[k] = (b+a)/2;
      }

      // Fill with the optimal value of the parameter
      mean = Math::Mean(opt);
      rms = Math::RMS(opt)/sqrt(npseudo);
      std::cerr << "Optimal par: " << j << "  " << mean << " +/- " << rms << std::endl;
      gopt[i]->SetPoint(id, j, mean);
      gopt[i]->SetPointError(id, 0, rms);
      id++;
    }

    // Fit the graph with a power law, add a legend entry with its result
    fpow->SetParameter(0, 1.);
    fpow->FixParameter(1., 2/(2.+functions[i]->Dimension()));
    gopt[i]->Fit(fpow);
    gopt[i]->GetFunction("fpow")->SetLineColor(gopt[i]->GetLineColor());
    a = gopt[i]->GetFunction("fpow")->GetParameter(0);
    exp = gopt[i]->GetFunction("fpow")->GetParameter(1);

    legopt->AddEntry(gopt[i], functions[i]->Title().c_str(), "LP");

/*    fleg = new TGraph(*gopt[i]);
    fleg->SetLineColorAlpha(0, 0);
    fleg->SetMarkerColorAlpha(0, 0);
    fleg->SetFillColor(0);
    fleg->SetTitle(TString::Format("%0.2fn^{%0.2f}", a, exp));
    mgopt->Add(fleg);*/
  }

  // Draw and save the optimised parameter
  TCanvas* canv = new TCanvas("c", "c", 1200, 800);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogx();
  gPad->SetLogy();

  mgopt->Draw("A");
  legopt->Draw("SAME");

  canv->SaveAs(TString::Format("%s_opt_%dD_tri.pdf", algo.c_str(), (int)functions[0]->Dimension()));

  TFile *outfile = new TFile(TString::Format(
	"%s_opt_%dD_tri.root", algo.c_str(), (int)functions[0]->Dimension()), "RECREATE");
  mgopt->Write();
  outfile->Close();
  delete canv;

  return 0;
}
