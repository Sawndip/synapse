// Cpp includes
#include <ctime>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <regex>
#include <map>

// Root includes
#include "TStyle.h"
#include "TObject.h"
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TColor.h"

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


// TEMPORARYYYYYYYY TODO
double gauslevel(double *x, double *par) {

  double R2 = TMath::ChisquareQuantile(x[0], par[0]);
  return exp(-R2/2)/pow(2*M_PI, par[0]/2);
}

int main(int argc, char ** argv) {

  gStyle->SetLabelSize(0.05, "XYZ");
  gStyle->SetTitleSize(0.05, "XYZ");

  std::vector<double> stops = {0., 0.125, .25, .375, .5, .625, .75, .875, 1.};
  std::vector<double> red = {0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764};
  std::vector<double> green = {0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832};
  std::vector<double> blue = {0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539};
  TColor::CreateGradientColorTable(9, &stops[0], &red[0], &green[0], &blue[0], 255, 1);

  ////////////////////////////////////////////
  //////////// GLOBAL PARAMETERS /////////////
  ////////////////////////////////////////////
  gStyle->SetOptStat(0);
  double alpha = .5;
  const size_t nsamples = 1e4;
  size_t npoints = 1e4;
  std::string algo = "knn";
  if ( argc < 2 ) {
    std::cerr << "No density estimation algorithm specified, using knn" << std::endl;
  } else {
    algo = argv[1];
  }

  ////////////////////////////////////////////
  /////////// TEST DISTRIBUTIONS /////////////
  ////////////////////////////////////////////
  // 1D distributions
  std::vector<DFunction*> f1d;	// List of 1D functions

  // Gaussian distribution
//  double gaus1m(0.), gaus1s(1.);
//  f1d.push_back(new DGaus(gaus1m, gaus1s));

  // Two-peak gaussian
//  f1d.push_back(new DMultiGaus(1, 2));

  // Chi squared distribution
//  f1d.push_back(new DChiSquared(4));

  // Cauchy
//  f1d.push_back(new DCauchy());

  // Exponential
//  f1d.push_back(new DExponential());

  // Uniform
  f1d.push_back(new DUniform());
  f1d.back()->SetRange(-1, 1);

  // Triangular
//  f1d.push_back(new DTriangular());
//  DTriangular tri1(1);
//  std::cerr << "Triangular radius: " << tri1.Radius(.5) << std::endl;
//  f1d.back()->SetRange(-tri1.Radius(.5), tri1.Radius(.5));

  // Maxwell-Boltzmann
//  f1d.push_back(new DMaxwell());

  ////////////////////////////////////////////
  // 2D distributions
  std::vector<DFunction*> f2d;	// List of 2D functions

  // Uncorrelated Gaussian distribution
//  Matrix<double> gaus2_cov({{1, 0}, {0, 1}});
//  std::vector<double> gaus2_means = {0, 0};
//  f2d.push_back(new DGaus(gaus2_means, gaus2_cov));

  // Correlated Gaussian distribution
/*  double gaus2_rho = 0.9;
  gaus2_cov = Matrix<double>({{1, gaus2_rho}, {gaus2_rho, 1}});
  gaus2_means = {0, 0};
  DGaus *fgaus2corr = new DGaus(gaus2_means, gaus2_cov);
  f2d.push_back(fgaus2corr);
  fgaus2corr->SetName("gaus2corr");
  fgaus2corr->SetTitle("2D Correlated Gaussian");*/

  // Three-peak Gaussian
//  f2d.push_back(new DMultiGaus(2, 3));

  // Exponential
  f2d.push_back(new DExponential(2));

  // Uniform
//  f2d.push_back(new DUniform(2));

  // Triangular
//  f2d.push_back(new DTriangular(2));

  // Maxwell-Boltzmann
  f2d.push_back(new DMaxwell({0., 0.}, Matrix<double>({{1, 0}, {0, 1}})));

  ////////////////////////////////////////////
  // 3D distributions
  std::vector<DFunction*> f3d;	// List of 3D functions

  // Gaussian distribution
/*  double Lgaus3 = 3;
  Matrix<double> gaus3_cov({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}});
  std::vector<double> gaus3_means({0, 0, 0}),
	gaus3_l({-Lgaus3, -Lgaus3, -Lgaus3}), gaus3_u({Lgaus3, Lgaus3, Lgaus3});
  DGaus *fgaus3 = new DGaus(gaus3_means, gaus3_cov);
  fgaus3->SetRange(gaus3_l, gaus3_u);
  f3d.push_back(fgaus3);*/

  // Four-peak Gaussian
//  f3d.push_back(new DMultiGaus(3, 4));

  // Uniform
//  f3d.push_back(new DUniform(3));

  // Exponential
//  f3d.push_back(new DExponential(3));

  // Triangular
//  f3d.push_back(new DTriangular(3));

  // Maxwell-Boltzmann
//  f3d.push_back(new DMaxwell(3));

  ////////////////////////////////////////////
  // Radial distributions (for higher dimensions, only way to represent functions is radially)
  std::vector<DFunction*> frad;	// List of Radial functions

  Pitch::setAnOutput(Pitch::debug, std::cerr);
//  frad.push_back(new DGaus(1));
//  frad.push_back(new DGaus(2));
//  frad.push_back(new DGaus(3));
  frad.push_back(new DGaus(4));
//  frad.push_back(new DGaus(6));

//  frad.push_back(new DTriangular(1));
//  frad.push_back(new DTriangular(2));
//  frad.push_back(new DTriangular(3));
//  frad.push_back(new DTriangular(4));

//  frad.push_back(new DUniform(1));
//  frad.push_back(new DUniform(4));
//  frad.back()->SetRange({-1,-1,-1,-1},{1,1,1,1});
//  frad.push_back(new DMaxwell(1));
//  frad.push_back(new DMaxwell(2));
//  frad.push_back(new DMaxwell(3));
//  frad.push_back(new DMaxwell(4));

  // Gaussian distribution
/*  double Lgaus4 = 3;
  Matrix<double> gaus4_cov({{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}});
  std::vector<double> gaus4_means({0, 0, 0, 0}),
	gaus4_l({-Lgaus4, -Lgaus4, -Lgaus4, -Lgaus4}), gaus4_u({Lgaus4, Lgaus4, Lgaus4, Lgaus4});
  DGaus *fgaus4 = new DGaus(gaus4_means, gaus4_cov);
//  fgaus4->SetRange(gaus4_l, gaus4_u);
  frad.push_back(fgaus4);*/

  ////////////////////////////////////////////
  //////////////// SAMPLING //////////////////
  ////////////////////////////////////////////
  // 1D samples
  std::vector<std::vector<std::vector<double>>> s1d(f1d.size());
  size_t i, j;
  for (i = 0; i < f1d.size(); i++) 
    for (j = 0; j < nsamples; j++)
        s1d[i].push_back(f1d[i]->RandomVector());

  ////////////////////////////////////////////
  // 2D samples
  std::vector<std::vector<std::vector<double>>> s2d(f2d.size());
  for (i = 0; i < f2d.size(); i++) 
    for (j = 0; j < nsamples; j++)
        s2d[i].push_back(f2d[i]->RandomVector());

  ////////////////////////////////////////////
  // 3D samples
  std::vector<std::vector<std::vector<double>>> s3d(f3d.size());
  for (i = 0; i < f3d.size(); i++) 
    for (j = 0; j < nsamples; j++)
        s3d[i].push_back(f3d[i]->RandomVector());

  ////////////////////////////////////////////
  // Radial samples
  std::vector<std::vector<std::vector<double>>> srad(frad.size());
  for (i = 0; i < frad.size(); i++) 
    for (j = 0; j < nsamples; j++)
        srad[i].push_back(frad[i]->RandomVector());

  DensityEstimator de(srad[0], algo, false, 1.);
  std::vector<double> levels(nsamples);
  for (i = 0; i < nsamples; i++)
      levels[i] = de(srad[0][i]);

  std::sort(levels.begin(), levels.end());
  TGraph* gpro = new TGraph();
  for (i = 0; i < nsamples; i++)
      gpro->SetPoint(i, (double)(nsamples-i)/(nsamples+1), levels[i]);

  TF1 *flevel = new TF1("level", gauslevel, 0, 1, 1);
  flevel->SetParameter(0, 4);

  TCanvas *ct = new TCanvas("c", "c", 1200, 800);
  gpro->Draw("APL");
  flevel->Draw("SAME");
  ct->SaveAs("profile.pdf");
  delete ct;

  // Plot the 3D exponential samples in xy, xz, yz
  /*TH2F* hxy = new TH2F("xy", "", 100, -3, 3, 100, -3, 3);
  TH2F* hxz = new TH2F("xz", "", 100, -3, 3, 100, -3, 3);
  TH2F* hyz = new TH2F("yz", "", 100, -3, 3, 100, -3, 3);
  TH3F* hxyz = new TH3F("xyz", "", 10, -3, 3, 10, -3, 3, 10, -3, 3);
  std::vector<std::vector<double>> test = srad[3];
  for (i = 0; i < nsamples; i++) {
    hxy->Fill(test[i][0], test[i][1]);
    hxz->Fill(test[i][0], test[i][2]);
    hyz->Fill(test[i][1], test[i][2]);
    hxyz->Fill(test[i][0], test[i][1], test[i][2]);
  }

  TCanvas *ctest = new TCanvas("c", "c", 1200, 800);
  hxy->Draw("COLZ");
  ctest->SaveAs("xy.pdf");
  delete ctest;

  ctest = new TCanvas("c", "c", 1200, 800);
  hxz->Draw("COLZ");
  ctest->SaveAs("xz.pdf");
  delete ctest;

  ctest = new TCanvas("c", "c", 1200, 800);
  hyz->Draw("COLZ");
  ctest->SaveAs("yz.pdf");
  delete ctest;

  ctest = new TCanvas("c", "c", 1200, 800);
  hxyz->Draw("");
  ctest->SaveAs("xyz.pdf");
  delete ctest;*/

  ////////////////////////////////////////////
  /////////// DENSITY ESTIMATION /////////////
  ////////////////////////////////////////////
  // 1D density estimators
  std::vector<TGraph*> de1d;	// Interpolated function
  std::vector<TGraph*> band1d;	// Error band with respect to true
  std::vector<TH1F*> cont1d;	// 1D contour
  std::vector<double> ise1d;	// Mean Integrated Square Error
  std::vector<double> vol1d;	// Volumes of the 1d contours

  std::vector<double> min, max;
  double x, ise, real, est;
  for (i = 0; i < f1d.size(); i++) {
    // Initialize the density estimator
    DensityEstimator de(s1d[i], algo, false);

    // Initialize graphs and fill it with points from the density estimator
    TGraph *graph = new TGraph(npoints);
    graph->SetTitle(f1d[i]->Title().c_str());
    TGraph *band = new TGraph(2*npoints);
    f1d[i]->Range(min, max);
    ise = 0.;
    for (j = 0; j < npoints; j++) {
      x = min[0]+(double)j*(max[0]-min[0])/(npoints-1);
      real = f1d[i]->Evaluate(x);
      est = de(x);
      ise += pow(est-real, 2)*(max[0]-min[0])/(npoints-1);

      graph->SetPoint(j, x, est);
      band->SetPoint(j, x, est);
      band->SetPoint(2*npoints-1-j, x, real);
    }

    de1d.push_back(graph);
    band1d.push_back(band);
    ise1d.push_back(ise);
//    vol1d.push_back(de.ContourVolume(alpha, "mc"));
//    cont1d.push_back(de.Contour1D());
  }

  ////////////////////////////////////////////
  // 2D density estimators
  npoints = 100;
  std::vector<TH2F*> de2d;	// Interpolated function
  std::vector<TH2F*> cont2d;	// TH2F at a certain contour level
  std::vector<double> ise2d;	// Mean Integrated Square Error
  std::vector<double> vol2d;	// Volumes of the 2d contours

  double y;
  size_t k;
  for (i = 0; i < f2d.size(); i++) {
    // Initialize the density estimator
    DensityEstimator de(s2d[i], algo, false, 1.);

    // Initialize graphs and fill it with points from the density estimator
    f2d[i]->Range(min, max);
    TH2F *graph = new TH2F((f2d[i]->Name()+"_de").c_str(), f2d[i]->Title().c_str(),
			   npoints, min[0], max[0], npoints, min[1], max[1]);
    ise = 0.;
    for (j = 0; j < (size_t)graph->GetNbinsX(); j++) {
      for (k = 0; k < (size_t)graph->GetNbinsY(); k++) {
      	x = graph->GetXaxis()->GetBinCenter(j+1);
      	y = graph->GetYaxis()->GetBinCenter(k+1);
      	real = f2d[i]->Evaluate({x, y});
      	est = de({x, y});
      	ise += pow(est-real, 2)*(max[0]-min[0])*(max[1]-min[1])/pow(npoints-1, 2);

      	graph->SetBinContent(j+1, k+1, est);
      }
    }

    de2d.push_back(graph);
    ise2d.push_back(ise);
    vol2d.push_back(de.ContourVolume(alpha, "mc"));
    cont2d.push_back(de.Contour2D());
  }

  ////////////////////////////////////////////
  // 3D density estimators
  std::vector<TH3F*> cont3d;	// TH3F at a certain contour level
  std::vector<double> ise3d;	// Mean Integrated Square Error
  std::vector<double> vol3d;	// Volumes of the 3d contours

  for (i = 0; i < f3d.size(); i++) {
    // Initialize the density estimator
    DensityEstimator de(s3d[i], algo, false, 1.);
    f3d[i]->Range(min, max);

    ise3d.push_back(0.); // TODO TODO TODO ???
    vol3d.push_back(de.ContourVolume(alpha, "mc"));
    cont3d.push_back(de.Contour3D());
  }

  ////////////////////////////////////////////
  // Radial density estimators
  npoints = 1000;
  std::vector<TGraph*> derad;	// Interpolated radial function
  std::vector<TGraph*> bandrad;	// Error band with respect to true
  std::vector<TH1F*> contrad;	// Radial contour
  std::vector<double> iserad;	// Mean Integrated Square Error
  std::vector<double> volrad;	// Volumes of the contours

  double r;
  std::vector<double> v;
  for (i = 0; i < frad.size(); i++) {
    // Initialize the density estimator
    DensityEstimator de(srad[i], algo, false, 1.);
    v.resize(de.GetDimension());

    // Initialize graphs and fill it with points from the density estimator
    TGraph *graph = new TGraph(npoints);
    graph->SetTitle(frad[i]->Title().c_str());
    TGraph *band = new TGraph(2*npoints);
    frad[i]->Range(min, max);
    ise = 0.;
    for (j = 0; j < npoints; j++) {
      r = (double)j*max[0]/(npoints-1);
      real = frad[i]->Radial(r);
      v[0] = r;
      if ( de.GetDimension() > 1 ) {
	v[0] = r;
	v[1] = 0;
      }
      est = de(v);
      ise += pow(r, de.GetDimension()-1)*pow(est-real, 2)*max[0]/(npoints-1);
      graph->SetPoint(j, r, est);
      band->SetPoint(j, r, est);
      band->SetPoint(2*npoints-1-j, r, real);
    }

    derad.push_back(graph);
    bandrad.push_back(band);
//    iserad.push_back(de.GetDimension()*pow(2, de.GetDimension())*ise);
    iserad.push_back(de.GetDimension()*Math::UnitBallVolume(de.GetDimension(), 2)*ise);
//    volrad.push_back(de.ContourVolume(alpha, "mc"));
//    contrad.push_back(de.Contour1D());
  }

  ////////////////////////////////////////////
  /////////////// COMPARAISON ////////////////
  ////////////////////////////////////////////
  // 1D distributions
  for (i = 0; i < f1d.size(); i++) {

    // Initialize the Canvas
    TCanvas *c = new TCanvas("c", "c", 1200, 800);
//    gPad->SetLogy();
    gPad->SetLeftMargin(.125);

    // Set the style
    de1d[i]->SetLineColor(4);

    band1d[i]->SetLineColor(0);
    band1d[i]->SetFillColorAlpha(4, 0.25);

    // Legend for the value of ISE and CRE
    TLegend *leg = new TLegend(.65, .8, .89, .89);
    leg->SetLineColorAlpha(0, 0);
    leg->SetFillStyle(3001);
    leg->SetFillColorAlpha(0, 0);
    leg->AddEntry(band1d[i], TString::Format("ISE: %0.6f", ise1d[i]), "f");
    double truevol = f1d[i]->ContourVolume(alpha);
//    leg->AddEntry(cont1d[i], TString::Format("CRE(%d%%): %0.4f%%",
//	(int)(100*alpha), 100*(vol1d[i]-truevol)/truevol), "l");

    // Draw
    de1d[i]->Draw("AL");
    de1d[i]->SetTitle(";x;#rho(x)");
    f1d[i]->Draw("L SAME");
    band1d[i]->Draw("F SAME");
//    for (TLine* line : f1d[i]->Contour1D(alpha))
//	line->Draw("SAME");
//    cont1d[i]->Draw("E SAME");
    leg->Draw("SAME");
    c->SaveAs(f1d[i]->Name()+TString(".pdf"));
    delete c;
    delete leg;
  }


  ////////////////////////////////////////////
  // 2D distributions
  for (i = 0; i < f2d.size(); i++) {

    // Save the true and estimated distributions separately
    TCanvas *cdist = new TCanvas("c", "c", 1200, 800);
    gPad->SetRightMargin(.15);
    de2d[i]->GetXaxis()->SetRangeUser(-3, 3);
    de2d[i]->GetYaxis()->SetRangeUser(-3, 3);
    de2d[i]->Draw("CONTZ");
    de2d[i]->SetTitle(";x;y;#rho(x,y)");
    cdist->SaveAs(f2d[i]->Name()+TString("_de.pdf"));
    delete cdist;

    cdist = new TCanvas("c", "c", 1200, 800);
    f2d[i]->Draw("CONTZ");
    cdist->SaveAs(f2d[i]->Name()+TString("_true.pdf"));
    delete cdist;

    // Initialize the Canvas
    TCanvas *c = new TCanvas("c", "c", 1200, 800);

    // Set the style
    de2d[i]->SetLineWidth(2);
    de2d[i]->SetLineColor(4);

    // Legend for the value of ISE and CRE
    TLegend *leg = new TLegend(.6, .72, .89, .89);
    leg->SetLineColorAlpha(0, 0);
    leg->SetFillStyle(3001);
    leg->SetFillColorAlpha(0, 0);
    leg->AddEntry(de2d[i], TString::Format("ISE: %0.6f", ise2d[i]), "");
    double truevol = f2d[i]->ContourVolume(alpha);
    leg->AddEntry(cont2d[i], TString::Format("CRE(%d%%): %0.2f%%",
	(int)(100*alpha), 100*(vol2d[i]-truevol)/truevol), "l");

    // Draw
    de2d[i]->Draw("CONT2");
    de2d[i]->SetLineColor(0);
//    f2d[i]->Draw("CONT2 SAME");
    leg->Draw("SAME");
    for (TObject* obj : f2d[i]->Contour2D(alpha))
	obj->Draw("SAME");
    cont2d[i]->Draw("CONT3 SAME");
    c->SaveAs(f2d[i]->Name()+TString(".pdf"));
    delete c;
    delete leg;
  }

  ////////////////////////////////////////////
  // 3D distributions
  for (i = 0; i < f3d.size(); i++) {
    // Initialize the Canvas
    TCanvas *c = new TCanvas("c", "c", 1200, 800);

    // Legend for the value of ISE and CRE
    TLegend *leg = new TLegend(.6, .75, .89, .89);
    leg->SetLineColorAlpha(0, 0);
    leg->SetFillStyle(3001);
    leg->SetFillColorAlpha(0, 0);
//    leg->AddEntry(de2d[i], TString::Format("ISE: %0.6f", ise2d[i]), "l");
    double truevol = f3d[i]->ContourVolume(alpha);
    leg->AddEntry(cont3d[i], TString::Format("CRE(%d%%): %0.4f%%",
	(int)(100*alpha), 100*(vol3d[i]-truevol)/truevol), "lf");

    // Draw
    f3d[i]->Contour3D(alpha)->Draw("FB");
    cont3d[i]->Draw("SAME ISO FBBB");
    leg->Draw("SAME");
    c->SaveAs(f3d[i]->Name()+TString(".pdf"));
    delete c;
    delete leg;
  }

  ////////////////////////////////////////////
  // Radial distributions
  for (i = 0; i < frad.size(); i++) {

    // Initialize the Canvas
    TCanvas *c = new TCanvas("c", "c", 1200, 800);
//    gPad->SetLogy();

    // Set the style
    derad[i]->SetLineColor(4);

    bandrad[i]->SetLineColor(0);
    bandrad[i]->SetFillColorAlpha(4, 0.25);

    // Legend for the value of ISE and CRE
    TLegend *leg = new TLegend(.6, .82, .89, .89);
    leg->SetLineColorAlpha(0, 0);
    leg->SetFillStyle(3001);
    leg->SetFillColorAlpha(0, 0);
    leg->AddEntry(bandrad[i], TString::Format("ISE: %0.6f", iserad[i]), "f");
//    double truevol = frad[i]->ContourVolume(alpha);
//    leg->AddEntry(contrad[i], TString::Format("CRE(%d%%): %0.4f%%",
//	(int)(100*alpha), 100*(volrad[i]-truevol)/truevol), "l");

    // Draw
    derad[i]->Draw("AL");
    derad[i]->SetTitle(";r;#hat{#rho}(r)");
    gPad->SetLeftMargin(0.12);
    frad[i]->DrawRadial("L SAME");
    bandrad[i]->Draw("F SAME");
//    for (TLine* line : frad[i]->ContourRadial(alpha))
//	line->Draw("SAME");
//    contrad[i]->Draw("E SAME");
    leg->Draw("SAME");
    c->SaveAs(frad[i]->Name()+TString("r.pdf"));
    delete c;
    delete leg;
  }
}
