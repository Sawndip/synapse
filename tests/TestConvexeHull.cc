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

// Additional modules
#include "Bunch.hh"
#include "Statistics.hh"
#include "Geometry.hh"
#include "DensityEstimator.hh"
#include "AlphaComplex.hh"
#include "DGaus.hh"
#include "DMultiGaus.hh"
#include "DChiSquared.hh"
#include "DCauchy.hh"
#include "DExponential.hh"
#include "DUniform.hh"
#include "DTriangular.hh"
#include "DMaxwell.hh"

// Lebesgue p-norm to the power of p
double LpMagp(const double& p, const std::vector<double>& x) {

  double mag(0.);
  for (size_t i = 0; i < x.size(); i++)
      mag += pow(fabs(x[i]), p);
  return mag;
}

double ChiSquaredDist(double x, double d) {

  return  exp(-x/2)*pow(x, d/2-1)/pow(2, d/2)/tgamma(d/2);
}

int main(int argc, char ** argv) {

  ////////////////////////////////////////////
  //////////// GLOBAL PARAMETERS /////////////
  ////////////////////////////////////////////
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(.05, "XY");
  gStyle->SetLabelSize(.05, "XY");
  size_t n = 1e5;		// Number of points in the hull
  size_t d = 1;			// Dimensionality of the space
  size_t npseudo = 4e2;		// Number of pseudo experiments
  TRandom3 rdmzer(time(NULL));

  Beam::Bunch beam;	// Requested for it to compile, undefined reference otherwise (???) TODO TODO TODO

  ////////////////////////////////////////////
  //////////// GLOBAL PARAMETERS /////////////
  ////////////////////////////////////////////
  // Marker and color styles to be used
//  std::vector<int> colors = {kBlack, kMagenta+1, kBlue+1, kCyan+1, 
//				kGreen+1, kOrange+1, kRed+1, kViolet+1, kPink+1};
  std::vector<int> colors = {kBlack, kRed+2, kGreen+2, kBlue+2, kMagenta+2, kCyan+2, kYellow+2};
  std::vector<int> markers = {24, 25, 26, 32, 20, 21, 22, 23};
  size_t styleid(0);

  ////////////////////////////////////////////
  /////////// ARBITRARY SHAPES ///////////////
  ////////////////////////////////////////////
  // General multi graph to plot the different results on the same canvas
  TMultiGraph *mgmean = new TMultiGraph();
  TMultiGraph *mgrms = new TMultiGraph();

  // Set up the legends
  TLegend *legmean = new TLegend(.15, .15, .4, .4);
  legmean->SetLineColor(0);
  TLegend *legrms = new TLegend(.15, .15, .4, .4);
  legrms->SetLineColor(0);

  std::cerr << TMath::ChisquareQuantile(.1, 3) << std::endl;

  // Sample points in a p-unit n-ball, compute the volume of their convexe hull, measure the spread
  std::string dist = "tgaus"; // uni:Uniform, tgaus:Truncated Gaussian
  double alpha = .09;
  if ( dist == "uni" )
      alpha = 1.;
  size_t mind(1), maxd(4);
  std::vector<double> points;
  std::vector<double> vols(npseudo);
  std::vector<double> pvalues = {2.};
  std::vector<std::map<double, TGraphErrors*>> gball(maxd-mind+1);
  std::vector<std::map<double, TGraphErrors*>> gball_rms(maxd-mind+1);
  std::vector<double> point;
  double truevol(1.);
  size_t i, j;
  size_t id;
  for (d = mind; d < maxd+1; d++) {

    point.resize(d);
    for (const double& p : pvalues) {
      // Compute the true volume of the p-unit n-ball
      if ( dist == "uni" ) {
        truevol = pow(2*std::tgamma(1.+1./p), d)/std::tgamma(1.+(double)d/p);
      } else if ( dist == "tgaus" ) {
	double R2 = TMath::ChisquareQuantile(alpha, d);
	truevol = pow(M_PI, (double)d/2)*pow(R2, (double)d/2)/std::tgamma(1.+(double)d/2);
      }

      // Set the style
      gball[d-mind][p] = new TGraphErrors();
      gball_rms[d-mind][p] = new TGraphErrors();
      for (TGraphErrors* graph : {gball[d-mind][p], gball_rms[d-mind][p]}) {
        graph->SetMarkerStyle(markers[styleid]);
        graph->SetMarkerSize(1.8);
        graph->SetLineColor(colors[styleid]);
        graph->SetLineWidth(2);
        graph->SetFillColor(0);
	if ( dist == "uni" ) {
          graph->SetTitle(TString::Format("%d-unit %d-ball", (int)p, (int)d));
	} else if ( dist == "tgaus" ) {
//          graph->SetTitle(TString::Format("%d-Gaussian, %0.0f %%", (int)d, 100*alpha));
          graph->SetTitle(TString::Format("%d-Gaussian", (int)d));
	}
      }
      ++styleid;
      mgmean->Add(gball[d-mind][p], "EP");
      legmean->AddEntry(gball[d-mind][p], gball[d-mind][p]->GetTitle(), "lp");
      mgrms->Add(gball_rms[d-mind][p], "EP");
      legrms->AddEntry(gball_rms[d-mind][p], gball_rms[d-mind][p]->GetTitle(), "lp");

      id = 0;
      for (n = 100; n < 1.5e5; n *= 2) {
        for (i = 0; i < npseudo; i++) {
	  points.resize(0);
	  if ( dist == "uni" ) {
 	    // Generate uniform points in a ball of Lebesgue metric p
            while ( points.size() < n*d ) {
              for (j = 0; j < d; j++)
	          point[j] = rdmzer.Uniform(-1, 1);
	      if ( LpMagp(p, point) < 1 )
		  points.insert(points.end(), point.begin(), point.end());
            }
	  }

	  if ( dist == "tgaus" ) {
	    // Generate Gaussian points, keep a fraction alpha of them only
	    std::vector<std::pair<Vector<double>, double>> extpoints;
	    std::vector<double> mean(d, 0.);
	    Matrix<double> covmat(d, d, 0.);
            while ( extpoints.size() < n ) {
              for (j = 0; j < d; j++)
	          point[j] = rdmzer.Gaus(0, 1);
//	          point[j] = rdmzer.Uniform(-1, 1);
//	      if ( LpMagp(p, point) > 1 )
//		  continue;
	      Math::IncrementCovarianceMatrix(extpoints.size(), covmat, mean, point);
              for (j = 0; j < d; j++)
		  Math::IncrementMean(extpoints.size(), mean[j], point[j]);
	      extpoints.push_back(std::pair<Vector<double>, double>(point, 0.));
            }

	    // Compute the amplitudes and retain the core alpha points
	    Vector<double> centred;
	    covmat.Invert();

	    for (j = 0; j < extpoints.size(); j++) {
	      centred = extpoints[j].first-Vector<double>(mean);
	      extpoints[j].second = centred*(covmat*centred);
	    }
 
	    std::sort(extpoints.begin(), extpoints.end(),
		[] (const std::pair<Vector<double>, double>& a,
		    const std::pair<Vector<double>, double>& b) {
		return a.second < b.second; });

	    for (j = 0; j < alpha*(extpoints.size()+1)-1; j++)
		points.insert(points.end(), extpoints[j].first.std().begin(),
					    extpoints[j].first.std().end());

	    if ( n == 1600 && i == 0 && d == 1 ) {
		TH1F* htotal = new TH1F("htotal", "", 100, -5, 5);
		TH1F* hselec = new TH1F("htotal", "", 100, -5, 5);
		for (j = 0; j < extpoints.size(); j++)
		    htotal->Fill(extpoints[j].first[0]);
		for (j = 0; j < points.size(); j++)
		    hselec->Fill(points[j]);
		
		std::cerr << htotal->GetEntries() << "  " << hselec->GetEntries() << std::endl;
		TCanvas *c = new TCanvas("c", "c", 1200, 800);
		htotal->Draw();
		hselec->Draw("SAME");
		c->SaveAs("selec.pdf");
		delete c;
	    }
	    if ( n == 1600 && i == 0 && d == 2 ) {
		TH2F* htotal = new TH2F("htotal", "", 100, -5, 5, 100, -5, 5);
		TH2F* hselec = new TH2F("hselec", "", 100, -5, 5, 100, -5, 5);
		hselec->SetMarkerColor(2);
		for (j = points.size()/2; j < extpoints.size(); j++)
		    htotal->Fill(extpoints[j].first[0], extpoints[j].first[1]);
		for (j = 0; j < points.size()/2; j++) {
		    hselec->Fill(points[2*j], points[2*j+1]);
		}
		
		TCanvas *c = new TCanvas("c", "c", 1200, 800);
		htotal->Draw("SCAT");
		hselec->Draw("SCAT SAME");
		TEllipse *ell = new TEllipse(0, 0, 1., 1.);
		ell->SetFillStyle(0);
		ell->Draw("SAME");
		c->SaveAs("selec2.pdf");
		delete c;
	    }
	  }

          // Compute the convex hull
          if ( d == 1 ) {
            std::sort(points.begin(), points.end());
            vols[i] = (points.back()-points.front());
            continue;
          }

          orgQhull::Qhull qhull("", d, points.size()/d, &points[0], "");
          vols[i] = qhull.volume();
        }

        gball[d-mind][p]->SetPoint(id, n*alpha, -(Math::Mean(vols)-truevol)/truevol);
        gball[d-mind][p]->SetPointError(id, 0, Math::MeanError(vols)/truevol);
	std::cerr << d << "  " << p << "  " << n << "  " << (Math::Mean(vols)-truevol)/truevol 
		  << "  " << Math::MeanError(vols)/truevol << std::endl;
        gball_rms[d-mind][p]->SetPoint(id, n, Math::RMS(vols)/Math::Mean(vols));
        gball_rms[d-mind][p]->SetPointError(id, 0, Math::RMSError(vols)/Math::Mean(vols));
        ++id;
      }
    }
  }

  // Expected trend
  TF1 *fn[maxd-mind+1];

  // Draw all the graphs
  TCanvas *c = new TCanvas("c", "c", 1200, 800);
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  mgmean->Draw("A");
  mgmean->SetTitle(";#alphan;-E(#hat{V}-V)/V");
//    fn[d-mind] = new TF1("fn", "pow(1+[0]*x, [1])", 10, 1e6);
  for (d = mind; d < maxd+1; d++) {
//    fn[d-mind] = new TF1("fn", "[0]*pow(x, [1])", 10, 1e6);
    fn[d-mind] = new TF1("fn", "pow(1+pow([0], 1./[1])*x, [1])", 1, 1e6);
    double C = Math::HullVolumeUniformFactor(d, 2);
    fn[d-mind]->FixParameter(0, C);
    fn[d-mind]->FixParameter(1, -2./(d+1));
    if ( dist == "uni" ) {
      fn[d-mind]->SetTitle(TString::Format("%0.2fn^{%0.2f}",
      	fn[d-mind]->GetParameter(0), fn[d-mind]->GetParameter(1)));
    } else if ( dist == "tgaus" ) {
      fn[d-mind]->SetTitle(TString::Format("%0.2f(#alphan)^{%0.2f}",
	fn[d-mind]->GetParameter(0)/alpha, fn[d-mind]->GetParameter(1)));
    }
    gball[d-mind][2]->Fit("fn");
    fn[d-mind]->SetLineColor(gball[d-mind][2]->GetLineColor());
    fn[d-mind]->Draw("SAME");
  }
  legmean->Draw("SAME");
  c->SaveAs("hull_mean.pdf");
  delete c;

  c = new TCanvas("c", "c", 1200, 800);
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  mgrms->Draw("A");
  mgrms->SetTitle(";n;-RMS(#hat{V})/#hat{V}");
  for (d = mind; d < maxd+1; d++) {
    fn[d-mind] = new TF1("fn", "[0]*pow(x, [1])", 1, 1e6);
    if ( dist == "uni" ) {
      fn[d-mind]->FixParameter(1, -(d+3.)/(2*(d+1.)));
      fn[d-mind]->FixParameter(0, sqrt(2));
    } else if ( dist == "tgaus" ) {
      double R2 = TMath::ChisquareQuantile(alpha, d);
      fn[d-mind]->FixParameter(0, ((double)d/2)*sqrt(alpha*(1-alpha))/ChiSquaredDist(R2, d)/R2);
      fn[d-mind]->FixParameter(1, -.5);
    }
    gball_rms[d-mind][2]->Fit("fn");
    fn[d-mind]->SetTitle(TString::Format("%0.2fN^{%0.2f}",
	fn[d-mind]->GetParameter(0), fn[d-mind]->GetParameter(1)));
    fn[d-mind]->SetLineColor(gball_rms[d-mind][2]->GetLineColor());;
    fn[d-mind]->Draw("SAME");
  }
  legmean->Draw("SAME");
  c->SaveAs("hull_rms.pdf");
  delete c;

  return 0;
}
