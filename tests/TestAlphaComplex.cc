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
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TGraphErrors.h"

// Additional modules
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

#include "DTFE.hh"
#include "TText.h"
#include "TGraph2D.h"
#include "TColor.h"

TRandom3 rdmzerr;

double ran() {

  return rdmzerr.Uniform(-2, 2);
//  return rdmzerr.Gaus(0.5, 0.1);
}

double prefactor(double d) {

  double num = 2*pow(d, d-2)*pow(M_PI, (d-1)/2)*tgamma((d*d+1)/2)*pow(tgamma(d/2), d);
  double den = (d+1)*tgamma(d*d/2)*pow(tgamma((d+1)/2), d);
  return num/den;
}

/*double gengamma(double* x, double* par) {

  double k = par[0];
  double th = par[1];

  return (1./(tgamma(k)*pow(th, k)))*pow(x[0], k-1)*exp(-pow(x[0],par[2])/th);
}*/

double gengamma(double* x, double* par) {

  double a = par[0];
  double b = par[1];
  double c = par[2];

  return ((a*pow(b, c/a))/(tgamma(c/a)))*pow(x[0], c-1)*exp(-b*pow(x[0], a));
}

double invgengamma(double* x, double* par) {

  double a = par[0];
  double b = par[1];
  double c = par[2];

  return ((a*pow(b, c/a))/(tgamma(c/a)))*pow(x[0], 1.-c)*exp(-b*pow(x[0], -a))/x[0]/x[0];
}

double agengamma(double* x, double* par) {

  double p = 2*par[0];

  return (pow(p, p)/tgamma(p))*pow(x[0], p-1)*exp(-p*x[0]);
}

double ainvgengamma(double* x, double* par) {

  double p = 2*par[0];

  return (pow(p, p)/tgamma(p))*pow(x[0], -p-1)*exp(-p/x[0]);
}

int main(int argc, char ** argv) {
  gStyle->SetLabelSize(0.04, "XY");
  gStyle->SetTitleSize(0.04, "XY");
  gStyle->SetOptStat(0);

  std::vector<double> stops = {0., 0.125, .25, .375, .5, .625, .75, .875, 1.};
  std::vector<double> red = {0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764};
  std::vector<double> green = {0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832};
  std::vector<double> blue = {0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539};
  TColor::CreateGradientColorTable(9, &stops[0], &red[0], &green[0], &blue[0], 255, 1);


  // Create grid
  std::vector<Vertex> simpvertices;
  for (size_t i = 0; i < 20; i ++) {
    double x = rdmzerr.Uniform(-3, 3);
    simpvertices.push_back(Vertex({x}, TMath::Gaus(x, 0, 1)/sqrt(2*M_PI)));
  }

  Interpolator sinterp(simpvertices);

  std::vector<Vertex> simpvertices2;
  std::vector<std::vector<double>> vorvertices;
  for (size_t i = 0; i < 100; i ++) {
    double x = rdmzerr.Uniform(-2, 2);
    double y = rdmzerr.Uniform(-2, 2);
    vorvertices.push_back({ran(),ran()});
    simpvertices2.push_back(Vertex({x, y}, TMath::Gaus(x, 0, 1)*TMath::Gaus(y, 0, 1)/(2*M_PI)));
  }

  Interpolator sinterp2(simpvertices2);
  

  // Create grid
  std::vector<Vertex> intvertices;
  for (size_t i = 0; i < 11; i ++) {
    double x = -3 + (double)i*6./10;
    intvertices.push_back(Vertex({x}, TMath::Gaus(x, 0, 1)/sqrt(2*M_PI)));
  }
  std::vector<Vertex> intvertices2;
  for (size_t i = 0; i < 5; i++) {
    double x = (double)i/4;
    for (size_t j = 0; j < 5; j++) {
      double y = (double)j/4;
      intvertices2.push_back(Vertex({x, y}, ran()));
    }
  }
  Grid grid(intvertices);
  Interpolator interp(grid);

  Grid grid2(intvertices2);
  Interpolator interp2(grid2);

  // Print the interp
  TGraph* ginterp = new TGraph();
  ginterp->SetLineColor(4);
  ginterp->SetLineWidth(2);
  ginterp->SetLineStyle(2);
  for (size_t i = 0; i < 1001; i++) {
    double x = -3 + (double)i*6./1000;
    ginterp->SetPoint(i, x, sinterp(x));
  }

  TF1 *fgaus = new TF1("gaus", "pow(2*TMath::Pi(), -0.5)*TMath::Gaus(x, 0, 1)", -3.5, 3.5);
  TCanvas *cint = new TCanvas("c", "c", 1200, 800);
  fgaus->SetTitle(";x;#rho(x)");
  fgaus->Draw();
  ginterp->Draw("L SAME");
  cint->SaveAs("gaus_lint.pdf");
  delete cint;

  // Print the interp
  TH2F* ginterp2 = new TH2F("ginterp2", ";x;y;#rho(x,y)", 100, -2, 2, 100, -2, 2);
  for (size_t i = 0; i < ginterp2->GetXaxis()->GetNbins(); i++) {
    for (size_t j = 0; j < ginterp2->GetYaxis()->GetNbins(); j++) {
      double x = ginterp2->GetXaxis()->GetBinCenter(i+1);
      double y = ginterp2->GetYaxis()->GetBinCenter(j+1);
      ginterp2->SetBinContent(i+1, j+1, sinterp2({x, y}));
    }
  }

  cint = new TCanvas("c", "c", 1200, 800);
  gPad->SetRightMargin(.15);
  ginterp2->Draw("COLZ");
  for (TPolyLine* poly : sinterp2.Meshing())
	poly->Draw("SAME");
  cint->SaveAs("random_lint.pdf");
  delete cint;

  Voronoi vor(vorvertices, true, 1);
  vor.Draw("", true, true, true);
  cint = new TCanvas("c", "c", 1200, 800);
  delete cint;


  /////////////////////////////////////////////////////////////////////////////////
  std::vector<Vertex> vertices;
  TGraph *gdtfe = new TGraph();
  gdtfe->SetMarkerStyle(21);
  size_t nv = 1000;
  std::vector<TText*> texts(nv);
  for (size_t i = 0; i < nv; i++) {
     vertices.push_back(Vertex({ran(), ran()}, 1));
     gdtfe->SetPoint(i, vertices[i][0], vertices[i][1]);
     texts[i] = new TText(vertices[i][0], vertices[i][1], TString::Format("%d", (int)i));
  }
  DTFE dtfe(vertices);

  TCanvas* cdtfe = new TCanvas("c", "c", 800, 800);
  gdtfe->Draw("AP");
  for (size_t i = 0; i < 10; i++)
    texts[i]->Draw("SAME");
  for (TPolyLine* poly : dtfe.Polygons())
    poly->Draw("SAME");
  cdtfe->SaveAs("dtfe.pdf");
  delete cdtfe;

  TGraph2D* graph2d = new TGraph2D(nv);
  for (size_t i = 0; i < nv; i++)
      graph2d->SetPoint(i, dtfe.GetVertexArray()[i][0], dtfe.GetVertexArray()[i][1],
						dtfe.GetVertexArray()[i].GetValue());
     
  cdtfe = new TCanvas("c", "c", 800, 800);
  graph2d->Draw("PCOL");
  cdtfe->SaveAs("dtfe2.pdf");
  delete cdtfe;

  // Draw the Tanemura functions up to dimension 4, and its inverse

/*  TF1 *ftane1 = new TF1("tane1", gengamma, 0, 5, 3);
  ftane1->SetParameters(1, 1, 1);
  TF1 *ftane2 = new TF1("tane2", gengamma, 0, 5, 3);
  ftane2->SetParameters(0.85, 1.83, 1.50);
  TF1 *ftane3 = new TF1("tane3", gengamma, 0, 5, 3);
  ftane3->SetParameters(0.74, 2.47, 1.81);
  TF1 *ftane4 = new TF1("tane4", gengamma, 0, 5, 3);
  ftane4->SetParameters(0.70, 2.97, 1.97);
  TF1 *ftane5 = new TF1("tane5", gengamma, 0, 5, 3);
  ftane5->SetParameters(0.61, 3.49, 2.08);
  std::vector<TF1*> ftanes = {ftane1, ftane2, ftane3, ftane4, ftane5};
*/
  std::vector<TF1*> ftanes(6);
  for (size_t i = 0; i < 6; i++) {
    ftanes[i] = new TF1(TString::Format("tane%d", (int)i), agengamma, 0, 3, 1);
    ftanes[i]->SetParameter(0, 2*(i+1));
  }

/*  TF1 *finvtane1 = new TF1("invtane1", invgengamma, 0, 5, 3);
  finvtane1->SetParameters(1, 1, 1);
  TF1 *finvtane2 = new TF1("invtane2", invgengamma, 0, 5, 3);
  finvtane2->SetParameters(0.85, 1.83, 1.50);
  TF1 *finvtane3 = new TF1("invtane3", invgengamma, 0, 5, 3);
  finvtane3->SetParameters(0.74, 2.47, 1.81);
  TF1 *finvtane4 = new TF1("invtane4", invgengamma, 0, 5, 3);
  finvtane4->SetParameters(0.70, 2.97, 1.97);
  TF1 *finvtane5 = new TF1("invtane5", invgengamma, 0, 5, 3);
  finvtane5->SetParameters(0.61, 3.49, 2.08);
*/
  std::vector<TF1*> finvtanes(6);
  for (size_t i = 0; i < 6; i++) {
    finvtanes[i] = new TF1(TString::Format("invtane%d", (int)i), ainvgengamma, 0, 3, 1);
    finvtanes[i]->SetParameter(0, 2*(i+1));
  }

  TLegend *legtane = new TLegend(.75, .45, 1., .89);
  legtane->SetFillStyle(0);
  legtane->SetLineColorAlpha(0, 0);

  std::vector<int> colors = {kBlack, kMagenta+1, kBlue+1, kCyan+1, kGreen+2, kOrange+1, kRed+1};
  for (size_t i = 0; i < ftanes.size(); i++) {
    ftanes[i]->SetLineWidth(2);
    finvtanes[i]->SetLineWidth(2);
    if ( i != 0 ) {
        ftanes[i]->SetLineWidth(3);
        finvtanes[i]->SetLineWidth(3);
    }
    ftanes[i]->SetLineColor(colors[i]);
    ftanes[i]->SetLineStyle(i+1);
    finvtanes[i]->SetLineColor(colors[i]);
    finvtanes[i]->SetLineStyle(i+1);
    legtane->AddEntry(ftanes[i], TString::Format("d=%d", (int)i+1), "l");
  }
  
  TCanvas* canv = new TCanvas("c", "c", 1200, 800);
  ftanes[5]->Draw();
  ftanes[5]->SetTitle(";#upsilon; f_{d} (#upsilon)");
  for (size_t i = 0; i < ftanes.size()-1; i++)
      ftanes[i]->Draw("SAME");
  legtane->Draw("SAME");
  TLine* line = new TLine(1, 0, 1, ftanes[5]->GetHistogram()->GetMaximum()*(1+gStyle->GetHistTopMargin()));
  line->SetLineColorAlpha(kBlack, .4);
  line->Draw("SAME");
  canv->SaveAs("del_area_pdf.pdf");
  delete canv;

  canv = new TCanvas("c", "c", 1200, 800);
  finvtanes[5]->Draw();
  finvtanes[5]->SetTitle(";#upsilon^{-1}; f_{d}^{-1} (#upsilon^{-1})");
  for (size_t i = 0; i < ftanes.size()-1; i++)
      finvtanes[i]->Draw("SAME");
  legtane->Draw("SAME");
  TLine* line2 = new TLine(1, 0, 1, finvtanes[5]->GetHistogram()->GetMaximum()*(1+gStyle->GetHistTopMargin()));
  line2->SetLineColorAlpha(kBlack, .4);
  line2->Draw("SAME");
  canv->SaveAs("del_invarea_pdf.pdf");
  delete canv;

  // Draw the Delaunay triangulation of some points
  gStyle->SetOptStat(2210);
  gStyle->SetOptFit(1);
  std::vector<std::vector<double>> s1, s2, s3, s4, s5;
  std::vector<double> point1(1), point(2), point3(3), point4(4), point5(5);
  TRandom3 rdmzerr;
  size_t nnn = 1e4;
  while ( s3.size() < nnn ) {

    double x = ran();
    double y = ran();
    if ( (x-.5)*(x-.5)+(y-.5)*(y-.5) < 0.25 ) {
//    if ( (x-.5)*(x-.5)-2*(x-.5)*(y-.5)+2*(y-.5)*(y-.5) < 0.25 ) {
      point[0] = x;
      point[1] = y;
    } else {
      continue;
    }

//    point[0] = ran();
//    point[1] = ran();
    point1[0] = ran();
    point3[0] = ran();
    point3[1] = ran();
    point3[2] = ran();
    point4[0] = ran();
    point4[1] = ran();
    point4[2] = ran();
    point4[3] = ran();
    point5[0] = ran();
    point5[1] = ran();
    point5[2] = ran();
    point5[3] = ran();
    point5[4] = ran();
    s1.push_back(point1);
    s2.push_back(point);
    s3.push_back(point3);
    s4.push_back(point4);
    s5.push_back(point5);
  }

  // Drift
  for (size_t i = 0; i < s2.size(); i++)
      s2[i][0] += 10*s2[i][1];


  // Voronoi tesselation
  TGraph* graph = new TGraph();
  if ( true ) {
    graph->SetMarkerStyle(20);
    size_t i;
    for (i = 0; i < s2.size(); i++)
	graph->SetPoint(i, s2[i][0], s2[i][1]);
  }
//  gStyle->SetPalette(kGreyScale);
  Voronoi tess2(s2, 1);
//  TCanvas *cv = new TCanvas("c", "c", 1800, 600);
  tess2.Draw("", true, true);
//  graph->Draw("P");
//  for (TPolyLine* poly : tess2.Polygons() )
//      poly->Draw("SAME");
//  cv->SaveAs("vor2.pdf");
//  delete cv;


  

  std::cerr << "n=2: " << prefactor(2) << std::endl;
  std::cerr << "n=3: " << prefactor(3) << std::endl;
  std::cerr << "n=4: " << prefactor(4) << std::endl;
  std::cerr << "n=5: " << prefactor(5) << std::endl;
  std::cerr << "n=6: " << prefactor(6) << std::endl;


  std::cerr << s2[0][0] << "  " << s2[1][0] << std::endl;
  AlphaComplex del2(s2, 0);
//  del2.Paint("del2", false, true);

/*  std::vector<double> areas = del2.GetCellVolumeArray();
  std::vector<double> logareas(areas.size());
  std::vector<double> invareas(areas.size());
  std::vector<double> normareas(areas.size());
  double totarea = 0;
  for (size_t i = 0; i < areas.size(); i++) {
     logareas[i] = log10(areas[i]);
//     normareas[i] = sqrt(areas[i]*nnn*prefactor(2));
     normareas[i] = areas[i]*nnn;
     totarea += areas[i];
     if ( areas[i] > 2*1e-4 )
	std::cerr << areas[i] << std::endl;
  }*/

  std::vector<double> areas;
  Voronoi vorr(s2, false);
  for (const Hull& cell : vorr.GetCellArray())
    if ( cell.GetVolume() >= 0 )
        areas.push_back(cell.GetVolume());

  std::vector<double> logareas(areas.size());
  std::vector<double> invareas(areas.size());
  std::vector<double> normareas(areas.size());
  double totarea = 0;
  for (size_t i = 0; i < areas.size(); i++) {
     logareas[i] = log10(areas[i]);
//     normareas[i] = sqrt(areas[i]*nnn*prefactor(2));
     normareas[i] = areas[i]*nnn*4/M_PI;
     totarea += areas[i];
  }


  std::cerr << "points size: " << s2.size() << std::endl;
  std::cerr << "areas size: " << areas.size() << std::endl;
  std::cerr << "mean area: " << Math::Mean(areas) << std::endl;
  std::cerr << "tota area: " << totarea << std::endl;

//  TF1* fgamma = new TF1("fgamma", gengamma, 0, 5, 3);
//  fgamma->SetParameters(1, 2, 2);
  TF1* fgamma = new TF1("fgamma", "pow([1], [0])*pow(x, [0]-1)*exp(-[1]*x)/TMath::Gamma([0])", 0, 5);
  fgamma->FixParameter(0, 4);
  fgamma->FixParameter(1, 4);
  TF1* flognorm = new TF1("flognorm", "[0]*TMath::Gaus(log(x), [1], [2])/x", 0, 5);
  flognorm->SetParameters(1, 0, 1);

  TH1F* harea = new TH1F("area", ";A", 100, 0, 5);
  harea->FillN(logareas.size(), &normareas[0], NULL);
  harea->Sumw2();
  TCanvas *cc = new TCanvas("c", "c", 1200, 800);
//  gPad->SetLogy();
  harea->Scale(100/harea->GetEntries()/5);
  harea->Draw();
  harea->Fit(fgamma);
  cc->SaveAs("areas.pdf");
  delete cc;

  AlphaComplex del1(s1, 0);
  areas = del1.GetCellVolumeArray();
  logareas.resize(areas.size());
  std::cerr << areas.size() << std::endl;
  normareas.resize(areas.size());
  for (size_t i = 0; i < areas.size(); i++) {
     logareas[i] = log10(areas[i]);
     normareas[i] = areas[i]*areas.size();
  }
  std::cerr << Math::Mean(normareas) << std::endl;
  std::cerr << Math::RMS(normareas) << std::endl;

  TH1F* harea1 = new TH1F("area1", ";A", 100, 0, 5);
  harea1->FillN(logareas.size(), &normareas[0], NULL);
  harea1->Sumw2();
//  harea1->Fit(fgamma);
  cc = new TCanvas("c", "c", 1200, 800);
  harea1->Scale(100/harea1->GetEntries()/5);
  harea1->Draw();
  cc->SaveAs("areas1.pdf");
  delete cc;



  AlphaComplex del3(s3, 0);
//  del3.Paint("del3", false, true);

  areas = del3.GetCellVolumeArray();
  logareas.resize(areas.size());
  std::cerr << areas.size() << std::endl;
  normareas.resize(areas.size());
  for (size_t i = 0; i < areas.size(); i++) {
     logareas[i] = log10(areas[i]);
     normareas[i] = pow(areas[i]*nnn*prefactor(3), 1./3);
  }

  TH1F* harea3 = new TH1F("area3", ";A", 100, 0, 5);
  harea3->FillN(logareas.size(), &normareas[0], NULL);
  harea3->Sumw2();
  cc = new TCanvas("c", "c", 1200, 800);
  harea3->Scale(100/harea3->GetEntries()/5);
  harea3->Draw();
  harea3->Fit(fgamma);
  cc->SaveAs("areas3.pdf");
  delete cc;

  AlphaComplex del4(s4, 0);
  areas = del4.GetCellVolumeArray();
  logareas.resize(areas.size());
  std::cerr << areas.size() << std::endl;
  normareas.resize(areas.size());
  for (size_t i = 0; i < areas.size(); i++) {
     logareas[i] = log10(areas[i]);
     normareas[i] = pow(areas[i]*nnn*prefactor(4), 1./4);
  }


  TH1F* harea4 = new TH1F("area4", ";A", 100, 0, 5);
  harea4->FillN(logareas.size(), &normareas[0], NULL);
  harea4->Sumw2();
  cc = new TCanvas("c", "c", 1200, 800);
  harea4->Scale(100/harea4->GetEntries()/5);
  harea4->Draw();
  harea4->Fit(fgamma);
  cc->SaveAs("areas4.pdf");
  delete cc;

  AlphaComplex del5(s5, 0);
  areas = del5.GetCellVolumeArray();
  logareas.resize(areas.size());
  std::cerr << areas.size() << std::endl;
  normareas.resize(areas.size());
  for (size_t i = 0; i < areas.size(); i++) {
     logareas[i] = log10(areas[i]);
     normareas[i] = pow(areas[i]*nnn*prefactor(5), 1./5);
  }
  std::cerr << Math::Mean(normareas) << std::endl;
  std::cerr << Math::RMS(normareas) << std::endl;

  TH1F* harea5 = new TH1F("area4", ";A", 100, 0, 5);
  harea5->FillN(logareas.size(), &normareas[0], NULL);
  harea5->Sumw2();
  cc = new TCanvas("c", "c", 1200, 800);
  harea5->Scale(100/harea5->GetEntries()/5);
  harea5->Draw();
  harea5->Fit(fgamma);
  cc->SaveAs("areas5.pdf");
  delete cc;


/*  TCanvas* cdel = new TCanvas("c", "c", 800, 800);
//  cdel->DrawFrame(-1, -1, 1, 1);
  for ( TPolyLine* pol : del2.Polygons() ) {
      pol->SetLineWidth(2);
      pol->Draw("SAME");
  }
  cdel->SaveAs("del2.pdf");*/
  

  // Draw the PDFs and CDFs
  std::vector<DFunction*> funcs;
  funcs.push_back(new DTriangular(0, 0.5));
  funcs.push_back(new DExponential(1));
  funcs.push_back(new DCauchy());
  funcs.push_back(new DGaus(1));
  funcs.push_back(new DMaxwell(1));
  funcs.push_back(new DUniform(1));

  std::vector<TGraph*> graphs;
  std::vector<TGraph*> cdfs;
  std::vector<std::string> names =
	{"Triangular (#times0.5)", "Exponential", "Cauchy", "Gaus", "Maxwell", "Uniform"};
  
  TLegend *leg = new TLegend(.65, .6, .89, .89);
  TLegend *legcdf = new TLegend(.65, .15, .89, .44);
  leg->SetFillStyle(0);
  leg->SetLineColorAlpha(0, 0);
  legcdf->SetFillStyle(0);
  legcdf->SetLineColorAlpha(0, 0);
  for (size_t i = 0; i < funcs.size(); i++) {
    funcs[i]->SetRange(0, 3);
    graphs.push_back(funcs[i]->GraphRadial());
    cdfs.push_back(funcs[i]->GraphCDFRadial());
    if ( i != 0 ) {
      graphs[i]->SetLineWidth(3);
      cdfs[i]->SetLineWidth(3);
    }
    graphs[i]->SetLineColor(colors[i]);
    graphs[i]->SetLineStyle(i+1);
    cdfs[i]->SetLineColor(colors[i]);
    cdfs[i]->SetLineStyle(i+1);
    leg->AddEntry(graphs[i], names[i].c_str(), "l");
    if ( i == 0 ) {
      legcdf->AddEntry(cdfs[i], "Triangular", "l");
      continue;
    }
    legcdf->AddEntry(cdfs[i], names[i].c_str(), "l");
  }

  TCanvas* c = new TCanvas("c", "c", 1200, 800);
  graphs[0]->Draw();
  graphs[0]->SetTitle(";R;#rho(R)");
  graphs[0]->GetXaxis()->SetRangeUser(0, 3);
  for (TGraph* graph : graphs)
      graph->Draw("SAME");
  leg->Draw("SAME");
  c->SaveAs("PDF.pdf");
  delete c;

  c = new TCanvas("c", "c", 1200, 800);
  cdfs[0]->Draw();
  cdfs[0]->SetTitle(";R;F(R)");
  cdfs[0]->GetXaxis()->SetRangeUser(0, 3);
  for (TGraph* cdf : cdfs)
      cdf->Draw("SAME");
  legcdf->Draw("SAME");
  c->SaveAs("CDF.pdf");
  delete c;

  TGraph *graph1 = new TGraph();
  graph1->SetTitle("L^{1}");
  graph1->SetFillColor(0);
  graph1->SetLineColor(0);
  graph1->SetMarkerStyle(20);
  graph1->SetMarkerSize(2);
  graph1->SetMarkerColor(colors[0]);
  TGraph *graph2 = new TGraph();
  graph2->SetTitle("L^{2}");
  graph2->SetFillColor(0);
  graph2->SetLineColor(0);
  graph2->SetMarkerStyle(21);
  graph2->SetMarkerSize(2);
  graph2->SetMarkerColor(colors[1]);
  TGraph *graphinf = new TGraph();
  graphinf->SetTitle("L^{#infty}");
  graphinf->SetFillColor(0);
  graphinf->SetLineColor(0);
  graphinf->SetMarkerStyle(22);
  graphinf->SetMarkerSize(2);
  graphinf->SetMarkerColor(colors[2]);

  TMultiGraph* glp = new TMultiGraph("glp", ";n;V_{p}/R^{n}");
  glp->Add(graph1, "P");
  glp->Add(graph2);
  glp->Add(graphinf);

  for (size_t i = 0; i < 6; i++) {
    graph1->SetPoint(i, i+1, pow(2, i+1)/TMath::Gamma(i+2));
    graph2->SetPoint(i, i+1, pow(M_PI, (double)(i+1)/2)/TMath::Gamma((double)(i+1)/2+1));
    graphinf->SetPoint(i, i+1, pow(2, i+1));
  }
  c = new TCanvas("c", "c", 1200, 800);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  glp->Draw("AP");
  gPad->BuildLegend(.15, .15, .25, .3);
  c->SaveAs("contour_volume.pdf");
  delete c;  


  ////////////////////////////////////////////
  //////////// GLOBAL PARAMETERS /////////////
  ////////////////////////////////////////////
  gStyle->SetOptStat(0);
  double alpha = .5;
  bool sorted = false;
  std::string algo = "knn";
  if ( argc < 2 ) {
    std::cerr << "No density estimation algorithm specified, using knn" << std::endl;
  } else {
    algo = argv[1];
  }

  ////////////////////////////////////////////
  ///////////// VISUAL EXAMPLE ///////////////
  ////////////////////////////////////////////
  // Import the points from the "MICE" string samples if available
  gSystem->Load("libTree");
  TFile data_file("samples_MICE.root");
  if ( data_file.IsOpen() ) {
    double rho, R;
    TTree *T = (TTree*)data_file.Get("tree");
    std::vector<double>* v(NULL);
    T->SetBranchAddress("points", &v);
    std::vector<std::vector<double>> points_mice;
    for (size_t i = 0; i < (size_t)T->GetBranch("points")->GetEntries(); i++) {
      T->GetBranch("points")->GetEntry(i);
      points_mice.push_back(*v);
    }

    // Flatten the imported points to 2d and draw their alpha-complex
    std::vector<std::vector<double>> points_mice_2d = points_mice;
    for (size_t i = 0; i < points_mice.size(); i++)
        points_mice_2d[i].resize(2);

    rho = 1./4;
    R = 1.25/sqrt(rho*points_mice.size());
    AlphaComplex alphac(points_mice_2d, 1./R);
    alphac.Paint("MICE", true, true);

    // Draw the 3D alpha-complex of the MICE samples for the optimal radius
    rho = 1./0.8;
    R = 1./pow(rho*points_mice.size(), 1./3);
    AlphaComplex alphac3d(points_mice, 1./R);
    alphac3d.Paint("MICE_3d", false, false, {"pdf", "root"});

    // Close the file
    data_file.Close();
  } else {
    Pitch::print(Pitch::warning, "String samples not found, will proceed with rest of the test");
  }

  ////////////////////////////////////////////
  /////////// TEST DISTRIBUTIONS /////////////
  ////////////////////////////////////////////
  std::vector<DFunction*> functions;
//  functions.push_back(new DGaus(1));
  functions.push_back(new DMultiGaus(1, 2));
//  functions.push_back(new DChiSquared(4));
//  functions.push_back(new DCauchy());
//  functions.push_back(new DExponential(1));
//  functions.push_back(new DUniform(1));
//  functions.push_back(new DTriangular(1));
//  functions.push_back(new DMaxwell(1));

//  functions.push_back(new DGaus(2));
  functions.push_back(new DMultiGaus(2, 3));
//  functions.push_back(new DExponential(2));
//  functions.push_back(new DUniform(2));
//  functions.push_back(new DTriangular(2));
//  functions.push_back(new DMaxwell(2));

//  functions.push_back(new DGaus(3));

//  functions.push_back(new DGaus(4));


  ////////////////////////////////////////////
  //////////// OPTIMIZE RADIUS ///////////////
  ////////////////////////////////////////////
  // Initialize a function to fit the optimal radius with, should go as (N*rho)^{-1/n}
  TF1* ffit = new TF1("ffit", "[0]*pow([1]*x, -1./[2])", 0, 1e6);
  ffit->SetParNames("C", "#rho_{C}", "n");
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);

  // Loop over all the test functions
  double dim, vol, level;
  double frac, tempr, tempvol, temperr;
  std::vector<std::vector<double>> tempsample;
  AlphaComplex tempcomp;
  size_t nmax = 1e5;
  size_t i, j;
  for (i = 0; i < functions.size(); i++) {  
    // Fetch the parameters of the current distribution under test
    dim = functions[i]->Dimension();
    vol = functions[i]->ContourVolume(alpha);
    level = functions[i]->Level(alpha);

    // Initialize a TGraphErrors to fill with the points
    TGraphErrors *graph = new TGraphErrors();
    graph->SetMarkerStyle(21);
    graph->SetLineWidth(2);
    graph->SetTitle(TString::Format("Optimal radius for %s", (functions[i]->Title()).c_str()));

    // Loop over fractions varying from 100% to 5% of the total sample
    for (frac = 1; frac > 0.05-1e-6; frac -= .05) {
      // Sample the distribution and sort if requested
      tempsample.resize(0);
      for (j = 0; j < frac*nmax; j++)
	  tempsample.push_back(functions[i]->RandomVector());
      if ( sorted ) {
        DensityEstimator est(tempsample, algo);
        std::vector<Vertex> vertices(tempsample.size());
        for (j = 0; j < tempsample.size(); j++)
	    vertices[j] = Vertex(tempsample[j], est(tempsample[j]));
        std::sort(vertices.begin(), vertices.end(),
		[] (const Vertex& a, const Vertex& b) { return a.GetValue() < b.GetValue(); });

        tempsample.resize(0);
        for (j = alpha*frac*nmax; j < frac*nmax; j++)
	    tempsample.push_back(vertices[j].GetCoordinates());
      }

      // Compute the alpha-complex and optimize R
      tempr = 1./pow(level*tempsample.size(), 1./dim); // Starting point
      tempcomp = AlphaComplex(tempsample, 1./tempr);
      tempvol = tempcomp.GetVolume();
      size_t ite(0);
      if ( tempvol < vol ) {
        while ( tempvol < vol && ite < 1000 ) {
          tempr += .01*tempr;
          tempcomp.SetAlpha(1./tempr);
          tempvol = tempcomp.GetVolume();
	  ite++;
        }
      } else if ( tempvol > vol ) {
        while ( tempvol > vol && ite < 1000 ) {
          tempr -= .01*tempr;
          tempcomp.SetAlpha(1./tempr);
          tempvol = tempcomp.GetVolume();
	  ite++;
        }
      }

      // If it is 2D and at 10%, draw the meshing for the optimal R
      if ( dim == 2 && fabs(frac-.1) < 1e-6)
	  tempcomp.Paint(functions[i]->Name());

      // If it is 3D and at 5%, draw the meshing for the optimal R
      if ( dim == 3 && fabs(frac-.05) < 1e-6)
	  tempcomp.Paint(functions[i]->Name(), false, false, {"pdf", "root"});

      // Relative distance to true volume is relative error on tempr
      temperr = tempr*fabs(tempvol-vol)/vol;
      graph->SetPoint((int)((frac-.05+1e-6)/.05), tempsample.size(), tempr);
      graph->SetPointError((int)((frac-.05+1e-6)/.05), 0., temperr);
    }

    // Fix the level and dimension for the fit
    ffit->FixParameter(1, level);
    ffit->FixParameter(2, dim);
    graph->Fit(ffit);

    // Draw the graph and the fit
    TCanvas *c = new TCanvas("c", "c", 1200, 800);
    graph->Draw("AP");
    graph->GetXaxis()->SetTitle("N");
    graph->GetYaxis()->SetTitle("R");
    if ( !sorted ) {
      c->SaveAs(TString::Format("%s_optr.pdf", (functions[i]->Name()).c_str()));
    } else {
      c->SaveAs(TString::Format("%s_optr_sorted.pdf", (functions[i]->Name()).c_str()));
    }
    delete c;
  }

  return 0;
}
