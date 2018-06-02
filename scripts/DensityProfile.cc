// Cpp includes
#include <iostream>
#include <vector>
#include <map>

// Additional modules
#include "Globals.hh"
#include "Extractor.hh"

/** @file  DensityProfile.cc
 *
 *  @brief Produces non-parametric density profiles of the beam
 *
 *  	   Algorithm that produces non-parametric density profiles of the beam upstream and
 *   	   downstream of the MICE absorber. It uses the requested density estimator to reconstruct
 *	   the minimum local density as a function of the fraction of the beam included.
 **/

double therr(double* x, double* par) {

  DGaus gaus4(4);
  double R = gaus4.Radius(x[0]);
  return sqrt(exp(R*R/2)*(1-x[0])/(x[0]*par[0]));
}

/** @brief	Main function
 *
 *  @param	argc		Number of command line arguments
 *  @param	argv		List of command line aguments
 */
int main(int argc, char ** argv) {

  // Reroot the info messages to the terminal for now
  Pitch::setAnOutput(Pitch::info, std::cerr);

  // Amplitude algorithm global parameters (stored in globals)
  Globals &globals = Globals::GetInstance(argc, argv);

  // Types and plane ids of the data to extract. Empty means import all (recmc and data)
  size_t inid((size_t)globals["tku_vid"]), outid((size_t)globals["tkd_vid"]);
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
  std::map<std::string, Beam::Stream> streams = ext.GetStreams(types, plane_ids);

  // Test the reconstruction with an ideal gaussian beam
  // Calculate the LiH parameters in MICE
/*  double depth = 6.5;		// [cm], thickness of the absorber
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

  Transport transport(4);
  if ( true )
      transport.AddDriftSpace(1000, 140);
  if ( false )
      transport.AddSolenoid(1, 0.5);

  in["data"] = GaussianBunch(105.66, 140, 1e4, 10, 500, 0.);
  out["data"] = in["data"];

  eloss.IonizeBunch(out["data"], 140, 105.66);
  scat.ScatterBunch(out["data"], 140, 105.66);
  transport.TransportBunch(out["data"]);*/

  /*std::vector<std::string> tvars = {"px", "py", "x", "y"};
  THStack *htest = new THStack();
  TH1F* hres[39];
  for (size_t i = 0; i < 39; i++) {
     hres[i] = new TH1F(TString::Format("res%d", (int)i), "", 100, 0, 1);
     htest->Add(hres[i]);
  }
  TMultiGraph* gtest = new TMultiGraph();
  for (size_t ite = 0; ite < 1e3; ite++) {
    std::cerr << ite << std::endl;
    size_t d = 2;
//    std::map<std::string, std::vector<double>> samples = GaussianBunch(105.66, 140, 1e4, 10, 500, 0.);
    int time = std::chrono::duration_cast<std::chrono::nanoseconds>			
			(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
    TRandom3 rdmzer(time);
    std::map<std::string, std::vector<double>> samples;
    for (size_t i = 0; i < 1e4; i++) {
	samples["x"].push_back(rdmzer.Gaus(0, 1));
//	samples["y"].push_back(rdmzer.Gaus(0, 1));
	samples["px"].push_back(rdmzer.Gaus(0, 1));
//	samples["py"].push_back(rdmzer.Gaus(0, 1));
    }

//    Bunch(samples).CovarianceMatrix().Print();
    double scale = 1.;
//    RescaleBunch(samples, scale);

    // Convert the data to simple arrays to be fed to the estimator
    std::vector<std::string> tvars = {"px", "x"};
    double size = samples["x"].size();
    std::vector<std::vector<double>> points(size);
    std::vector<double> point(d);
    for (size_t i = 0; i < size; i++) {
      for (size_t j = 0; j < d; j++)
	  point[j] = samples[tvars[j]][i];
      points[i] = point;
    }

    // Initialiaze the estimator
    DensityEstimator de(points, "knn", false, 1./scale);

    // Compute the density in the points, find the level of the relevant quantile
    std::vector<double> levels(size);
    for (size_t i = 0; i < size; i++)
	levels[i] = de(points[i]);
    std::sort(levels.rbegin(), levels.rend());
    
    double alpha, level;
    for (size_t i = 0; i < 39; i++) {
      alpha = 0.025+i*0.025;
      level = levels[alpha*(size-1)];
      hres[i]->Fill(level);
    }

    // Fill a graph, add it to the list
    size_t id = 0;
    TGraph* g = new TGraph();
    g->SetLineColorAlpha(kBlack, .05);
    g->SetLineWidth(2);
    size_t step = 100;
    for (size_t i = 0; i < size; i += step) {
      alpha = (i+1)*1./(size+1);
      g->SetPoint(id, alpha, levels[i]);
      id++;
    }
    gtest->Add(g, "L");
  }

  TCanvas *ct = new TCanvas("c", "c", 1200, 800);
  htest->Draw("NOSTACK");
  ct->SaveAs("res.pdf");
  delete ct;

  TGraphErrors* gfinal = new TGraphErrors();
  TF1* ferr = new TF1("ferr", therr, 0.01, 0.99, 1);
  ferr->SetParameter(0, 1e4);
  gfinal->SetLineColor(kBlack);
  gfinal->SetLineWidth(2);
  gfinal->SetFillStyle(1000);
  gfinal->SetFillColorAlpha(kBlack, .25);
  for (size_t i = 0; i < 39; i++) {
    double alpha = 0.025+i*0.025;
    double ratio = hres[i]->GetRMS()/hres[i]->GetMean();
    double err = ratio*sqrt(pow(hres[i]->GetRMSError()/hres[i]->GetRMS(), 2) +
				pow(hres[i]->GetMeanError()/hres[i]->GetMean(), 2));
    gfinal->SetPoint(i, alpha, ratio);
    gfinal->SetPointError(i, 0, err);
    std::cerr << alpha << ": " << hres[i]->GetMean() << "  " << hres[i]->GetRMS() << "  " << ratio << std::endl;
  }

  TFile* ofile = new TFile("gfinal.root", "RECREATE");
  gfinal->Write("graph");
  ofile->Close();

  ct = new TCanvas("c", "c", 1200, 800);
  gfinal->Draw("ALE3");
  ferr->Draw("SAME");
  ct->SaveAs("resg.pdf");
  delete ct;

  ct = new TCanvas("c", "c", 1200, 800);
  gtest->Draw("A");
  ct->SaveAs("resmultig.pdf");
  delete ct;*/

  // Initialize the density estimators
  Pitch::print(Pitch::info, "Initializing the estimators");
  std::vector<std::string> tvars = {"x", "y", "px", "py"};
  Beam::Bunch set;
  std::map<std::string, std::map<std::string, std::vector<std::vector<double>>>> points;
  std::map<std::string, std::map<std::string, DensityEstimator>> des;
  std::map<std::string, std::map<std::string, double>> scales;
  std::map<std::string, std::map<std::string, double>> trans;
  std::vector<double> point(4);
  size_t size;
//  int nn, k;
  for (const std::string& type : types) {
    // Get the transmission
    trans[type]["in"] = 1.;
    trans[type]["out"] = (double)streams[type].Back().Size()/streams[type].Front().Size();

    // Get the amount of neighbours from the size of the upstream sample
//    nn = streams[type].Front().Size();
//    k = 1.584*pow(nn, 0.45);

    for (const std::string loc : {"in", "out"}) {
      // Rescale the beam to the normal
      set = (loc == "in") ? streams[type].Front() : streams[type].Back();

      // Convert the data to simple arrays to be fed to the estimator
      size = set.Size();
      points[type][loc].resize(size);
      for (size_t i = 0; i < size; i++) {
 	for (size_t j = 0; j < 4; j++)
	    point[j] = set.Samples(tvars[j])[i];
        points[type][loc][i] = point;
      }

      // Initialiaze the estimator
      des[type][loc] = DensityEstimator(points[type][loc], "knn",
							false, trans[type][loc]);
      des[type][loc].SetName(loc+"_"+type);
    }
  }

  // Show the density profiles and change in the sample
  TFile *outfile = new TFile(TString::Format("%s_de.root", run_name.c_str()), "RECREATE");
  std::map<std::string, std::map<std::string, std::vector<double>>> levels;
  std::map<std::string, std::map<std::string, TGraphErrors*>> profiles;
  std::map<std::string, TGraphErrors*> ratios;
  double alpha, err, ratio;
  size_t id, step;
  TMultiGraph* mg = new TMultiGraph();
  mg->SetTitle(";Fraction #alpha;#rho_{#alpha}  [mm^{-2}(MeV/c)^{-2}]");
  for (const std::string& type : types) {

    step = points[type]["in"].size()/1e3;
    for (const std::string loc : {"in", "out"}) {  
      // Evaluate the density in all of the training points
      Pitch::print(Pitch::info, "Profiling "+loc+" DE for "+type);
      size = points[type][loc].size();
      levels[type][loc].resize(size);
      for (size_t i = 0; i < size; i++)
	  levels[type][loc][i] = des[type][loc](points[type][loc][i]);

      // Sort them by decreasing order of density
      std::sort(levels[type][loc].rbegin(), levels[type][loc].rend());

      // Fill the graph with the levels as a function of the sample fraction	
      id = 0;
      profiles[type][loc] = new TGraphErrors();
      profiles[type][loc]->SetName(TString::Format("de_profile_%s_%s", type.c_str(), loc.c_str()));
      for (size_t i = 0; i < size; i += step) {
	alpha = (i+1)*trans[type][loc]/(size+1);
	profiles[type][loc]->SetPoint(id, alpha, levels[type][loc][i]);
	err = pow(alpha, -1./3)*pow(levels[type][loc][i]/levels[type][loc][0], -.25)/sqrt(size);
	profiles[type][loc]->SetPointError(id, 0., err*levels[type][loc][i]);
	id++;
      }
      profiles[type][loc]->SetPoint(id, trans[type][loc], 0.);

      // Add them to the multigraph to draw, write them to the file to be save
      mg->Add(profiles[type][loc], "LE3");
      profiles[type][loc]->Write(profiles[type][loc]->GetName());
      profiles[type][loc]->SetLineWidth(2);
      profiles[type][loc]->SetFillStyle(1000);
    }

    // Set the style
    profiles[type]["in"]->SetLineColor(kBlack);
    profiles[type]["in"]->SetFillColorAlpha(kBlack, .5);
    profiles[type]["out"]->SetLineColor(kBlue+1);
    profiles[type]["out"]->SetFillColorAlpha(kBlue+2, .5);

    // Draw
    TCanvas *c = new TCanvas("c", "c", 1200, 800);
    mg->Draw("A");
    mg->SetMinimum(0);
    mg->GetXaxis()->SetLimits(0., 1.);
    c->SaveAs(TString::Format("de_profile_%s.pdf", type.c_str()));
    delete c;

    // Draw the ratio of the output over the input
    ratios[type] = new TGraphErrors();
    ratios[type]->SetLineColor(kBlue+1);
    ratios[type]->SetLineWidth(2);
    ratios[type]->SetFillStyle(1000);
    ratios[type]->SetFillColorAlpha(kBlue+2, .5);
    ratios[type]->SetTitle(";Fraction #alpha;#rho_{#alpha}^{u}/#rho_{#alpha}^{d}");
    for (size_t i = 0; i < (size_t)profiles[type]["in"]->GetN(); i++) {
      alpha = profiles[type]["in"]->GetX()[i];
      if ( i+1 > (size_t)profiles[type]["out"]->GetN() || !profiles[type]["in"]->GetY()[i] ) {
	ratios[type]->SetPoint(i, alpha, 0.);
	continue;
      }

      ratio = profiles[type]["out"]->GetY()[i]/profiles[type]["in"]->GetY()[i];
      ratios[type]->SetPoint(i, alpha, ratio);
      if ( !profiles[type]["out"]->GetY()[i] )
	  continue;

      err = ratio*sqrt(pow(profiles[type]["in"]->GetEY()[i]/profiles[type]["in"]->GetY()[i], 2) +
			pow(profiles[type]["out"]->GetEY()[i]/profiles[type]["out"]->GetY()[i], 2));
      ratios[type]->SetPointError(i, 0., err);
    }

    c = new TCanvas("c", "c", 1200, 800);
    ratios[type]->Draw("ALE3");
    ratios[type]->SetMinimum(0);
    ratios[type]->SetMaximum(1.5);
    ratios[type]->GetXaxis()->SetLimits(0., 1.);
    TLine *baseline = new TLine(0, 1., 1., 1.);
    baseline->SetLineColor(kBlack);
    baseline->SetLineWidth(2);
    baseline->Draw("SAME");
/*    TLine *line = new TLine(0, 1./(cool*cool), 1., 1./(cool*cool));
    line->SetLineColor(kRed);
    line->SetLineWidth(2);
    line->Draw("SAME");*/
    c->SaveAs(TString::Format("de_ratio_%s.pdf", type.c_str()));
    delete c;
  }
  outfile->Close();

  if ( globals["poincare"]  ) {
    Pitch::print(Pitch::info, "Saving the Poincaré sections");
    TFile* outfile = new TFile(TString::Format("%s_poincare_de.root", run_name.c_str()), "RECREATE");
    for (const std::string& type : types) {
      for (const std::string loc : {"in", "out"}) {
        // Get the sections and save them to a file, move them to a directory
	for (size_t i = 0; i < 4; i++) {
	  for (size_t j = i+1; j < 4; j++) {
            Pitch::print(Pitch::info, "Saving "+loc+" "+tvars[i]+tvars[j]+" for "+type);
	    TH2F* graph = des[type][loc].Graph2D(-100, 100, -100, 100, i, j, {0, 0, 0, 0});
	    graph->Write(graph->GetName());
	  }
	}
      }  
    }
    outfile->Close();
  }

  Pitch::print(Pitch::info, "Moving the density graphs to "+run_name+"/de");
  std::string sysCmd = "mkdir -p "+run_name+"/de; mv *de*.pdf *_de.root "+run_name+"/de; done";
  if ( std::system(sysCmd.c_str()) )
      Pitch::print(Pitch::error, "Couldn't move density plots");
}
