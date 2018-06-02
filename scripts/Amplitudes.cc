// Cpp includes
#include <iostream>
#include <vector>
#include <map>

// Additional modules
#include "Globals.hh"
#include "DPoisson.hh"
#include "Extractor.hh"

/** @file  Amplitudes.cc
 *
 *  @brief Reconstructs the amplitude distributions of the beam
 **/

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

  // List of types of amplitudes to consider
  std::vector<std::string> locs = {"up", "down", "scrap"};
  std::map<std::string, std::string> loc_names =
	{{"up","Upstream"}, {"down","Downstream"}, {"scrap","Scraped"}};

  // Initialize the beam objects with the samples
  Pitch::print(Pitch::info, "Computing the amplitude distributions");
  for (const std::string& type : types) {
    for (const size_t& i : streams[type].GetPlaneIds()) {

      // Switch to corrected amplitudes
      if ( globals["corrected"] )
          streams[type][i].SetCorrectedAmplitudes();

      // Switch to corrected amplitudes
      if ( globals["mcd"] )
          streams[type][i].SetMCDAmplitudes();

      // Switch to generalised amplitudes
      if ( globals["generalised"] )
          streams[type][i].SetGeneralisedAmplitudes(1./streams[type].Transmission(i));

      // If the Voronoi volumes are requested, reconstruct them
      if ( globals["voronoi"] )
          streams[type][i].SetVoronoiVolumes();
    }
  }

  // If requested, save the Poincaré sections of the beam along with their amplitdues
  if ( globals["poincare"] ) {
    std::vector<std::string> vars = {"x", "y", "px", "py", "pz"};
    TFile* outfile = new TFile(TString::Format("%s_poincare.root", run_name.c_str()), "RECREATE");
    Beam::Bunch temp;
    for (const std::string& type : types) {
      for (const std::string loc : {"in", "out"}) {
        // Get the sections and save them to a file, move them to a directory
	temp = (loc == "in") ? streams[type].Front() : streams[type].Back();
        temp.SetName(loc+"_"+type);
	for (size_t i = 0; i < 4; i++) {
	  for (size_t j = i+1; j < 4; j++) {
	    TGraph2D* graph = temp.AmplitudeScatter(vars[i], vars[j])->GetGraph2D();
	    graph->Write(graph->GetName());
	  }
	}
      }  
    }
    outfile->Close();
  }

  // Intialize the drawer and the info box
  Beam::Drawer drawer;
  if ( globals["mice"] ) {
    std::string data_type = globals["data"] ? "Preliminary" : "[simulation]";
    std::string maus_version = ext.Get("MausVersion");
    InfoBox* info = new InfoBox(data_type, maus_version,
  	globals["run_name"].AsString(), globals["user_cycle"].AsString());
    drawer.SetInfoBox(info);
  }

  // Show the amplitude distributions and change in the sample
  TFile *outfile = new TFile(TString::Format("%s_amps.root", run_name.c_str()), "RECREATE");
  Pitch::print(Pitch::info, "Filling the amplitude distributions");
  std::map<std::string, THStack*> amp_stacks;
  std::map<std::string, std::map<std::string, TH1F*>> amp_hists;	// Transverse amplitudes
  std::map<std::string, TH2F*> vamp_hists;
  std::map<std::string, TH1F*> ramp_hists, famp_hists;
  size_t nbins(20);
  size_t minamp(0), maxamp(100);
//  size_t minamp(0), maxamp(5e7);
  for (const std::string& type : types) {

    // Initialize the amplitude change histogram and the histogram stack
    vamp_hists[type] = new TH2F(TString::Format("vamp_%s", type.c_str()),
	";A_{#perp}^{U} [mm];A_{#perp}^{D}-A_{#perp}^{U} [mm]",
	nbins, minamp, maxamp, 25, -25, 25);

    ramp_hists[type] = new TH1F(TString::Format("ramp_%s", type.c_str()),
	";A_{#perp}^{U} [mm];N_{S}/N_{U} [%]",
	nbins, minamp, maxamp);

    famp_hists[type] = new TH1F(TString::Format("famp_%s", type.c_str()),
	";A_{#perp}^{U} [mm];(N_{D}-N_{U})/N_{U} [%]",
	nbins, minamp, maxamp);

    amp_stacks[type] = new THStack(TString::Format("amp_stack_%s", type.c_str()),
	";A_{#perp}  [mm]");
//	"Transverse amplitude;V ln^{2}(#rho/#rho(0)) [mm^{2}MeV^{2}/c^{2}]");


    // Fill the amplitude vectors
    std::map<std::string, Vector<double>> amps;
    amps["up"] = Vector<double>(streams[type].Front().Amplitudes());
    amps["down"] = Vector<double>(streams[type].Back().Amplitudes());
    amps["scrap"] = Vector<double>(0);
    amps["through"] = Vector<double>(0);
    std::vector<double> vamps(0);
    size_t i, id(0);
    for (i = 0; i < amps["up"].size(); i++)
      if ( streams[type].GetReaches()[i] != streams[type].GetPlaneIds().back() ) {
	amps["scrap"].push_back(amps["up"][i]);
      } else {
        amps["through"].push_back(amps["up"][i]);
        vamps.push_back(amps["down"][id]-amps["up"][i]);
	id++;
      }

    // Initialize the amplitude histograms
    for (const std::string& loc : locs) {
      if ( !globals["rebin"] ) {
        amp_hists[type][loc] = new TH1F(TString::Format("amp_%s_%s", loc.c_str(), type.c_str()),
				        "", nbins, minamp, maxamp);
      } else {
        // Get the edges of bins of equal size in R^4 (A_perp is R^2)
        std::vector<double> edges(nbins+1);
        for (i = 0; i < nbins+1; i++)
            edges[i] = minamp + pow((double)i/nbins, .5)*(maxamp-minamp);
        amp_hists[type][loc] = new TH1F(TString::Format("amp_%s_%s", loc.c_str(), type.c_str()),
				        "", nbins, &edges[0]);        
      }

      amp_hists[type][loc]->FillN(amps[loc].size(), &(amps[loc][0]), NULL);
      amp_stacks[type]->Add(amp_hists[type][loc]);
    }

    vamp_hists[type]->FillN(amps["through"].size(), &(amps["through"][0]), &(vamps[0]), NULL, 1);
  }

  // Plot the true and reconstructed amplitude evolution
  std::vector<int> colors = {kBlack, kBlue+1, kRed+1};
  std::vector<int> fills = {0, 1000, 3144};
  Pitch::print(Pitch::info, "Printing the amplitude distributions");
  size_t id;
  for (const std::string& type : types) {
    TCanvas* canv = new TCanvas("c", "c", 1200, 800);
    TLegend* leg = new TLegend(.7, .5, .89, .69);
//    TLegend* leg = new TLegend(.11, .55, .3, .7);
    leg->SetLineColorAlpha(0, 0);
    leg->SetFillStyle(3001);
    leg->SetFillColorAlpha(0, 0);

    id = 0;
    for (const std::string& loc : locs) {
      leg->AddEntry(amp_hists[type][loc],
	TString::Format("%s", loc_names[loc].c_str()), "f");
//	TString::Format("%s: %d", loc_names[loc].c_str(), (int)amp_hists[type][loc]->GetEntries()), "ple");
      amp_hists[type][loc]->SetLineWidth(2);
      amp_hists[type][loc]->SetLineColor(colors[id]);
      amp_hists[type][loc]->SetFillStyle(fills[id]);
      amp_hists[type][loc]->SetFillColorAlpha(colors[id]+1, .5);
//      amp_hists[type][loc]->SetMarkerStyle(24+id);
//      amp_hists[type][loc]->SetMarkerSize(1.5);
      amp_hists[type][loc]->Write(amp_hists[type][loc]->GetName());
      id++;
    }

    // Set the bin labels according to the excess and significance
    if ( globals["significance"] ) {
      int nup, ndown, ntotup(0), ntotdown(0);
      double frac, pval, signi;
      for (size_t i = 0; i < (size_t)amp_hists[type]["up"]->GetNbinsX(); i++) {
        nup = amp_hists[type]["up"]->GetBinContent(i+1);
        ntotup += nup;
        ndown = amp_hists[type]["down"]->GetBinContent(i+1);
        ntotdown += ndown;
        frac = (double)(ndown-nup)/nup;
        DPoisson pdf(nup);
        pval= 1.-pdf.CDF(ndown);
        signi = sqrt(2)*TMath::ErfInverse(2*TMath::Gamma(nup+1, ndown)-1);
        amp_hists[type]["up"]->GetXaxis()->SetBinLabel(i+1,
	TString::Format("#splitline{#splitline{#DeltaN: %d}{%2.1f %%}}{%2.1f#sigma}",
			ndown-nup, 100*frac, signi));
      }

      frac = (double)(ntotdown-ntotup)/ntotup;
      DPoisson pdf(ntotup);
      pval= 1.-pdf.CDF(ntotdown);
      signi = sqrt(2)*TMath::ErfInverse(1.-2*pval);
      Pitch::print(Pitch::info, "Total excess: "+std::to_string(100*frac)+" %");
      Pitch::print(Pitch::info, "Exess p-value: "+std::to_string(pval));
      Pitch::print(Pitch::info, "Exess significance: "+std::to_string(signi)+" sigma");
    }

    // Draw everything
    amp_stacks[type]->Draw("HISTE NOSTACK");
    drawer.GetInfoBox()->Draw();

    leg->Draw("SAME");
    canv->SaveAs(std::string("amps_"+type+".pdf").c_str());
    delete leg;
    delete canv;

    // Draw the amplitude change distribution
    gStyle->SetOptStat(0);
    canv = new TCanvas("c", "c", 1200, 800);
//    gPad->SetLogx();
    vamp_hists[type]->Write(vamp_hists[type]->GetName());
    vamp_hists[type]->Draw("COLZ");
    TH1F* vamp_x = drawer.ProjectionMeanX(vamp_hists[type]);
    vamp_x->SetLineColor(2);
    vamp_x->SetLineWidth(3);
    vamp_x->SetMarkerStyle(21);
    vamp_x->Write(vamp_x->GetName());
    vamp_x->Draw("P SAME");
 
    drawer.GetInfoBox()->Draw();
    canv->SaveAs(std::string("vamps_"+type+".pdf").c_str());
    delete vamp_x;
    delete canv;

    // Draw the relative amplitude change plot
    gStyle->SetOptStat(0);
    canv = new TCanvas("c", "c", 1200, 800);
    amp_hists[type]["down"]->Sumw2();
    amp_hists[type]["up"]->Sumw2();
    famp_hists[type]->Add(amp_hists[type]["down"], amp_hists[type]["up"], 1, -1);
    famp_hists[type]->Divide(famp_hists[type], amp_hists[type]["up"]);
    famp_hists[type]->SetLineWidth(2);
    famp_hists[type]->Draw("PE");
    famp_hists[type]->Write(famp_hists[type]->GetName());

    drawer.GetInfoBox()->Draw();
    canv->SaveAs(std::string("famps_"+type+".pdf").c_str());
    delete canv;

    // Draw the transmission plot
    gStyle->SetOptStat(0);
    canv = new TCanvas("c", "c", 1200, 800);
    gPad->SetLogy();
    amp_hists[type]["scrap"]->Sumw2();
    ramp_hists[type]->Divide(amp_hists[type]["scrap"], amp_hists[type]["up"]);
    ramp_hists[type]->SetLineWidth(2);
    ramp_hists[type]->Write(ramp_hists[type]->GetName());
    ramp_hists[type]->Draw("PE");

    drawer.GetInfoBox()->Draw();
    canv->SaveAs(std::string("ramps_"+type+".pdf").c_str());
    delete canv;
  }
  outfile->Close();

  Pitch::print(Pitch::info, "Moving the amplitude graphs to "+run_name+"/amps");
  std::string sysCmd = "mkdir -p "+run_name+"/amps; mv *amps*.pdf *_amps.root "+run_name+"/amps; done";
  if ( std::system(sysCmd.c_str()) )
      Pitch::print(Pitch::error, "Couldn't move amplitude plots");
}
