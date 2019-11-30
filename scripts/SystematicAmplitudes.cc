// Cpp includes
#include <iostream>
#include <vector>
#include <map>

// Additional modules
#include "Globals.hh"
#include "DPoisson.hh"
#include "Extractor.hh"

/** @file  SystematicErrors.cc
 *
 *  @brief Evaluates systematics errors in the MICE data.
 *
 *	   Algorithm that evaluates systematic errors to the MICE data provided with
 *	   with an array of MAUS simulation of the beam.
 **/

/** @brief	Main function
 *
 *  @param	argc		Number of command line arguments
 *  @param	argv		List of command line aguments
 */
int main(int argc, char ** argv) {

  // Reroot the info messages to the terminal for now
  Pitch::setAnOutput(Pitch::info, std::cerr);

  // Systematics algorithm global parameters (stored in globals)
  Globals &globals = Globals::GetInstance(argc, argv);

  // Get the run name from the first file name (remove superfluous with regex)
  std::string run_name = std::regex_replace(globals.GetDataFiles()[0],
				std::regex("import_|(.*/201)|(/tk.*)"), std::string(""));
  run_name.insert(0, "201");
  Pitch::print(Pitch::info, "Run name: "+run_name);

  // Set up the requested plane IDs
  std::map<std::string, size_t> vids = {{"tku", globals["tku_vid"]}, {"tkd", globals["tkd_vid"]}};
  std::map<std::string, size_t> stids = {{"tku", 0}, {"tkd", 5}};
  std::vector<std::string> types = {"utruth", "truth", "recmc"};
  std::vector<std::string> dets;
  std::map<std::string, std::vector<size_t>> plane_ids;
  std::map<std::string, std::vector<size_t>> default_plane_ids = 
	{{"utruth",{vids["tku"], vids["tkd"]}}, {"truth",{vids["tku"], vids["tkd"]}}, {"recmc",{0, 5}}};

  // Extract the types of systematics
  size_t nsets = globals.GetDataFiles().size();
  std::vector<std::string> set_names(nsets);
  std::string set_name;
  for (size_t i = 0; i < nsets; i++) {
    set_name = globals.GetDataFiles()[i];
    set_name = std::regex_replace(set_name, std::regex("(.*/tk)"), std::string(""));
    set_name = std::regex_replace(set_name, std::regex("(/.*)"), std::string(""));
    set_names[i] = "tk"+set_name;
  }

  // Loop over the simulation files, evaluate and store the corrections for each
  // amplitude bin upstream and downstream of the absorber
  TFile *outfile = new TFile(TString::Format("%s_syst.root", run_name.c_str()), "RECREATE");
  std::map<std::string, std::map<std::string, std::vector<double>>> amps;
  std::map<std::string, TH2F*> migrations;
  std::map<std::string, TH1F*> efficiencies;
  std::vector<std::map<std::string, TH1F*>> corrections(nsets);
  std::string set_det;
  for (size_t iFile = 0; iFile < nsets; iFile++) {
    // Set up the extractor
    Beam::Extractor ext({globals.GetDataFiles()[iFile]});

    // Get the setting, determine which stations are required
    set_name = set_names[iFile];
    set_det = set_name.substr(0, 3);
    set_name = set_name.substr(4);
    Pitch::print(Pitch::info, "Importing the "+set_name+" systematics set for "+set_det);

    if ( set_name != "base" ) {
      plane_ids = {{"utruth",{vids[set_det]}}, {"truth",{vids[set_det]}}, {"recmc",{stids[set_det]}}};
      dets = {set_det};
    } else {
      plane_ids = default_plane_ids;
      set_det = "";
      dets = {"tku", "tkd"};
    }

    // Extract the upstream and downstream stations only
    std::map<std::string, Beam::Stream> streams = ext.GetStreams(types, plane_ids);

    // Get the amplitudes
    for (const std::string& type : types) {
      if ( globals["corrected"] )
        for (const size_t& i : streams[type].GetPlaneIds())
            streams[type][i].SetCorrectedAmplitudes();

      for (size_t i = 0; i < dets.size(); i++)
          amps[type][dets[i]] = streams[type][plane_ids[type][i]].Amplitudes();
    }

    // Compute the migration matrices and efficiency vectors
    for (const std::string det : dets) {
      // Fill the migration matrix
      if ( migrations[det] )
	  delete migrations[det];
      migrations[det] = new TH2F(TString::Format("syst_mig_%s_%s", det.c_str(), set_name.c_str()),
			    	 ";A_{#perp}^{r | t}  [mm];A_{#perp}^{r | r}  [mm]",
	  		    	 20, 0, 100, 20, 0, 100);
      migrations[det]->FillN(amps["truth"][det].size(),
		&amps["truth"][det][0], &amps["recmc"][det][0], NULL, 1);
      migrations[det]->Sumw2();

      // Normalize each horizontal bin to 1 
      std::vector<size_t> nentries(20);
      for (size_t j = 0; j < (size_t)migrations[det]->GetNbinsY(); j++) {
	nentries[j] = 0;
        for (size_t i = 0; i < (size_t)migrations[det]->GetNbinsX()+1; i++) // Include outliers
            nentries[j] += (size_t)migrations[det]->GetBinContent(i+1, j+1);
      }

      for (size_t j = 0; j < (size_t)migrations[det]->GetNbinsY(); j++)
        if ( nentries[j] )
          for (size_t i = 0; i < (size_t)migrations[det]->GetNbinsY(); i++)
              migrations[det]->SetBinContent(i+1, j+1,
			migrations[det]->GetBinContent(i+1, j+1)/nentries[j]);

      // Compute the efficiency vectors
      TH1F* amps_rt = new TH1F(TString::Format("syst_eff_tr_%s_%s",
					det.c_str(), set_name.c_str()),	"", 20, 0, 100);
      efficiencies[det] = new TH1F(TString::Format("syst_eff_%s_%s", det.c_str(), set_name.c_str()),
			  	   ";True A_{#perp}  [mm];N^{t | t}/N^{r | t}", 20, 0, 100);
      amps_rt->FillN(amps["truth"][det].size(), &amps["truth"][det][0], NULL);
      efficiencies[det]->FillN(amps["utruth"][det].size(), &amps["utruth"][det][0], NULL);

      // Compute the ratio
      amps_rt->Sumw2();
      efficiencies[det]->Sumw2();
      efficiencies[det]->Divide(amps_rt);

      // For each bin, compute the full correction factor
      TH1F* amps_rr = new TH1F(TString::Format("syst_eff_rr_%s_%s",
						det.c_str(), set_name.c_str()), "", 20, 0, 100);
      amps_rr->FillN(amps["recmc"][det].size(), &amps["recmc"][det][0], NULL);
      corrections[iFile][det] = new TH1F(TString::Format("syst_corr_%s_%s",
				det.c_str(), set_name.c_str()),
				";True A_{#perp}  [mm];N_{#perp}  /#tilde{N}_{#perp}", 20, 0, 100);

      double Ntt, Nrt, Nrr;
      for (size_t i = 0; i < (size_t)migrations[det]->GetNbinsX(); i++) {	
        Nrt = 0.;
        for (size_t j = 0; j < (size_t)migrations[det]->GetNbinsY(); j++)
	    Nrt += migrations[det]->GetBinContent(i+1, j+1)*amps_rr->GetBinContent(j+1);

	Nrr = amps_rr->GetBinContent(i+1);
	Ntt = efficiencies[det]->GetBinContent(i+1)*Nrt;

	corrections[iFile][det]->SetBinContent(i+1, Ntt/Nrr);
	corrections[iFile][det]->SetBinError(i+1, efficiencies[det]->GetBinError(i+1));
      }

      // Set the style of the histograms
      efficiencies[det]->SetMarkerStyle(20);
      efficiencies[det]->SetLineWidth(2);
      efficiencies[det]->SetLineColor(kBlack);

      corrections[iFile][det]->SetMarkerStyle(20);
      corrections[iFile][det]->SetLineWidth(2);
      corrections[iFile][det]->SetLineColor(kBlack);

      // If the run corresponds to the baseline, draw the histograms
      if ( set_name == "base" ) {
	// Migration matrix
        TCanvas* c = new TCanvas("c", "c", 1200, 800);
        migrations[det]->Draw("COLZ");
        c->SaveAs(TString::Format("syst_mig_%s.pdf", det.c_str()));
        delete c;

        // Efficiency histogram
        c = new TCanvas("c", "c", 1200, 800);
        efficiencies[det]->Draw("EP");
        efficiencies[det]->SetMinimum(0.8);
        efficiencies[det]->SetMaximum(1.2);

        TLine *line = new TLine(0, 1, 100, 1);
        line->SetLineColorAlpha(kBlack, .5);
        line->Draw("SAME");

        c->SaveAs(TString::Format("syst_eff_%s.pdf", det.c_str()));
        delete c;

	//Correction histogram
        c = new TCanvas("c", "c", 1200, 800);
        corrections[iFile][det]->Draw("EP");
        corrections[iFile][det]->SetMinimum(0.8);
        corrections[iFile][det]->SetMaximum(1.2);

        line->Draw("SAME");

        c->SaveAs(TString::Format("syst_corr_%s.pdf", det.c_str()));
        delete c;
      }

      // Save histograms to file
      outfile->cd();
      migrations[det]->Write(migrations[det]->GetName());
      efficiencies[det]->Write(efficiencies[det]->GetName());
      corrections[iFile][det]->Write(corrections[iFile][det]->GetName());
      corrections[iFile][det]->SetDirectory(0);
    }
  } // End of the list of files

  outfile->Close();

  // Loop over the systematics types, add the corrections to a stack
  dets = {"tku", "tkd"};
  std::map<std::string, int> colors = {{"base", kBlack}, {"pos_plus", kRed+1},
	{"rot_plus", kGreen+2}, {"scale_C_plus",kMagenta+1}, {"scale_E1_plus", kCyan+1},
        {"scale_E2_plus", kYellow+1}, {"density_plus", kGray}};
  std::map<std::string, int> markers = {{"base", 24}, {"pos_plus", 25},
	{"rot_plus", 26}, {"scale_C_plus",20}, {"scale_E1_plus", 21},
        {"scale_E2_plus", 22}, {"density_plus", 28}};
  std::map<std::string, std::string> names = {{"base", "Baseline"}, {"pos_plus", "+/- 1mm"},
	{"rot_plus", "+/- 1mrad"}, {"scale_C_plus", "CC +/- 1%"}, {"scale_E1_plus", "E1 +/- 5%"},
        {"scale_E2_plus", "E2 +/- 1%"}, {"density_plus", "Density +/- 25%"}};
  std::map<std::string, TLegend*> legends;
  for (const std::string& det : dets) {
    legends[det] = new TLegend(.1, .7, .3, .9);
  }

  // Extract the baseline, throw if absent
  std::map<std::string, TH1F*> hbase;
  std::map<std::string, TH1F*> hfull;
  for (size_t i = 0; i < nsets; i++)
    if ( set_names[i].find("base") != std::string::npos ) {
      hbase["tku"] = corrections[i]["tku"];
      hbase["tkd"] = corrections[i]["tkd"];
      break;
    }

  for (const std::string det : dets)
      hfull[det] = new TH1F(TString::Format("full_%s", det.c_str()), "", 20, 0, 100);

  if ( !hbase["tku"] || !hbase["tkd"] ) {
    Pitch::print(Pitch::error, "Baseline missing, abort");
    return 1;
  }

  // Extract the influence of the variations
  std::map<std::string, std::vector<TH1F*>> stacks;
  double error;
  for (size_t i = 0; i < nsets; i++) {
    set_name = set_names[i];
    set_name = set_name.substr(4);

    // Leave the baseline alone
    if ( set_name.find("base") != std::string::npos )
	continue;

    // Set the style
    std::string det = (set_names[i].find("tku") != std::string::npos) ? "tku" : "tkd";
    corrections[i][det]->SetLineColor(colors[set_name]);
    corrections[i][det]->SetMarkerStyle(markers[set_name]);
    corrections[i][det]->SetMarkerSize(2);
    corrections[i][det]->Divide(hbase[det]);

    // Scale to the baseline
    for (size_t j = 0; j < (size_t)corrections[i][det]->GetNbinsX(); j++)
	corrections[i][det]->SetBinContent(j+1, fabs(corrections[i][det]->GetBinContent(j+1)-1.));

    // Increment the total error
    for (size_t j = 0; j < (size_t)hfull[det]->GetNbinsX(); j++) {
      error = pow(hfull[det]->GetBinContent(j+1), 2);
      error += pow(corrections[i][det]->GetBinContent(j+1), 2);
      hfull[det]->SetBinContent(j+1, sqrt(error));
    }

    // Stack the summary plot
    stacks[det].push_back(corrections[i][det]);
    legends[det]->AddEntry(corrections[i][det], names[set_name].c_str(), "l");
  }

  // Add the full error to the stack
  for (const std::string& det : dets) {
    hfull[det]->SetLineWidth(2);
    hfull[det]->SetLineColor(kBlack);
    stacks[det].push_back(hfull[det]);
    legends[det]->AddEntry(hfull[det], "Full", "l");
  }

  // Draw the stacks
  gStyle->SetOptStat(0);
  for (const std::string& det : dets) {
    if ( !stacks[det].size() )
	continue;

    TCanvas *c = new TCanvas("c", "c", 1200, 800);
    for (TH1F* hist : stacks[det]) {
      hist->SetMinimum(0.);
      hist->SetMaximum(.04);
      hist->SetTitle(";A_{#perp}  [mm];");
      hist->Draw("HIST SAME");
    }

    legends[det]->Draw("SAME");
    c->SaveAs(TString::Format("syst_stack_%s.pdf", det.c_str()));
    delete c;
  }

  Pitch::print(Pitch::info, "Moving the systematics graphs to "+run_name+"/syst");
  std::string sysCmd = "mkdir -p "+run_name+"/syst; mv syst*.pdf *_syst.root "+run_name+"/syst; done";
  if ( std::system(sysCmd.c_str()) )
      Pitch::print(Pitch::error, "Couldn't move systematics plots");
}
