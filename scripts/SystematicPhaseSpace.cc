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

  // Set up the plane IDs mapping and other containers
  // Only need to correct the truth in the downstream tracker
  // Need to correct the reconstruction in every tracker station
  std::vector<size_t> vids = {64, 60, 57, 53, 49, 80, 84, 87, 91, 94};
  std::map<size_t, size_t> map_vids =
	{{0, 64}, {1, 60}, {2, 57}, {3, 53}, {4, 49}, {5, 80}, {6, 84}, {7, 87}, {8, 91}, {9, 94}};
  std::vector<std::string> types = {"utruth", "truth", "recmc"};
  std::vector<std::string> dets;
  size_t min_vid((size_t)globals["min_vid"]),
		tkd_vid((size_t)globals["tkd_vid"]), max_vid((size_t)globals["max_vid"]);
  std::map<std::string, std::vector<size_t>> plane_ids;
  std::vector<Beam::SumStat> stats = globals["de"] ?
		std::vector<Beam::SumStat>({Beam::den, Beam::vol}) :
		std::vector<Beam::SumStat>({Beam::amp, Beam::subeps, Beam::vol});

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

  // Loop over the simulation files, evaluate and store the corrections for summary statistic
  TFile *outfile = globals["de"] ?
	new TFile(TString::Format("%s_ps_de_syst.root", run_name.c_str()), "RECREATE") :
	new TFile(TString::Format("%s_ps_syst.root", run_name.c_str()), "RECREATE");
  std::map<Beam::SumStat, std::map<std::string, std::map<size_t, double>>> vals;
  std::map<Beam::SumStat, std::map<std::string, std::map<size_t, double>>> corrs;
  std::map<Beam::SumStat, std::map<std::string, std::map<std::string, TH1F*>>> hcorrs;
  std::string set_det;
  for (size_t iFile = 0; iFile < nsets; iFile++) {
    // Set up the extractor
    Beam::Extractor ext({globals.GetDataFiles()[iFile]});

    // Get the setting, determine which stations are required
    set_name = set_names[iFile];
    set_det = set_name.substr(0, 3);
    set_name = set_name.substr(4);
    Pitch::print(Pitch::info, "Importing the "+set_name+" systematics set for "+set_det);

    plane_ids.clear();
    if ( set_name != "base" ) {

      // Add the first simulation plane as the transmission reference
      for (const std::string type : {"truth", "utruth"})
	plane_ids[type].push_back(min_vid); // Reference for transmission

      // Add the upstream tracker planes if tku is varried
      if ( set_det == "tku" )
        for (size_t i = 0; i < 5; i++) {
	  plane_ids["recmc"].push_back(i);
	  plane_ids["utruth"].push_back(map_vids[i]);
        }

      // Add the downstream tracker planes if tkd is varried
      if ( set_det == "tkd" ) {
	plane_ids["recmc"].push_back(0); // Reference for transmission
        for (size_t i = 5; i < 10; i++) {
	  plane_ids["recmc"].push_back(i);
	  plane_ids["truth"].push_back(map_vids[i]);
	  plane_ids["utruth"].push_back(map_vids[i]);
        }
     }

      dets = {set_det};
    } else {
      // Add uncut truth at the upstream tracker planes
      for (size_t i = 0; i < 5; i++) 
  	  plane_ids["utruth"].push_back(map_vids[i]);

      // Add the truth and uncut truth for all the planes downstream of the tkd ref. plane
      for (const std::string type : {"truth", "utruth"}) {
	plane_ids[type].push_back(min_vid); // Reference for transmission
        for (size_t i = tkd_vid; i <= max_vid; i++)
	    plane_ids[type].push_back(i);
      }

      // Add all the stations of both trackers
      for (size_t i = 0; i < 10; i++)
	  plane_ids["recmc"].push_back(i);

      set_det = "";
      dets = {"tku", "tkd"};
    }

    // Extract the beam ensembles
    std::map<std::string, Beam::Stream> streams = ext.GetStreams(types, plane_ids, globals["maxn"]);

    // Set the corrected amplitudes if requested, set the core fraction
    for (const std::string& type : types) {
      if ( globals["de"] ) {
        // Set the levels
        for (const size_t& i : streams[type].GetPlaneIds())
            streams[type][i].SetDensityLevels(streams[type].Transmission(i));

        // Set the core fraction of the beam to de
        streams[type].SetCoreFraction(globals["frac"], "de");

        // Set the volumes
        for (const size_t& i : streams[type].GetPlaneIds()) {
            streams[type][i].SetCoreVolume("de", streams[type].Transmission(i));
	    std::cerr << streams[type][i].FracEmittance() << std::endl;
        }

        // Proceed
        continue;
      }

      if ( globals["corrected"] )
        for (const size_t& i : streams[type].GetPlaneIds())
            streams[type][i].SetCorrectedAmplitudes();

      // Set the core fraction of the beam to amplitudes
      streams[type].SetCoreFraction(globals["frac"]);   

      // Set the volumes
      for (const size_t& i : streams[type].GetPlaneIds())
          streams[type][i].SetCoreVolume();
    }

    // Compute the requested summary statistics in each of the planes
    for (const Beam::SumStat& stat : stats)
      for (const std::string& type : types)
        for (const size_t& i : streams[type].GetPlaneIds())
	    vals[stat][type][i] = streams[type][i].SummaryStatistic(stat).GetValue();

    // Compute the correction factors
    for (const Beam::SumStat& stat : stats) {
      for (const size_t& i : streams["recmc"].GetPlaneIds()) {
          corrs[stat]["recmc"][i] = vals[stat]["utruth"][map_vids[i]]/vals[stat]["recmc"][i];
      }

      for (const size_t& i : streams["truth"].GetPlaneIds())
          corrs[stat]["truth"][i] = vals[stat]["utruth"][i]/vals[stat]["truth"][i];
    }

    // Store the correction factors into histograms to be fetched later
    for (const Beam::SumStat& stat : stats)
      for (const std::string& type : {"truth", "recmc"}) {
	std::string name = "syst_corr_"+Beam::SumStatDict[stat].name+"_"+type+"_"+set_names[iFile];
	hcorrs[stat][type][set_name] = new TH1F(name.c_str(), "", 100, 0, 100);
        for (size_t i = 0; i < (size_t)hcorrs[stat][type][set_name]->GetNbinsX(); i++)
	    hcorrs[stat][type][set_name]->SetBinContent(i+1, 1.);

        for (const size_t& i : streams[type].GetPlaneIds()) {
	  size_t binid = (type == "recmc") ? i : i - min_vid;
	  hcorrs[stat][type][set_name]->SetBinContent(binid+1, corrs[stat][type][i]);
	}
	outfile->cd();
	hcorrs[stat][type][set_name]->Write(hcorrs[stat][type][set_name]->GetName());
      }
  } // End of the list of files

  outfile->Close();

  Pitch::print(Pitch::info, "Moving the systematics graphs to "+run_name+"/syst");
  std::string sysCmd = "mkdir -p "+run_name+"/syst; mv syst*.pdf *_syst.root "+run_name+"/syst; done";
  if ( std::system(sysCmd.c_str()) )
      Pitch::print(Pitch::error, "Couldn't move systematics plots");
}
