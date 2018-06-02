// Cpp includes
#include <iostream>
#include <vector>
#include <map>

// Additional modules
#include "Globals.hh"
#include "Extractor.hh"

/** @file  PhaseSpace.cc
 *
 *  @brief Main phase space density evolution code.
 *
 *	   Reconstructs a broad set of phase space statistics and outputs their evolution
 *	   along the beam line to ROOT TCanvases and ROOT TFiles.
 **/

/** @brief	Main function
 *
 *  @param	argc		Number of command line arguments
 *  @param	argv		List of command line aguments
 */
int main(int argc, char ** argv) {

  // Reroot the info messages to the terminal for now (TODO)
  Pitch::setAnOutput(Pitch::debug, std::cerr);

  // Emittance algorithm global parameters (stored in globals)
  Globals &globals = Globals::GetInstance(argc, argv);

  // Types and plane ids of the data to extract. Empty means import all (recmc and data)
  std::vector<std::string> types;
  for (const std::string type : {"truth", "recmc", "data"})
    if ( globals[type] )
	types.push_back(type);

  size_t min_vid((size_t)globals["min_vid"]), max_vid((size_t)globals["max_vid"]);
  std::map<std::string, std::vector<size_t>> plane_ids;
  if ( globals["truth"] )
    for (size_t i = min_vid; i <= max_vid; i++)
        plane_ids["truth"].push_back(i);

  // Get the requested streams from the beam extractor
  Beam::Extractor ext(globals.GetDataFiles());
  std::string run_name = ext.GetRunName();
  std::map<std::string, Beam::Stream> streams = ext.GetStreams(types, plane_ids);

  // If the corrected amplitudes are requested, set them
  // Set the core fraction of the streams, reconstruct the core volume
  for (const std::string& type : types) {
    streams[type].SetCoreFraction(globals["frac"]);    
    if ( globals["corrected"] )
      for (const size_t& i : streams[type].GetPlaneIds())
          streams[type][i].SetCorrectedAmplitudes();

    if ( !globals["de"] ) {
      for (const size_t& i : streams[type].GetPlaneIds())
          streams[type][i].SetCoreVolume();
    } else {
      for (const size_t& i : streams[type].GetPlaneIds())
          streams[type][i].SetCoreVolume("de", 1./streams[type].Transmission(i));
    }
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

  // List of requested beam summary statistics
  std::vector<Beam::SumStat> stats =
	{Beam::alpha, Beam::beta, Beam::gamma, Beam::mecl, Beam::neps,
	 Beam::amp, Beam::subeps, Beam::vol, Beam::mom, Beam::trans};
  std::vector<Beam::SumStat> fracstats = {Beam::amp, Beam::subeps, Beam::vol};

  // Draw the evolution graphs for a set of selected summary statistics
  TFile *outfile = new TFile(TString::Format("%s_plots.root", ext.GetRunName().c_str()), "RECREATE");
  std::map<std::string, TGraphErrors*> graphs;
  for (const Beam::SumStat& stat : stats) {
    for (const std::string& type : types) {
      graphs[type] = streams[type].EvolutionGraph(stat);
      graphs[type]->SetName(TString::Format("%s_%s", graphs[type]->GetName(), type.c_str()));
      graphs[type]->Write();
    }

    if ( !Beam::SumStatDict[stat].frac ) {
      drawer.SaveMultiGraph(graphs, Beam::SumStatDict[stat].name);
    } else {
      drawer.SaveMultiGraph(graphs, Beam::SumStatDict[stat].name,
				streams[types[0]].FractionalFunction(stat));
    }
  }
  outfile->Close();
  
  // Move all the output to an appropriate directory
  Pitch::print(Pitch::info, "Moving the phase-space evolution graphs to "+ext.GetRunName());
  std::string sysCmd = "mkdir -p "+ext.GetRunName()+"; for name in";
  for (const Beam::SumStat& stat : stats)
      sysCmd += " "+Beam::SumStatDict[stat].name;
  sysCmd += "; do mv ${name}.pdf "+ext.GetRunName()+"; done; mv *_plots.root "+ext.GetRunName();
  if ( std::system(sysCmd.c_str()) )
      Pitch::print(Pitch::error, "Couldn't move plots");

  // Draw the fractional graphs for a set of selected summary statistics
  outfile = new TFile(TString::Format("%s_frac.root", ext.GetRunName().c_str()), "RECREATE");
  std::map<Beam::SumStat, TGraphErrors*> buff;
  std::map<Beam::SumStat, std::map<std::string, TGraphErrors*>> fgraphs;
  for (const std::string& type : types) {
    buff = (type != "truth") ? streams[type].FractionalGraphs(0, 5) :
			       streams[type].FractionalGraphs(globals["tku_vid"], globals["tkd_vid"]);
    for (const Beam::SumStat& stat : fracstats) {
      fgraphs[stat][type] = buff[stat];
      buff[stat]->SetName(TString::Format("%s_%s", buff[stat]->GetName(), type.c_str()));
      buff[stat]->Write();
    }
  }
  outfile->Close();

  for (const Beam::SumStat& stat : fracstats)
      drawer.SaveMultiGraph(fgraphs[stat], Beam::SumStatDict[stat].name+"_frac");

  // Move all the output to an appropriate directory
  Pitch::print(Pitch::info, "Moving the fractional graphs to "+ext.GetRunName()+"/frac");
  sysCmd = "mkdir -p "+ext.GetRunName()+"/frac; for name in";
  for (const Beam::SumStat& stat : fracstats)
      sysCmd += " "+Beam::SumStatDict[stat].name;
  sysCmd += "; do mv ${name}_frac.pdf "+ext.GetRunName()+
	"/frac; done; mv *_frac.root "+ext.GetRunName()+"/frac";
  if ( std::system(sysCmd.c_str()) )
      Pitch::print(Pitch::error, "Couldn't move fractional plots");
}
