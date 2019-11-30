// Cpp includes
#include <iostream>
#include <map>

// Additional modules
#include "InfoBox.hh"
#include "Globals.hh"
#include "Bunch.hh"
#include "Extractor.hh"
#include "Drawer.hh"

/** @file  Profiles.cc
 *
 *  @brief Produces phase space profiles.
 *
 *  	   Algorithm that produces phase space profiles of the beam in the reference plane
 *   	   upstream and downstream of the MICE absorber. It also stores the profiles at every
 *         plane to a ROOT file. 
 **/

/** @brief	Main function
 *
 *  @param	argc		Number of command line arguments
 *  @param	argv		List of command line aguments
 */
int main(int argc, char ** argv) {

  // Reroot the info messages to the terminal for now
  Pitch::setAnOutput(Pitch::info, std::cerr);

  // SynAPSE global parameters (stored in globals)
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
  std::map<std::string, Beam::Stream> streams = ext.GetStreams(types, plane_ids, globals["maxn"]);

  // Intialize the drawer and the info box
  Beam::Drawer drawer;
  if ( globals["mice"] ) {
    std::string data_type = globals["data"] ? "Preliminary" : "[simulation]";
    std::string maus_version = ext.Get("MausVersion");
    InfoBox* info = new InfoBox(data_type, maus_version,
  	globals["run_name"].AsString(), globals["user_cycle"].AsString());
    info->SetPosition("tl", .02, .24);
    drawer.SetInfoBox(info);
  }

  // Save the profiles
  std::vector<std::string> vars = {"x", "y", "px", "py", "pz"};
  TFile *outfile = new TFile("profiles_histograms.root", "RECREATE");
  outfile->cd();
  size_t ix, iy;
  for (const std::string& type : types) {
    for (const size_t& id : streams[type].GetPlaneIds()) {
      for (ix = 0; ix < vars.size(); ix++) {
	streams[type][id].Histogram(vars[ix])->Write();
        for (iy = ix+1; iy < vars.size(); iy++)
	    streams[type][id].Histogram(vars[ix], vars[iy])->Write();
      }
    }
  }
  outfile->Close();

  // Draw the profiles in the upstream and downstream reference planes
  size_t tempid;
  for (const std::string& var : vars) {
    for (size_t id : {0, 5}) {
      // Get the tracker name and station id
      std::string det = id/5 ? "tkd" : "tku";
      size_t st = id%5;

      // Initialiaze the stack
      std::string name = "profile_"+det+"_"+std::to_string((int)st)+"_"+var;

      // Get different types of data histogram
      std::map<std::string, TH1F*> hists;
      for (const std::string& type : types) {
	if ( type.find("truth") != std::string::npos ) {
	  tempid = (size_t)globals[det+"_vid"];
        } else {
	  tempid = id;
        }
  	hists[type] = streams[type][tempid].Histogram(var);
	hists[type]->Sumw2();
        hists[type]->Scale(1./hists[type]->GetEntries());
      }

      // Draw
      drawer.SaveStack(hists, name, true);
    }
  }
  outfile->Close();

  // Move the plots to the appropriate directory
  Pitch::print(Pitch::info, "Moving the beam profiles to "+run_name+"/profiles");
  std::string sysCmd = "mkdir -p "+run_name+"/profiles; "
 	"mv profile_*.pdf profiles_histograms.root "+run_name+"/profiles; done";
  if ( std::system(sysCmd.c_str()) )
      Pitch::print(Pitch::error, "Couldn't move the profile plots");
}
