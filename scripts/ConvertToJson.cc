// Cpp includes
#include <iostream>
#include <vector>
#include <map>

// Root includes
#include "TSystem.h"
#include "TNtuple.h"
#include "TFile.h"

// Additional modules
#include "ProgressBar.hh"
#include "Globals.hh"
#include "Pitch.hh"

/** @file  ConvertToJson.cc
 *
 *  @brief Converts MICE phase space data to text files
 *
 *	   Algorithm that converts the imported MICE data into a set of 10 text files (one
 *	   for each tracker station) that contain a list of of lines that each corresponds
 *	   to the phase space of a single selected particle: [x, px, y, py, pz]
 **/

/** @brief	Main function
 *
 *  @param	argc		Number of command line arguments
 *  @param	argv		List of command line aguments
 */
int main(int argc, char ** argv) {

  // Reroot the info messages to the terminal for now
  Pitch::setAnOutput(Pitch::info, std::cerr);

  // Converter algorithm global parameters (stored in globals)
  Globals &globals = Globals::GetInstance(argc, argv);
  if ( !globals.GetDataFiles().size() ) {
    Pitch::print(Pitch::error, "No data file specified, please execute in the following fashion:");
    Pitch::print(Pitch::error, "./convert_to_json [options] import0.root [... importN.root]");
    return 2;
  } else {
    for (const std::string& file : globals.GetDataFiles())
      if ( std::string(file).find(".root") > std::string(file).size() ) {
        Pitch::print(Pitch::error, "Invalid data file type: "+std::string(file));
	return 2;
      }
  }

  // List of detectors
  std::vector<std::string> dets = {"tku", "tkd"};	// Detectors

  // Load libTree
  gSystem->Load("libTree");

  // Output files
  std::string out_dir = "beams/"+(std::string)globals["run_number"];
  std::string sysCmd = "mkdir -p "+out_dir;
  if ( std::system(sysCmd.c_str()) ) {
    Pitch::print(Pitch::error, "Couldn't create a directory for the beam files");
    return 2;
  }

  std::map<std::string, std::vector<std::ofstream*>> out;
  for (const std::string& det : dets) {
    out[det].resize(5);
    for (int st = 0; st < 5; st++)
        out[det][st] = new std::ofstream(out_dir+"/"+det+"_"+std::to_string(st+1)+".json");
  }

  // Fill the sample maps from the imported ROOT file
  std::string maus_version;
  std::map<std::string, double> vpoint;
  for (const std::string& file : globals.GetDataFiles()) {
    // Set up the ROOT file and data pointer
    TFile data_file(file.c_str());		// Load the MAUS output file
    if ( !data_file.IsOpen() ) {
      Pitch::print(Pitch::error, "Data file not found: "+file);
      return 2;
    }
    Pitch::print(Pitch::info, "Processing "+file);

    // Get the MAUS version if the file was produced in MAUS
    if ( data_file.GetListOfKeys()->Contains("MausVersion") ) {
      maus_version = data_file.Get("MausVersion")->GetTitle();
      maus_version = maus_version.substr(maus_version.size()-5);
      Pitch::print(Pitch::info, "MAUS version: "+maus_version);
    }

    // Import the reconstructed data if it is provided
    for (const std::string& type : {"Data"}) {
      std::string type_lc = type;
      std::transform(type_lc.begin(), type_lc.end(), type_lc.begin(), ::tolower);
      if ( data_file.GetListOfKeys()->Contains(type.c_str()) ) {
        TNtuple* rec_samples = (TNtuple*)data_file.Get(type.c_str());
        float* ntuple;
	std::string det;
        size_t st;

        Pitch::print(Pitch::info, "Processing the "+type_lc);
        ProgressBar pbar;
        for (size_t i = 0; i < (size_t)rec_samples->GetEntries(); ++i) {

	  // Fetch the variables from the Ntuple
	  rec_samples->GetEntry(i);
	  ntuple = rec_samples->GetArgs();
	  det = ntuple[2] ? "tkd" : "tku";
	  st = ntuple[3];

	  vpoint["x"] = ntuple[4];
	  vpoint["y"] = ntuple[5];
	  vpoint["px"] = ntuple[7];
	  vpoint["py"] = ntuple[8];
	  vpoint["pz"] = ntuple[9] += .6; // Compensate for eloss in the station

	  // Fill the true samples
	  *out[det][st] << "[" << vpoint["x"];
	  for (const std::string var : {"px", "y", "py", "pz"})
	      *out[det][st] << ", " << vpoint[var];
	  *out[det][st] << "]\n";

          // Display the progress in %
          pbar.GetProgress(i, (size_t)rec_samples->GetEntries());
        }
      }
    }

    data_file.Close();
  } // End of the list of files

  // Close the output files
  for (const std::string& det : dets)
    for (int st = 0; st < 5; st++) {
      out[det][st]->close();
      delete out[det][st];
    }
}
