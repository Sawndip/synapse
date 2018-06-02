// Cpp includes
#include <iostream>
#include <vector>
#include <map>

// Root includes
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"

// Additional modules
#include "ProgressBar.hh"
#include "Pitch.hh"
#include "Globals.hh"

/** @file ImportG4BL.cc
 *
 *  @brief Imports G4BL ASCII output.
 *
 *	   Algorithm that imports the G4BeamLine simulation ASCII output to a single TNtuple:
 *	    - Truth: true phase space information of each track at each virtual plane.
 **/

/** @brief	Main function
 *
 *  @param	argc		Number of command line arguments
 *  @param	argv		List of command line aguments
 */
int main(int argc, char ** argv) {

  // Importer algorithm global parameters (stored in globals)
  Globals &globals = Globals::GetInstance(argc, argv);
  if ( !globals.GetDataFiles().size() ) {
    Pitch::print(Pitch::error, "No data file specified, please execute in the following fashion:");
    Pitch::print(Pitch::error, "./import_g4bl [options] data0.txt [... dataN.txt]");
    return 2;
  } else {
    for (const std::string& file : globals.GetDataFiles())
      if ( std::string(file).find(".txt") > std::string(file).size() ) {
        Pitch::print(Pitch::error, "Invalid data file type: "+std::string(file));
	return 2;
      }
  }

  // Loop over the input data files, extract the relevant to the beam evolution (x, y, px, py, pz)
  std::vector<double> virtuals_z;
  std::string out_file = globals["import_filename"];
  if ( !out_file.size() )
      out_file = "import.txt";
  TFile *out = new TFile(out_file.c_str(), "RECREATE");
  TNtuple* truth_samples = new TNtuple("truth_ntuple", "", 
	"SpillID:EventID:VirtualPlaneID:x:y:z:px:py:pz");

  for (const std::string& file : globals.GetDataFiles()) {
    // Set up the input file
    std::ifstream data_file(file);		// Load the G4BL output file
    if ( !data_file.is_open() ) {
      Pitch::print(Pitch::error, "Data file not found: "+file);
      return 2;
    }
    Pitch::print(Pitch::info, "Processing "+file);

    // Extract the variable names from the header line
    std::string line;
    std::vector<std::string> vars;
    while ( !data_file.eof() ) {
      std::getline(data_file, line);
      if ( line.find("#x") != line.npos ) {
	line = line.substr(1);
        std::istringstream f(line);
    	std::string s;
    	while (std::getline(f, s, ' '))
            vars.push_back(s);

        std::replace_if(line.begin(), line.end(), [] (const char& c) { 
		return std::isspace(c) ;
	}, ':');
	break;
      }
    }

    // Fill an Ntuple with all the variables available
    TNtuple* tmp_samples = new TNtuple("tmp_ntuple", "", line.c_str());
    tmp_samples->ReadStream(data_file);

    // Get the list of IDs, in order, of the variables to extract
    std::vector<std::string> target =
	{"EventID", "TrackID", "VirtualPlaneID", "x", "y", "z", "Px", "Py", "Pz"};
    std::vector<size_t> ids(9);
    size_t i, j;
    for (i = 0; i < vars.size(); i++)
      for (j = 0; j < target.size(); j++)
	if ( vars[i] == target[j] ) {
	  ids[j] = i;
	  break;
	}

    // Loop and fill and the ones that need to be kept only
    size_t vid;
    float* ntuple, outntuple[9];
    bool matched;
    for (i = 0; i < (size_t)tmp_samples->GetEntries(); i++) {
      tmp_samples->GetEntry(i);
      ntuple = tmp_samples->GetArgs();

      // From the z position, figure out if it is a new vplane or not
      if ( !virtuals_z.size() ) {
	virtuals_z.push_back(ntuple[ids[5]]);
	vid = 0;
      } else {
	matched = false;
	for (j = 0; j < virtuals_z.size(); j++) {
	  if ( fabs(virtuals_z[j]-ntuple[ids[5]]) < 1 ) {	// 1 mm tolerance
	    vid = j;
	    matched = true;
	    break;
	  }
	}

	if ( !matched ) {
	  virtuals_z.push_back(ntuple[ids[5]]);
	  vid = virtuals_z.size()-1;
	}
      }

      // Fill the truth samples
      for (j = 0; j < target.size(); j++) {
	if ( j == 2 ) {
	  outntuple[j] = vid;
	  continue;
	}
	outntuple[j] = ntuple[ids[j]];
      }

      truth_samples->Fill(outntuple);
    }

    delete tmp_samples;
    data_file.close();
  }

  // Fill the output file with the Ntuples and information about the sample
  out->cd();
  if ( truth_samples->GetEntries() )
      truth_samples->Write("Truth");
  delete truth_samples;
}
