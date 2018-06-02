#include "Extractor.hh"

namespace Beam {

Extractor::Extractor() :
  _data_files(), _run_name(), _dictionary() {
}

Extractor::Extractor(const std::vector<std::string>& data_files) :
  _data_files(data_files), _run_name(), _dictionary() {

  try {
    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "Extractor::Extractor"));
  }
}

Extractor::Extractor(const Extractor& ext) {
  *this = ext;
}

Extractor& Extractor::operator=(const Extractor& ext) {
  if ( this == &ext )
      return *this;

  _data_files = ext._data_files;
  _run_name = ext._run_name;
  _dictionary = ext._dictionary;

  return *this;
}

Extractor::~Extractor () {}

const std::vector<std::string> Extractor::GetListOfKeys() const {

  std::vector<std::string> keys;
  for (const std::pair<std::string, std::string> key: _dictionary)
      keys.push_back(key.first);

  return keys;
}

Bunch Extractor::GetBunch(const std::string& data_type,
		          const size_t plane_id,
		          const bool measerr) const {
  try {
    Pitch::print(Pitch::debug, "Calling GetStream() with one plane id", "Extractor::GetBunch");
    return GetStream(data_type, {plane_id}, measerr)[plane_id];
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not get the bunch at "+std::to_string((int)plane_id)+std::string(e.what()),
	  "Extractor::GetBunch"));
  }
}

Stream Extractor::GetStream(const std::string& data_type,
	    	 	    const std::vector<size_t> plane_ids,
			    const bool measerr) const {

  try {
    Pitch::print(Pitch::debug, "Calling GetStreams() with one data type", "Extractor::GetStreams");
    std::map<std::string, std::vector<size_t>> plane_ids_map;
    plane_ids_map[data_type] = plane_ids;
    return GetStreams({data_type}, plane_ids_map, measerr).at(data_type);
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not get the stream for "+DataTypeName[data_type]+std::string(e.what()),
	  "Extractor::GetBunch"));
  }
}

std::map<std::string, Stream>
	Extractor::GetStreams(const std::vector<std::string>& data_types,
	    	 	      const std::map<std::string, std::vector<size_t>> plane_ids,
			      const bool measerr) const {

  // Sort the plane IDs for the binary search
  std::map<std::string, std::vector<size_t>> ids;
  for (const std::string type : data_types)
    if ( plane_ids.find(type) != plane_ids.end() ) {
      ids[type] = plane_ids.at(type);
      std::sort(ids[type].begin(), ids[type].end());
    }

  // Loop over the data files and extract the requested data
  std::map<std::string, std::map<size_t, BunchMap>> samples, errors;
  std::map<std::string, std::map<size_t, double>> positions;
  std::map<std::string, std::vector<size_t>> reaches;
  for (const std::string& file : _data_files) {
    // Set up the ROOT file and data pointer
    Pitch::print(Pitch::debug, "Processing "+file, "Extractor::GetStreams");
    TFile data_file(file.c_str());		// Load the MAUS output file

    // Import the truth if it is provided
    for (const std::string type : data_types) {
      if ( data_file.GetListOfKeys()->Contains(DataTypeName[type].c_str()) ) {
        Pitch::print(Pitch::debug, "Extracting the "+DataTypeName[type], "Extractor::GetStreams");
        TNtuple* tree = (TNtuple*)data_file.Get(DataTypeName[type].c_str());
        float* ntuple;
        size_t i, plane_id, particle_id, buff_id(0), n((size_t)tree->GetEntries());
        double pos_buff(-DBL_MAX);
        bool mctruth = type.find("truth") != std::string::npos;
        size_t off = mctruth ? 0 : 1;

        ProgressBar pbar(Pitch::debug);
        for (i = 0; i < n; ++i) {

          // Display the progress in %
          pbar.GetProgress(i, n);

	  // Fetch the variables from the Ntuple
	  tree->GetEntry(i);
	  ntuple = tree->GetArgs();
	  plane_id = mctruth ? ntuple[2] : ntuple[2]*5+ntuple[3];
	  if ( ids[type].size() )
	    if ( !std::binary_search(ids[type].begin(), ids[type].end(), plane_id) )
	        continue;

	  // Fill the samples
          samples[type][plane_id]["x"].push_back(ntuple[off+3]);
          samples[type][plane_id]["y"].push_back(ntuple[off+4]);
          samples[type][plane_id]["px"].push_back(ntuple[off+6]);
          samples[type][plane_id]["py"].push_back(ntuple[off+7]);
          samples[type][plane_id]["pz"].push_back(ntuple[off+8]);

	  // Fill the errors if relevant
	  if ( !mctruth && measerr ) {
            errors[type][plane_id]["x"].push_back(ntuple[10]);
            errors[type][plane_id]["y"].push_back(ntuple[11]);
            errors[type][plane_id]["px"].push_back(ntuple[13]);
            errors[type][plane_id]["py"].push_back(ntuple[14]);
            errors[type][plane_id]["pz"].push_back(ntuple[15]);
	  }

	  // Set the position of the requested planes (faster to overwrite each time)
	  positions[type][plane_id] = ntuple[off+5];

 	  // Record the reach of the particle (plane id furthest downstream)
	  particle_id = 1e6*ntuple[0]+ntuple[1];
	  if ( particle_id != buff_id )
	      pos_buff = -DBL_MAX;

          if ( positions[type][plane_id] > pos_buff ) {
	    if ( particle_id != buff_id ) {
	      reaches[type].push_back(plane_id);
	    } else {
	      reaches[type].back() = plane_id;
	    }
            pos_buff = positions[type][plane_id];
          } 

	  buff_id = particle_id;
        }
      }
    }
  }

  // If samples is empty, the requested type of data was not provided, throw
  for (const std::string& type : data_types)
    if ( !samples[type].size() )
        throw(Exceptions::Exception(Exceptions::nonRecoverable,
	      DataTypeName[type]+" is missing from the files",
	      "Extractor::GetStreams"));

  // If one of the requested plane ids is missing, throw
  for (const std::string& type : data_types)
    for (const size_t& id : ids[type])
      if ( samples[type].find(id) == samples[type].end() )
          throw(Exceptions::Exception(Exceptions::nonRecoverable,
	        "Plane "+std::to_string((int)id)+" is missing from the files for "+DataTypeName[type],
	        "Extractor::GetStreams"));

  // Initialize the stream objects
  std::map<std::string, Stream> streams;
  std::string name;
  for (const std::string& type : data_types) {
    std::map<size_t, Bunch> bunches;
    bool mctruth = type.find("truth") != std::string::npos;
    for (const auto& el : samples[type]) {
      int id = el.first;
      name = type+"_"+std::to_string((int)id);
      if ( !mctruth && measerr ) {
        bunches[id] = Bunch(samples[type][id], errors[type][id], positions[type][id], name);
      } else {
        bunches[id] = Bunch(samples[type][id], positions[type][id], name);
      }
    }

    streams[type] = Stream(bunches, reaches[type]);
  }

  return streams;
}

void Extractor::Initialize() {

  // Check that there is at list one data file in the list
  if ( !_data_files.size() )
      throw(Exceptions::Exception(Exceptions::recoverable,
	    "No data file specified, please execute in the following fashion:\n"
	    "./amplitudes [options] import0.root [... importN.root]",
	    "Extractor::Initialize"));

  // Check that the files are of the right type and exist
  for (const std::string& file : _data_files) {
    if ( std::string(file).find(".root") > std::string(file).size() )
        throw(Exceptions::Exception(Exceptions::nonRecoverable,
	      "Invalid data file type: "+file,
	      "Extractor::Initialize"));
    if ( !Exists(file) )
        throw(Exceptions::Exception(Exceptions::nonRecoverable,
	      "Data file not found at: "+file,
	      "Extractor::Initialize"));
  }

  // Get the MAUS version, check that they are consistent, issue a warning otherwise
  std::string buffer;
  for (const std::string& file : _data_files) {
    // Set up the ROOT file and data pointer
    TFile data_file(file.c_str());		// Load the MAUS output file
    if ( !data_file.IsOpen() )
        throw(Exceptions::Exception(Exceptions::nonRecoverable,
	      "Data file could not be opened: "+file,
	      "Extractor::Initialize"));

    // Loop over the list of keys, save them and their title
    TIter next(data_file.GetListOfKeys());
    TKey *key;
    while ( (key = (TKey*)next()) ) {
      // If the key already exists in the dictionary, check that the values of the
      // key match. Warn if they do not and do not overwrite.
      if ( _dictionary.find(key->GetName()) != _dictionary.end() )
	if ( key->GetTitle() != _dictionary.at(key->GetName()) ) {
	  Pitch::print(Pitch::warning, 
	      std::string(key->GetName())+" in\n\n\t"+file+" ("+std::string(key->GetTitle())+")\n\n"
	      "does not match the one in\n\n\t"+buffer+" ("+_dictionary.at(key->GetName())+")\n",
	      "Extractor::Initialize");
	  continue;
 	}

      _dictionary[key->GetName()] = key->GetTitle();
    }
    buffer = file;
  }

  // Get the run name from the first file of the list
  _run_name = std::regex_replace(_data_files[0],
	std::regex("import_|(.*/)|(.root)"), std::string(""));
  Pitch::print(Pitch::debug, "Run name: "+_run_name, "Extractor::Initialize");
}

bool Extractor::Exists(const std::string& file) const {

  // std::filesystem::exists does not exist prior to C++17, too new.
  // Rather use posix stat() that serves the same purpose but is compatible with C++11.
  struct stat buffer;   
  return (stat (file.c_str(), &buffer) == 0); 
}
} // namespace Beam
