#include "GeometryHandler.hh"

DetectorGlobals::DetectorGlobals() :
  _pos(0., 0., 0.), _angles(0., 0., 0.), _subpos() {}

DetectorGlobals::DetectorGlobals(const std::string& filename, const std::string& name) {

  try {
    Load(filename, name);
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::recoverable,
	  "Could not initialize"+std::string(e.what()),
	  "DetectorGlobals::DetectorGlobals"));
  }
}

DetectorGlobals::DetectorGlobals(const DetectorGlobals& detector_globals) {
  *this = detector_globals;
}

DetectorGlobals& DetectorGlobals::operator=(const DetectorGlobals& detector_globals) {
  if ( this == &detector_globals )
      return *this;

  _pos = detector_globals._pos;
  _angles = detector_globals._angles;
  _subpos = detector_globals._subpos;

  return *this;
}

DetectorGlobals::~DetectorGlobals () {}

void DetectorGlobals::Load(const std::string& filename, const std::string& name) {
  std::ifstream infile(filename.c_str());
  if ( !infile.is_open() ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "The geometry file could not be found at "+filename,
	  "DetectorGlobals::Load"));    
    return;
  }

  // If a "/" character is present, there are more than one possibility
  std::vector<std::string> names;
  if ( name.find(" ") != std::string::npos ) {
    std::istringstream iss(name);
    std::copy(std::istream_iterator<std::string>(iss),
	std::istream_iterator<std::string>(),
	std::back_inserter(names));
  }

  // Try to get the position and rotations of the requested detector
  bool det_found(false);
  while ( !infile.eof() ) {
    std::string word;
    infile >> word;
    bool tracker(false), up(false);
    if ( word.find(name) != std::string::npos ||
	 (names.size() && (word.find(names[0]) != std::string::npos
		           || word.find(names[1]) != std::string::npos
		           || word.find(names[2]) != std::string::npos
		           || word.find(names[3]) != std::string::npos)) ) {
      if ( word.find("Tracker") != std::string::npos )
	  tracker = true;
      if ( word.find("Tracker0") != std::string::npos )
	  up = true;
      while ( word != "}" ) {
	infile >> word;
	if ( word == "Position" ) {
	  infile >> word;
	  _pos.SetX(std::atof(word.c_str()));
	  infile >> word;
	  _pos.SetY(std::atof(word.c_str()));
	  infile >> word;
	  _pos.SetZ(std::atof(word.c_str()));
	} else if ( word == "Rotation" ) {
	  infile >> word;
	  _angles.SetX(std::atof(word.c_str()));
	  infile >> word;
	  _angles.SetY(std::atof(word.c_str()));
	  infile >> word;
	  _angles.SetZ(std::atof(word.c_str()));
	}
      }

      // If the module is a tracker, add the position of the 5 stations
      if ( tracker ) {
	_subpos.resize(5);
	double factor = up ? -1 : 1;
	size_t i; 
	for (i = 0; i < 5; i++) {
	  _subpos[i] = _pos;
	  _subpos[i].SetZ(_pos.z()+factor*station_offsets[i]);
	}
      }

      det_found = true;
      break;
    }
  }
  infile.close();

  if ( !det_found )
      throw(Exceptions::Exception(Exceptions::recoverable,
	    name+" could not be found in "+filename,
	    "DetectorGlobals::Load"));  
}

void DetectorGlobals::SetPosition(const TVector3& pos) {

  _pos = pos;
  for (size_t i = 0; i < _subpos.size(); i++) {
    _subpos[i] = _pos;
    _subpos[i].SetZ(_pos.z()+station_offsets[i]);
  }
}

GeometryHandler::GeometryHandler() :
  _names(), _detectors() {}

GeometryHandler::GeometryHandler(const std::string& filename) :
  _names(), _detectors() {

  try {
    // Only keep the tags corresponding to the requested modules
    for (const std::pair<std::string, std::string>& names : geo_mapping) {
      _names.push_back(names.first);
      _detectors[names.first] = DetectorGlobals(filename, names.second);
    }
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::recoverable,
	  "Could not initialize"+std::string(e.what()),
	  "GeometryHandler::GeometryHandler"));
  }
}

GeometryHandler::GeometryHandler(const std::string& filename,
				 const std::vector<std::string>& names) :
  _names(names), _detectors() {

  try {
    // Only keep the tags corresponding to the requested modules
    for (const std::string& name : names)
	_detectors[name] = DetectorGlobals(filename, geo_mapping[name]);
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::recoverable,
	  "Could not initialize"+std::string(e.what()),
	  "GeometryHandler::GeometryHandler"));
  }
}

GeometryHandler::GeometryHandler(const GeometryHandler& geohandler) {
  *this = geohandler;
}

GeometryHandler& GeometryHandler::operator=(const GeometryHandler& geohandler) {
  if ( this == &geohandler )
      return *this;

  _names = geohandler._names;
  _detectors = geohandler._detectors;

  return *this;
}

GeometryHandler::~GeometryHandler () {}

const DetectorGlobals& GeometryHandler::operator[](const std::string& name) const {

  if ( _detectors.find(name) == _detectors.end() )
      throw(Exceptions::Exception(Exceptions::recoverable,
	    "Module not initialized in the GeometryHandler: "+name,
	    "GeometryHandler::operator[]"));
  return _detectors.at(name);
}

DetectorGlobals& GeometryHandler::operator[](const std::string& name) {

  if ( _detectors.find(name) == _detectors.end() )
      throw(Exceptions::Exception(Exceptions::recoverable,
	    "Module not initialized in the GeometryHandler: "+name,
	    "GeometryHandler::operator[]"));
  return _detectors[name];
}
