#include "Globals.hh"

Globals& Globals::GetInstance(int argc,
	  		      char** argv) {

  static Globals instance(argc, argv);
  return instance;
}

Globals::Globals(int argc,
		 char** argv) :
  _globals(), _data_files(), _geoh() {

  try {
    LoadFile("../ConfigurationDefaults.txt");
    LoadCommandLine(argc, argv);
    LoadGeometry();
    SetDefaultStyle();

    Pitch::print(Pitch::debug, "", "Globals::Initialize");
    Print(Pitch::debug);
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "Globals::Globals"));
  }
}

Globals::~Globals() {}

const Global& Globals::operator()(const std::string& var) const {

  // Throw if the variable has not been loaded
  if ( _globals.find(var) == _globals.end() )
      throw(Exceptions::Exception(Exceptions::recoverable,
	    "This entry has not been loaded in the Globals dictionary: "+var,
	    "Globals::operator[]"));

  return _globals.at(var);
}

const Global& Globals::operator[](const std::string& var) const {

  // Throw if the variable has not been loaded
  if ( _globals.find(var) == _globals.end() )
      throw(Exceptions::Exception(Exceptions::recoverable,
	    "This entry has not been loaded in the Globals dictionary: "+var,
	    "Globals::operator[]"));

  return _globals.at(var);
}

void Globals::Print(Pitch::errorLevel level) const {

  std::stringstream dictionary;
  dictionary << "Parameters loaded in the globals dictionary:\n";
  for (const std::pair<std::string, Global>& el : _globals)
      dictionary << std::left << std::setw(20) << el.first << el.second << "\n";        

  dictionary << "\nData files:\n";
  for (const std::string& file : _data_files)
      dictionary << file << "\n";        

  Pitch::mout(level) << dictionary.str() << "\n";
}

void Globals::LoadFile(std::string filename) {

  // Load the file in a file stream, check that it exists
  std::ifstream datacard(filename);
  if ( !datacard.is_open() ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Datacards not found, please verify path: "+std::string(filename),
	  "Globals::LoadFile"));
    return;
  } else {
    std::string line, variable, value, des;
    while ( !datacard.eof() ) {
      // Get the next line in the file
      std::getline(datacard, line);

      // Check that the line is not empty
      if ( !line.size() )
	  continue;

      // Check that the line is not a comment
      if ( line[0] == '#' || line[0] == '/' )
	  continue;

      // Get the variable name and its value from the line
      std::istringstream iss(line);
      iss >> variable >> value;
      _globals[variable] = Global(value);

      // Get its description from the line
      des = std::regex_replace(iss.str(), std::regex("(.*#)"), std::string(""));
      _globals[variable].SetDescription(des);
    }
  }
}

void Globals::LoadCommandLine(int argc, char** argv) {

  // If no argument are provided, simply return
  if ( !argc )
      return;

  // Build the option index to choose from
  std::vector<std::string> opt_names;
  for (const std::pair<std::string, Global>& el : _globals)
      opt_names.push_back(el.first);

  static std::vector<option> long_options;
  size_t optid(0);
  for (const std::pair<std::string, Global>& el : _globals) {
    long_options.push_back({opt_names[optid].c_str(), required_argument, 0, 0});
    optid++;
  }
  long_options.push_back({"configuration_file", required_argument, 0, 	 1});
  long_options.push_back({"help", 		no_argument, 	   0, 	 'h'});
  long_options.push_back({"verbose", 		no_argument, 	   0, 	 'v'});
  long_options.push_back({NULL, 		0, 		   NULL, 0});

  // Define the help message
  std::stringstream help_message;
  help_message << "Usage: ./program [options] files...\n"
		   "Options:\n";
  for (size_t i = 0; i < opt_names.size(); i++)
      help_message << "  --" << std::left << std::setw(20) << opt_names[i] 
		   <<  _globals[opt_names[i]].Description() << "\n";        

  // Look for recognized command line arguments. If the argument does not match any 
  // current entry in the dictionary, throw.
  int c;
  while ( true ) {
    int option_index = 0;

    // Get the current argument
    c = getopt_long(argc, argv, "hv", &long_options[0], &option_index);
    if ( c == -1 )
	break;

    switch ( c ) {
      // In this case, the variable was found in the cards, override
      case 0:
	_globals[long_options[option_index].name] = Global(optarg);
        break;

      // In this case, the datacards are overwritten by a custom data fie
      case 1:
	LoadFile(optarg);
        break;

      // In this case, print the help information
      case 'h':
	std::cerr << help_message.str() << std::endl;
	exit(EXIT_SUCCESS);
        break;

      // In this case, reroot the debug messages
      case 'v':
  	Pitch::setAnOutput(Pitch::info, std::cerr);
  	Pitch::setAnOutput(Pitch::debug, std::cout);
        break;

      // In this case, the variable is unknown, throw and print the list of options
      default:
        throw(Exceptions::Exception(Exceptions::recoverable,
	      "Global variable not defined in the default data cards.\n"
	      "Run ./program -h or ./program --help for more information.",
	      "Globals::Update"));	
    }
  }

  // The remaining arguments are the name of the runs to be processed
  while ( optind < argc )
      _data_files.push_back(argv[optind++]);
}

void Globals::LoadGeometry() {

  try {
    // Check that the geometry file is specified in the data cards
    if ( _globals.find("geometry_filename") == _globals.end() )
        throw(Exceptions::Exception(Exceptions::recoverable,
	      "The geoemtry filename is not specified in the data cards",
	      "Globals::LoadGeometry"));

    // Load the geometry file specified in the data cards
    _geoh = GeometryHandler(_globals["geometry_filename"]);

  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not load geometry"+std::string(e.what()),
	  "Globals::LoadGeometry"));
  }
}

void Globals::SetDefaultStyle() {

  // Load the libTree, prevents unwanted SegFaults
  gSystem->Load("libTree");

  // Set the default palette to kBird (ROOT 6 default, B&W proof);
  std::vector<double> stops = {0., 0.125, .25, .375, .5, .625, .75, .875, 1.};
  std::vector<double> red = {0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764};
  std::vector<double> green = {0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832};
  std::vector<double> blue = {0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539};
  TColor::CreateGradientColorTable(9, &stops[0], &red[0], &green[0], &blue[0], 255, 1);

  // Set the limits of the stat box and its default options
  gStyle->SetOptStat(2210);
  gStyle->SetStatX(.9);
  gStyle->SetStatY(.9);

  // Set the default label and title size to slightly larger than ROOT default
  gStyle->SetTitleSize(.04, "XYZ");
  gStyle->SetLabelSize(.04, "XYZ");
  gStyle->SetTitleOffset(1.25, "Y");
}
