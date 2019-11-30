#include "ParticleIdentification.hh"

ParticleIdentification::ParticleIdentification() :
  _method(""), _tof_ids(0), _tof_params(0), _tof_dz(0.) {
}

ParticleIdentification::ParticleIdentification(const std::vector<std::string>& data_files,
			 		       const std::vector<size_t>& ids,
			 		       const size_t nentries) :
  _method("tof"), _tof_ids(ids), _tof_params(0), _tof_dz(0.) {

  
  try {
    InitializeTOF(data_files, nentries);
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize the TOF method"+std::string(e.what()),
	  "ParticleIdentification::ParticleIdentification"));
  }
}

ParticleIdentification::ParticleIdentification(const std::vector<std::string>& data_files,
					       const std::vector<size_t>& ids,
			 		       const double& mumin,
			 		       const double& mumax,
			 		       const size_t nentries) :
  _method("tof"), _tof_ids(ids), _tof_params(0), _tof_dz(0.) {

  
  try {
    InitializeTOFPeak(data_files, mumin, mumax, nentries);
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize the TOF method"+std::string(e.what()),
	  "ParticleIdentification::ParticleIdentification"));
  }
}

ParticleIdentification::ParticleIdentification(const std::vector<size_t>& ids,
			 		       const double& mumin,
			 		       const double& mumax) :
  _method("tof"), _tof_ids(ids), _tof_params(0), _tof_dz(0.) {

  
  try {
    _tof_params.push_back(mumin);
    _tof_params.push_back(mumax);
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize the TOF method"+std::string(e.what()),
	  "ParticleIdentification::ParticleIdentification"));
  }
}

ParticleIdentification::ParticleIdentification(const ParticleIdentification& pid) {
  *this = pid;
}

ParticleIdentification& ParticleIdentification::operator=(const ParticleIdentification& pid) {
  if ( this == &pid )
      return *this;

  _method = pid._method;
  _tof_ids = pid._tof_ids;
  _tof_params = pid._tof_params;
  _tof_dz = pid._tof_dz;

  return *this;
}

ParticleIdentification::~ParticleIdentification () {}

void ParticleIdentification::InitializeTOF(std::vector<std::string> data_files,
					   size_t nentries) {

  // Get the distance between the TOFs from the geometry
  Globals &globals = Globals::GetInstance();
  const GeometryHandler& geoh = globals.GetGeometryHandler();
  _tof_dz = geoh["tof"+std::to_string((int)_tof_ids[1])].z()-
	    geoh["tof"+std::to_string((int)_tof_ids[0])].z();

  // Get the array of TOFs
  std::vector<double> tof = TOFs(data_files, 0, 50, nentries);

  // Find the time taken by a particle travelling at the speed of light
  double c = 299792458;			// [m/s]
  double ct = 1e9*_tof_dz*1e-3/c; 	// [ns]

  // Define the limits of the TOF histogram from the recorded data
  double mean = Math::Mean(tof);
  double sig = sqrt(Math::Covariance(tof, tof));
  double min = std::max(mean-3*sig, ct+.5);	// A bit slower than light or far enough from mus
  double max = mean+3*sig;			// Far enough from the mean to encompass everything
  size_t nbins = (max-min)/.12;			// 120 ns bins (2*resolution)

  // Define a histogram that will contain the time of flight profile
  // Bin width ~ 250 ps	(4*res)
  TH1F* htof = new TH1F(TString::Format("tof%d%d", (int)_tof_ids[0], (int)_tof_ids[1]),
			TString::Format(";TOF_{%d%d}",
			(int)_tof_ids[0], (int)_tof_ids[1]), nbins, min, max);
  htof->FillN(tof.size(), &(tof[0]), NULL);

  // Find the peaks using the root TSpectrum code
  const size_t npeaks = 2; // muons, pions
  TSpectrum *spc = new TSpectrum(npeaks);
  size_t nfound = spc->Search(htof);

  // Fit their widths if the right amount is found
  if ( npeaks == nfound ) {
    float *xpeaks = spc->GetPositionX(); // array of npeaks floats
    std::vector<double> posx;
    for (size_t i = 0; i < npeaks; i++)
        posx.push_back((double)xpeaks[i]);
    std::sort(posx.begin(), posx.end());

    // Define a two peaks gaussian to fit to the TOF distribution
    TF1 *fpeaks = new TF1("fpeaks",
			  "[0]*TMath::Gaus(x, [1], [2])+[3]*TMath::Gaus(x, [4], [5])",
			  25, 45);
    fpeaks->SetParNames("A_{#mu}","#mu_{#mu}","#sigma_{#mu}",
			"A_{#pi}","#mu_{#pi}","#sigma_{#pi}");

    // Set first approximation of parameters
    for(size_t p = 0; p < npeaks; p++) {
      int bin = htof->GetXaxis()->FindBin(posx[p]);
      double yp = htof->GetBinContent(bin);

      fpeaks->FixParameter(3*p, yp);
      fpeaks->FixParameter(3*p+1, posx[p]);
      fpeaks->SetParameter(3*p+2, .1);
    }

    // Fit peaks with the two peak gaussian and extract the boundaries between peaks
    htof->Fit("fpeaks", "RQ");

    double tmu(fpeaks->GetParameter(1)), tpi(fpeaks->GetParameter(4)),
		sigmu(fpeaks->GetParameter(2)), sigpi(fpeaks->GetParameter(5));
    _tof_params.resize(2);
    _tof_params[0] = ct+.5;					// Limit between e^+ and mu^+
    _tof_params[1] = (tmu*sigpi+tpi*sigmu)/(sigmu+sigpi);	// Limit between mu^+ and pi^+
  }

  // Plot the histogram
  gStyle->SetStatX(.9);
  gStyle->SetStatY(.9);
  gStyle->SetStatH(.2);
  gStyle->SetStatW(.3);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  TCanvas *canv = new TCanvas("c", "c", 1200, 800);
  htof->Draw();
  canv->SaveAs(TString::Format("tof%d%d.pdf", (int)_tof_ids[0], (int)_tof_ids[1]));
  delete canv;
  delete htof;

  // If the wrong amount of peak was found, throw
  if ( npeaks != nfound )
      throw(Exceptions::Exception(Exceptions::recoverable,
	    "Wrong number of peaks found ("+std::to_string(nfound)+"/2)",
	    "ParticleIdentification::InitializeTOF"));
}

void ParticleIdentification::InitializeTOFPeak(std::vector<std::string> data_files,
					       const double& mumin,
			   		       const double& mumax,
			 		       const size_t nentries) {

  // Get the distance between the TOFs from the geometry
  Globals &globals = Globals::GetInstance();
  const GeometryHandler& geoh = globals.GetGeometryHandler();
  _tof_dz = geoh["tof"+std::to_string((int)_tof_ids[1])].z()-
	    geoh["tof"+std::to_string((int)_tof_ids[0])].z();

  // Find the time theoretical time taken by a particle travelling at the speed of light
  double c = 299792458;			// [m/s]
  double ct = 1e9*_tof_dz*1e-3/c; 	// [ns]

  // Find the position of the electron peak
  std::vector<double> tof = TOFs(data_files, ct-1, ct+1, nentries);
  double peak = Math::Mean(tof);

  // Set the distance from measured electron peak
  _tof_dz = c*peak*1e-9;

  // Draw histogram of the peak
  TH1F* htof = new TH1F(TString::Format("tof%d%d", (int)_tof_ids[0], (int)_tof_ids[1]),
			TString::Format(";TOF_{%d%d}",
			(int)_tof_ids[0], (int)_tof_ids[1]), 20, ct-1, ct+1);
  htof->FillN(tof.size(), &(tof[0]), NULL);

  TCanvas *canv = new TCanvas("c", "c", 1200, 800);
  htof->Draw();
  canv->SaveAs(TString::Format("tof%d%d_epeak.pdf", (int)_tof_ids[0], (int)_tof_ids[1]));
  delete canv;
  delete htof;

  // Set the limits
  _tof_params.push_back(peak+mumin);
  _tof_params.push_back(peak+mumax);
}

int ParticleIdentification::GetID(const double& fom) const {

  if ( _method == "tof" ) {
    return GetTofID(fom);
  } 
    
  throw(Exceptions::Exception(Exceptions::recoverable,
	"Particle identification method not initialized: "+_method,
	"ParticleIdentification::GetID"));
  return 0;
}

int ParticleIdentification::GetTofID(const double& tof) const {

  // Return whichever the particle hypothesis for whichever peak it is
  if ( tof < _tof_params[0] ) {
    return 11;
  } else if ( tof < _tof_params[1] ) {
    return 13;
  }

  return 211;
}

std::vector<double> ParticleIdentification::TOFs(const std::vector<std::string>& data_files,
			   			 const double& min,
			   			 const double& max,
			   			 const size_t nentries) {

  // Set the method
  Pitch::print(Pitch::info, "Extracting the time-of-flight TOF"+std::to_string(_tof_ids[0])
			    +"->"+std::to_string(_tof_ids[1])+" for PID ("
			    +std::to_string(nentries)+" requested)");

  // Counter for the number of particles found so far
  ProgressBar pbar;
  size_t n(0);

  // Container for all the time of flights
  std::vector<double> tof;
  double t;

  for (const std::string& file : data_files) {
    // Set up the ROOT file and data pointer
    TFile data_file(file.c_str()); // Load the MAUS output file
    if ( !data_file.IsOpen() )
        throw(Exceptions::Exception(Exceptions::nonRecoverable,
	      "Data file not found: "+file,
	      "ParticleIdentification::TOFs"));

    TTree *T = (TTree*)data_file.Get("Track");	// Pull out the TTree
    MiceTrack *track = NULL;	 		// A variable to store the data from each track
    T->SetBranchAddress("Recon", &track);  	// Set the address of data_ptr

    for (size_t i = 0; i < (size_t)T->GetEntries(); ++i) {

      // Update the track pointed to by *track	
      T->GetEntry(i);

      // Skip if there is not exactly one SP in each TOF station
      if ( track->tof_nsp[_tof_ids[0]] != 1 || track->tof_nsp[_tof_ids[1]] != 1 )
	  continue;
      t = track->t[_tof_ids[1]]-track->t[_tof_ids[0]];

      // Skip times that do not fall inside the requested interval
      if ( t < min || t > max )
	  continue;

      // Display the progress in %
      pbar.GetProgress(n, nentries);

      // Fill the TOF vector and increment the count
      tof.push_back(t);
	
      if ( ++n > nentries-1 )
	  break;
    }

    data_file.Close();
    if ( n > nentries-1 )
	break;
  }

  // Warn if the amount of particles requested is not met
  if ( n < nentries ) {
    pbar.GetProgress(nentries-1, nentries);
    Pitch::print(Pitch::warning, "Requested number of particles requested not met ("+
		 std::to_string(n)+"/"+std::to_string(nentries)+")",
		 "ParticleIdentification::TOFs");
  }

  return tof;
}
