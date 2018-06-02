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

void ParticleIdentification::InitializeTOF(std::vector<std::string> data_files, size_t nentries) {

  // Set the method
  Pitch::print(Pitch::info, "Fitting the time-of-flight TOF"+std::to_string(_tof_ids[0])
			    +"->"+std::to_string(_tof_ids[1])+" for PID ("
			    +std::to_string(nentries)+" requested)");

  // Counter for the number of particles found so far
  ProgressBar pbar;
  size_t n(0);

  // Container for all the time of flights
  std::vector<double> tof;
  _tof_dz = 0.;

  for (size_t iFile = 0; iFile < data_files.size(); iFile++) {
    // Set up the ROOT file and data pointer
    TFile data_file(data_files[iFile].c_str()); // Load the MAUS output file
    if ( !data_file.IsOpen() )
        throw(Exceptions::Exception(Exceptions::nonRecoverable,
	      "Data file not found: "+data_files[iFile],
	      "ParticleIdentification::InitializeTOF"));
    TTree *T = (TTree*)data_file.Get("Spill");	// Pull out the TTree
    MAUS::Data *data_ptr = new MAUS::Data(); 	// A variable to store the Data from each spill
    T->SetBranchAddress("data", &data_ptr);  	// Set the address of data_ptr

    // Loop over the spills, get the recon events
    for (size_t i = 0; i < (size_t)T->GetEntries(); ++i) {

      // Update the spill pointed to by data_ptr
      T->GetEntry(i);  
      MAUS::Spill* spill = data_ptr->GetSpill();  
      if (spill == NULL || !(spill->GetDaqEventType() == "physics_event"))
	  continue;
      std::vector<MAUS::ReconEvent*>* revts = spill->GetReconEvents();

      // Loop over recon events in spill
      for ( size_t ev = 0; ev < revts->size(); ++ev ) {
        if ( !revts->at(ev) )
	      continue;

	// Etract TOF space points
        MAUS::TOFEvent *tofevent = revts->at(ev)->GetTOFEvent();
        MAUS::TOFEventSpacePoint tofsp = tofevent->GetTOFEventSpacePoint();
        std::vector<std::vector<MAUS::TOFSpacePoint>> tofsps = {tofsp.GetTOF0SpacePointArray(),
							   	tofsp.GetTOF1SpacePointArray(),
							  	tofsp.GetTOF2SpacePointArray()};

	// Get the times and positions in each TOF, return if not a single SP
  	std::vector<double> t(2);
	bool found(true);
  	size_t i;
  	for (i = 0 ; i < _tof_ids.size(); i++) {
    	  if ( tofsps[_tof_ids[i]].size() == 1 ) {
      	    t[i] = tofsps[_tof_ids[i]][0].GetTime();
    	  } else {
      	    found = false;
    	  }
  	}
	if ( !found )
	    continue;

 	// Skip absurd values of the time-of-flight (misreconstruction)
        if ( t[1] - t[0] < 0 || t[1] - t[0] > 50 )
	    continue;

	// If the distance between TOFs is not yet known, extract it
        if ( !_tof_dz )
	    _tof_dz = tofsps[_tof_ids[1]][0].GetGlobalPosZ()
			- tofsps[_tof_ids[0]][0].GetGlobalPosZ();

        // Display the progress in %
	pbar.GetProgress(n, nentries);

	// Fill the TOF histogram and increment the count
	tof.push_back(t[1] - t[0]);
	
	if ( ++n > nentries-1 )
	    break;
      }
      if ( n > nentries-1 )
	  break;
    }
    data_file.Close();
    if ( n > nentries-1 )
	break;
  }
  if ( n < nentries )
      Pitch::print(Pitch::warning, "Requested number of particles not met ("+
		   std::to_string(n)+"/"+std::to_string(nentries)+"), will attempt to fit");

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
			TString::Format("Time-of-flight;TOF_{%d%d}",
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

int ParticleIdentification::GetID(MAUS::ReconEvent* event) const {

  if ( _method == "tof" ) {
    return GetTofID(event);
  } 
    
  throw(Exceptions::Exception(Exceptions::recoverable,
	"Particle identification method not initialized: "+_method,
	"ParticleIdentification::GetID"));
  return 0;
}

int ParticleIdentification::GetTofID(MAUS::ReconEvent* event) const {

  // Extract TOF space points
  MAUS::TOFEvent *tofevent = event->GetTOFEvent();
  MAUS::TOFEventSpacePoint tofsp = tofevent->GetTOFEventSpacePoint();
  std::vector<std::vector<MAUS::TOFSpacePoint>> tofsps = {tofsp.GetTOF0SpacePointArray(),
							  tofsp.GetTOF1SpacePointArray(),
							  tofsp.GetTOF2SpacePointArray()};

  // Get the times and positions in each TOF, return if not a single SP
  std::vector<double> t(2);
  size_t i;
  for (i = 0 ; i < _tof_ids.size(); i++) {
    if ( tofsps[_tof_ids[i]].size() == 1 ) {
      t[i] = tofsps[_tof_ids[i]][0].GetTime();
    } else {
      return 0;
    }
  }

  // Return whichever the particle hypothesis for whichever peak it is
  double dt = t[1] - t[0]; 		// [ns], time-of-flight
  if ( dt < _tof_params[0] ) {
    return 11;
  } else if ( dt < _tof_params[1] ) {
    return 13;
  }

  return 211;
}
