#include "Scattering.hh"

Scattering::Scattering() :
  _rdmzer(time(NULL)), _thickness(0.), _factor(0.) {
}

Scattering::Scattering(const double& thickness) :
  _rdmzer(time(NULL)), _thickness(thickness), _factor(0.) {

  ComputeFactor(_thickness);
}

Scattering::Scattering(const Scattering& scat) {
  *this = scat;
}

Scattering& Scattering::operator=(const Scattering& scat) {
  if ( this == &scat )
      return *this;

  _rdmzer = scat._rdmzer;
  _thickness = scat._thickness;
  _factor = scat._factor;

  return *this;
}

Scattering::~Scattering () {}

double Scattering::Scatter(const double& p,
			   const double& m) {

  // Figure out the velocity from p and m
  double beta = (p/m)/sqrt(1.+pow(p/m, 2));

  // Get the RMS scattering angle theta0
  double theta0 = _factor/(p*beta);

  // Increment the seed and get a random Gaussian scatter
  return _rdmzer.Gaus(0, theta0);
}

void Scattering::ScatterBunch(std::map<std::string, std::vector<double>>& beam,
			     const double& p,
			     const double& m) {

  // Check weather it is a 2 or 4D beam to scatter
  std::vector<std::string> moms;
  for (const std::string mom : {"px", "py"})
    if ( beam.find(mom) != beam.end() )
        moms.push_back(mom);
  if ( !moms.size() || beam.find("pz") == beam.end() )
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	    "Not enough components of the momentum provided",
	    "Scattering::ScatterBunch"));

  // Scatter the particles in the beam
  double angi, fact;
  std::map<std::string, double> ango;
  size_t i;
  for (i = 0; i < beam["pz"].size(); i++) {
    // First get the angles in all the requested transverse dimensions
    for (const std::string& par : moms) {
      angi = beam[par][i]/beam["pz"][i];
      ango[par] = angi + Scatter(p, m);
      fact += pow(ango[par], 2);
    }

    // First compute the output longitudinal momentum from the angles
//    beam["pz"][i] /= sqrt(fact);

    // Then the output transverse components
    for (const std::string& par : moms)
	beam[par][i] = ango[par]*beam["pz"][i];
  }
}

void Scattering::ComputeFactor(const double& thickness) {

  _factor = 13.6*sqrt(thickness)*(1+0.038*log(thickness));
}
