#include "EnergyLoss.hh"

EnergyLoss::EnergyLoss() :
  _mat(), _stepsize(1.), _thickness(0.) {
}

EnergyLoss::EnergyLoss(const Material& mat,
		       const double& thickness) :
  _mat(mat), _stepsize(1.), _thickness(thickness) {

  try {
    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "EnergyLoss::EnergyLoss"));
  }
}

EnergyLoss::EnergyLoss(const EnergyLoss& eloss) {
  *this = eloss;
}

EnergyLoss& EnergyLoss::operator=(const EnergyLoss& eloss) {
  if ( this == &eloss )
      return *this;

  _mat = eloss._mat;
  _stepsize = eloss._stepsize;
  _thickness = eloss._thickness;

  return *this;
}

EnergyLoss::~EnergyLoss () {}

double EnergyLoss::GetMomentumLoss(const double& p,
				   const double& m) const {

  // Get the input energy
  double E = sqrt(p*p+m*m);

  // Step through the volume, apply eloss, until exit or negative energy is reached
  double delz = 0.;
  while ( delz < _thickness ) {

    // Get the kinematic parameters
    double pi = sqrt(E*E-m*m);
    double bg = pi/m;

    // Compute the eloss and step
    E -= BetheBloch(bg)*_stepsize;
    delz += _stepsize;

    if ( E < m )
	return p;
  }

  // Return the output momentum
  return p-sqrt(E*E-m*m);
}

double EnergyLoss::GetMostProbableMomentumLoss(const double& p,
				   	       const double& m) const {

  // Get the input energy
  double E = sqrt(p*p+m*m);

  // Apply the most probable energy loss
  double bg = p/m;
  double b = bg/sqrt(1.+bg*bg);
  double xi = .5*_mat.k1*_thickness/(b*b);

  E -= xi*(log(_mat.k2*bg*bg)+log(xi*1e6/_mat.I)+0.2000-b*b*-0.);
  if ( E < m )
      E = m;
  
  // Return the output momentum
  return p-sqrt(E*E-m*m);
}

void EnergyLoss::IonizeBunch(std::map<std::string, std::vector<double>>& beam,
			    const double& p,
			    const double& m) const {

  // Check weather it is a 2 or 4D beam to scatter
  std::vector<std::string> moms = {"pz"};
  for (const std::string mom : {"px", "py"})
    if ( beam.find(mom) != beam.end() )
        moms.push_back(mom);
  if ( !moms.size() || beam.find("pz") == beam.end() )
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	    "Not enough components of the momentum provided",
	    "EnergyLoss::IonizeBunch"));

  // Reduce the momentum of all the samples
  double scale = 1. - GetMostProbableMomentumLoss(p, m)/p;
  size_t i;
  for (const std::string& par : moms)
    for (i = 0; i < beam.at(par).size(); i++)
	beam[par][i] *= scale;
}

double EnergyLoss::BetheBloch(const double& bg) const {
  return _mat.k1*(log(_mat.k2*bg*bg)*(1.+bg*bg)/(bg*bg)-1);
}

void EnergyLoss::Initialize() {

  // Check that all the required parameters have been provided
  if ( !_mat.name.size() )
      Pitch::print(Pitch::warning, "No material name specified", "EnergyLoss::Initialize");
  if ( !_mat.Z || !_mat.A || !_mat.rho || !_mat.I )
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	    "Missing essential material parameters",
	    "EnergyLoss::Initialize"));

  // Set the two constants to be used in the Bethe-Bloch
  double K = 0.307;  				// [MeV*cm^2/g]
  _mat.k1 = .1*K*_mat.Z*_mat.rho/_mat.A;	// [Mev/mm]
  _mat.k2 = 2*510998/_mat.I;

  // Set the common Material structure
  _mat = _mat;
}
