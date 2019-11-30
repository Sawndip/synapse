#include "EnergyLoss.hh"

// ROOT style delta function for the density effects
// par[0]=a, par[1]=k, par[2]=X0, par[3]=X1, par[4]=C, par[5]=delta0
double RDeltaFunction(double* x, double* par) {

  double X = log(x[0])/log(10);
  if ( X < par[2] ) {
    return par[5]*pow(10., 2*(X-par[2]));
  } else if ( X < par[3] ) {
    return 2*log(10)*X - par[4] + par[0]*pow(par[3]-X, par[1]);
  } else {
    return 2*log(10)*X - par[4];
  }
}

// ROOT style Bethe-Bloch function
// par[0]=k1, par[1]=k2, par[2-7]=density function fitting parameters
double RBetheBloch(double* x, double* par) {

  double b2 = (x[0]*x[0])/(1.+x[0]*x[0]);
  return (par[0]/b2)*(log(par[1]*x[0]*x[0])-b2-RDeltaFunction(x, &par[2])/2);
}

// ROOT style stopping function
// par[0]=number of elements, par[(1+n*8)-(8+n*8)]=material n (0--N)
double RStoppingPower(double* x, double* par) {

  size_t n = par[0];
  double total(0.);
  for (size_t i = 0; i < par[0]; i++)
       total += par[1+i]*RBetheBloch(x, &par[1+n+i*8]);
  return total;
}

// ROOT style integrand of the CSDA function
// par[0]=mass [MeV/c^2], par[1:] stopping power parameters
double RIntegrand(double* x, double* par) {

  double b = x[0]/sqrt(1.+x[0]*x[0]);
  return par[0]*b/RStoppingPower(x, &par[1]);
}

// ROOT style CSDA function
// par[0]=number of parameters, par[1]=initial &beta;&gamma; par[2:]=integrand parameters
double RCSDA(double* x, double* par) {

  TF1 fint("int", RIntegrand, 0, 100, par[0]);
  fint.SetParameters(&par[2]);
  //return fint.Integral(x[0], par[1], &par[2], 1e-3);
  return fint.Integral(x[0], par[1], 1e-3);
}

EnergyLoss::EnergyLoss() :
  _mat(), _thickness(0.), _params(), _Kstar(0.) {
}

EnergyLoss::EnergyLoss(const Material& mat,
		       const double& thickness) :
  _mat(), _thickness(thickness), _params(), _Kstar(0.) {

  try {
    SetMaterial(mat);
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
  _thickness = eloss._thickness;
  _params = eloss._params;

  return *this;
}

EnergyLoss::~EnergyLoss () {}

double EnergyLoss::GetMomentumLoss(const double& p,
				   const double& m) const {

  // Get the initial &beta;&gamma;
  double bg0 = p/m;

  // Define the CSDA function, find the bg that suits
  std::vector<double> params = {(double)_params.size(), bg0, m};
  params.insert(params.end(), _params.begin(), _params.end());
  TF1 fcsda("csda", RCSDA, 0, 100, params.size());
  fcsda.SetParameters(&params[0]);

  // Find the final &beta;&gamma; that yields the required material thickness as the CSDA range
  double bgf = fcsda.GetX(_thickness, 0., bg0);

  // Return the output momentum
  return p-bgf*m;
}

double EnergyLoss::GetApproxMomentumLoss(const double& p,
			     	         const double& m) const {

  // Get the input energy and C factor
  double E = sqrt(p*p+m*m);
  double C = (E*E+m*m)/E;

  // Return the momentum is the particle is stopped
  double diff = C-_Kstar*_thickness;
  if ( diff < 2*m )
      return p;

  // Compute the final energy
  double Ex = diff*(1.+sqrt(1.-pow(2*m/diff, 2)))/2.;

  // Return the momentum loss
  return p-sqrt(Ex*Ex-m*m);
}

double EnergyLoss::GetMostProbableMomentumLoss(const double& p,
				   	       const double& m) const {

  // Most probable not implemented for composite materials, throw if mixture
  if ( _mat.number != 1 )
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	    "Not supported for composite materials",
	    "EnergyLoss::GetMostProbableMomentumLoss"));

  // Get the input energy
  double E = sqrt(p*p+m*m);

  // Apply the most probable energy loss
  const Element &el = _mat.elements[0];
  double bg = p/m;
  double b2 = bg*bg/(1.+bg*bg);
  double xi = .5*.1*0.307*el.ZoA*el.rho*_thickness/b2;

  E -= xi*(log(2*510998*bg*bg/el.I)+log(xi*1e6/el.I)+0.2000-b2-DeltaFunction(bg, el));
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
  double scale = 1. - GetMomentumLoss(p, m)/p;
  size_t i;
  for (const std::string& par : moms)
    for (i = 0; i < beam.at(par).size(); i++)
	beam[par][i] *= scale;
}

double EnergyLoss::StoppingPower(const double& bg) const {

  double total = 0.;
  for (size_t i = 0; i < _mat.number; i++)
      total += _mat.fractions[i]*BetheBloch(bg, _mat.elements[i]);

  return total;
}

double EnergyLoss::BetheBloch(const double& bg,
			      const Element& el) const {

  double b2 = (bg*bg)/(1.+bg*bg);
  return (.1*0.307*el.ZoA*el.rho/b2)*(log(2*510998*bg*bg/el.I)-b2-DeltaFunction(bg, el)/2);
}

double EnergyLoss::DeltaFunction(const double& bg,
				 const Element& el) const {

  double X = log(bg)/log(10);
  if ( X < el.X0 ) {
    return el.delta0*pow(10., 2*(X-el.X0));
  } else if ( X < el.X1 ) {
    return 2*log(10)*X - el.C + el.a*pow(el.X1-X, el.k);
  } else {
    return 2*log(10)*X - el.C;
  }
}

void EnergyLoss::SetMaterial(const Material& mat) {

  // Check that there is a material name (optional, simply warn)
  if ( !mat.name.size() )
      Pitch::print(Pitch::warning, "No material name specified", "EnergyLoss::SetMaterial");

  // Check that elements compose the material
  if ( !mat.number )
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	    "No element specified in the definition of the material",
	    "EnergyLoss::SetMaterial"));

  // Deep copy
  _mat = mat;

  // Set the ROOT parameters
  _params.push_back(_mat.number);
  for (const double& frac : _mat.fractions)
    _params.push_back(frac);

  for (const Element& el : _mat.elements) {
    _params.push_back(.1*0.307*el.ZoA*el.rho);
    _params.push_back(2*510998/el.I);
    _params.push_back(el.a);
    _params.push_back(el.k);
    _params.push_back(el.X0);
    _params.push_back(el.X1);
    _params.push_back(el.C);
    _params.push_back(el.delta0);
  }

  // Set the approximate prefactor
  _Kstar = 0.;
  for (size_t i = 0; i < _mat.number; i++) {
    const Element &el = _mat.elements[i];
    _Kstar += mat.fractions[i]*.1*0.307*el.ZoA*el.rho*log(2*510998/el.I);
  }
}

void EnergyLoss::DrawBetheBloch() const {

  // Define the function
  TF1 func("func", RStoppingPower, 0.1, 1e2-1, _params.size());
  func.SetTitle(";#beta#gamma;[MeV/mm]");
  func.SetParameters(&_params[0]);

  // Draw and save
  TCanvas *c = new TCanvas("c", "c", 1200, 800);
  gPad->SetLogx();
  gPad->SetLogy();
  func.Draw();
  c->SaveAs((_mat.name+"_energy_loss.pdf").c_str());
  delete c;
}
