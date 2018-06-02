#include "DPoisson.hh"

double RPoisson(double *x, double *par) {

  return pow(x[0], par[0])*exp(-x[0])/tgamma(par[0]);
}

double RPoissonPContent(double *x, double *par) {

  // For the input level, find the corresponding radii and compute the p-value
  TF1 f("", RPoisson, std::max(0., par[0]-10*sqrt(par[0])), par[0]+10*sqrt(par[0]), 1);
  f.SetParameter(0, par[0]);
  par[1] = f.GetX(x[0], std::max(0., par[0]-10*sqrt(par[0])), par[0]);
  par[2] = f.GetX(x[0], par[0], par[0]+10*par[0]);

  return TMath::Gamma(par[0]+1, par[2])-TMath::Gamma(par[0]+1, par[1]);
}

DPoisson::DPoisson() :
  _mean(0) {

  this->Initialize();
}

DPoisson::DPoisson(const double mu) :
  _mean(mu) {

  try {
    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "DPoisson::DPoisson"));
  }
}

DPoisson::DPoisson(const DPoisson& df) {
  *this = df;
}

DPoisson& DPoisson::operator=(const DPoisson& df) {
  if ( this == &df )
      return *this;

  _mean = df._mean;

  return *this;
}

DPoisson::~DPoisson () {}

double DPoisson::CDF(const double& v) const {

  if ( v <= 0 )
      return 0.;
  return TMath::Gamma(_mean+1, v);
}

double DPoisson::CDFRadial(const double& R) const {

  throw(Exceptions::Exception(Exceptions::recoverable,
	"No radial symmetry for any definition of the p-norm of a vector",
	"DPoisson::CDFRadial"));
  return 0.;
}

std::vector<TLine*> DPoisson::Contour1D(const double alpha) const {

  Assert::IsProbability("DPoisson::Contour1D", "alpha", alpha);
  double r, R, l;
  Radii(alpha, r, R, l);
  TLine* line = new TLine(r, l, R, l);
  line->SetLineWidth(2);
  return {line};
}

std::vector<TObject*> DPoisson::Contour2D(const double alpha) const {

  // 2D Poisson distributions are not normalisable and as such do not exist, throw
  // If handled, it is possible to skip the ignore the empty output, will not crash
  throw(Exceptions::Exception(Exceptions::recoverable,
	"2D extension of the Poisson distribution does not exist",
	"DPoisson::Contour2D"));
  return {new TObject()};
}

TF3* DPoisson::Contour3D(const double alpha) const {

  // 3D Poisson distributions are not normalisable and as such do not exist, throw
  // If handled, it is possible to skip the ignore the empty output, will not crash
  throw(Exceptions::Exception(Exceptions::recoverable,
	"3D extension of the Poisson distribution does not exist",
	"DPoisson::Contour3D"));
  return new TF3();
}

std::vector<TLine*> DPoisson::ContourRadial(const double alpha) const {

  throw(Exceptions::Exception(Exceptions::recoverable,
	"No radial symmetry for any definition of the p-norm of a vector",
	"DPoisson::ContourRadial"));
  TLine* line = new TLine();
  return {line};
}

double DPoisson::ContourVolume(const double alpha) const {

  Assert::IsProbability("DPoisson::ContourVolume", "alpha", alpha);
  double r, R, l;
  Radii(alpha, r, R, l);
  return R-r;
}

double DPoisson::Evaluate(const std::vector<double>& v) const {

  Assert::SameSize("DPoisson::Evaluate", "Function space and argument", _lower, v);
  if ( v[0] <= 0 )
      return 0.;

  return pow(v[0], _mean)*exp(-v[0])/tgamma(_mean);
}

double DPoisson::Level(const double alpha) const {

  Assert::IsProbability("DPoisson::Level", "alpha", alpha);
  double r, R, l;
  Radii(alpha, r, R, l);
  return l;
}

double DPoisson::Norm() const {

  if ( _lower[0] <= 0. && _upper[0] <= 0. )
      return 0.;

  if ( _lower[0] <= 0. )
      return TMath::Gamma(_mean+1, _upper[0]);

  return TMath::Gamma(_mean+1, _upper[0])-TMath::Gamma(_mean+1, _lower[0]);
}

double DPoisson::Radial(const double R) const {

  throw(Exceptions::Exception(Exceptions::recoverable,
	"No radial symmetry for any definition of the p-norm of a vector",
	"DPoisson::Radial"));
  return 0.;
}

void DPoisson::Radii(const double alpha, double& r, double& R, double& l) const {

  // Cannot compute analytically. Compute the radii by optimizing the p-value
  Assert::IsProbability("DPoisson::Radii", "alpha", alpha);
  std::vector<double> pars = {_mean, r, R};
  double max = !_mean ? 1 : exp(_mean*(log(_mean)-1))/tgamma(_mean+1);
  TF1 fcont("", RPoissonPContent, 0, max, 3);
  fcont.SetParameters(&(pars[0]));
  l = fcont.GetX(alpha, 0, max);
  r = fcont.GetParameter(1);
  R = fcont.GetParameter(2);
}

double DPoisson::Random() {

  // Use the ROOT randomizer
  return _rdmzer.Poisson(_mean);
}

std::vector<double> DPoisson::RandomVector() {

  // Generate a 1-vector of a Poisson random variable
  return {Random()};
}

bool DPoisson::Initialize() {

  // Assert that the mean is strictly greater than 0
  Assert::IsGreater("DPoisson::Initialize", "mean", _mean, 0.);

  // Get the dimension of the space, always 1
  _dim = 1;

  // Initialize the default name of the function and its title
  _name = "poisson";
  _title = "Poisson Distrubution";

  // Set up a default range around the mean (mean+/-5*sqrt(mean))
  _lower.push_back(std::max(0., _mean-5*sqrt(_mean)));
  _upper.push_back(_mean+5*sqrt(_mean));

  // Initialize the pseudorandom number generator
  time_t tt;
  _rdmzer = TRandom3(time(&tt));

  return true;
}
