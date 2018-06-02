#include "DChiSquared.hh"

double RChiSquared(double *x, double *par) {

  return pow(2, -par[0]/2)*pow(x[0], par[0]/2-1)*exp(-x[0]/2)/tgamma(par[0]/2);
}

double RChiSquaredPContent(double *x, double *par) {

  // For the input level, find the corresponding radii and compute the p-value
  double xmax = par[0] > 2 ? par[0]-2 : 0;
  TF1 f("", RChiSquared, std::max(0., xmax-10*par[0]), par[0]+10*par[0], 1);
  f.SetParameter(0, par[0]);
  par[1] = f.GetX(x[0], std::max(0., xmax-10*par[0]), xmax);
  par[2] = f.GetX(x[0], xmax, par[0]+10*par[0]);

  return TMath::Gamma(par[0]/2, par[2]/2)-TMath::Gamma(par[0]/2, par[1]/2);
}

DChiSquared::DChiSquared(const size_t ndf) :
  _ndf(ndf), _frand() {

  try {
    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "DChiSquared::DChiSquared"));
  }
}

DChiSquared::DChiSquared(const DChiSquared& df) {
  *this = df;
}

DChiSquared& DChiSquared::operator=(const DChiSquared& df) {
  if ( this == &df )
      return *this;

  _ndf = df._ndf;
  _frand = df._frand;

  return *this;
}

DChiSquared::~DChiSquared () {}

double DChiSquared::CDF(const double& v) const {

  return TMath::Gamma((double)_ndf/2, v/2);
}

double DChiSquared::CDFRadial(const double& R) const {

  throw(Exceptions::Exception(Exceptions::recoverable,
	"No radial symmetry for any definition of the p-norm of a vector",
	"DChiSquared::CDFRadial"));
  return 0.;
}

std::vector<TLine*> DChiSquared::Contour1D(const double alpha) const {

  Assert::IsProbability("DChiSquared::Contour1D", "alpha", alpha);
  double r, R, l;
  Radii(alpha, r, R, l);
  TLine* line = new TLine(r, l, R, l);
  line->SetLineWidth(2);
  return {line};
}

std::vector<TObject*> DChiSquared::Contour2D(const double alpha) const {

  // The Chi-squared distribution is only defined in 1D as it is a function of matrices in nD
  // If handled, it is possible to skip the ignore the empty output, will not crash
  throw(Exceptions::Exception(Exceptions::recoverable,
	"The Chi-squared distribution has no vectorial 2D extension",
	"DChiSquared::Contour2D"));
  return {new TObject()};
}

TF3* DChiSquared::Contour3D(const double alpha) const {

  // The Chi-squared distribution is only defined in 1D as it is a function of matrices in nD
  // If handled, it is possible to skip the ignore the empty output, will not crash
  throw(Exceptions::Exception(Exceptions::recoverable,
	"The Chi-squared distribution has no vectorial 3D extension",
	"DChiSquared::Contour3D"));
  return new TF3();
}

std::vector<TLine*> DChiSquared::ContourRadial(const double alpha) const {

  throw(Exceptions::Exception(Exceptions::recoverable,
	"No radial symmetry for any definition of the p-norm of a vector",
	"DChiSquared::ContourRadial"));
  TLine* line = new TLine();
  return {line};
}

double DChiSquared::ContourVolume(const double alpha) const {

  Assert::IsProbability("DChiSquared::ContourVolume", "alpha", alpha);
  double r, R, l;
  Radii(alpha, r, R, l);
  return R-r;
}

double DChiSquared::Evaluate(const std::vector<double>& v) const {

  Assert::SameSize("DChiSquared::Evaluate", "Function space and argument", _lower, v);

  if ( v[0] < 0 )
      return 0.;

  double n = _ndf;
  return pow(2, -n/2)*pow(v[0], n/2-1)*exp(-v[0]/2)/tgamma(n/2);
}

double DChiSquared::Level(const double alpha) const {

  Assert::IsProbability("DChiSquared::Level", "alpha", alpha);
  double r, R, l;
  Radii(alpha, r, R, l);
  return l;
}

double DChiSquared::Norm() const {

  if ( _lower[0] <= 0. && _upper[0] <= 0. )
      return 0.;

  double n = _ndf;
  if ( _lower[0] <= 0. )
      return TMath::Gamma(n/2, _upper[0]/2);

  return TMath::Gamma(n/2, _upper[0]/2)-TMath::Gamma(n/2, _lower[0]/2);
}

double DChiSquared::Radial(const double R) const {

  throw(Exceptions::Exception(Exceptions::recoverable,
	"No radial symmetry for any definition of the p-norm of a vector",
	"DChiSquared::ContourRadial"));
  return 0.;
}

void DChiSquared::Radii(const double alpha, double& r, double& R, double& l) const {

  // Cannot compute analytically. Compute the radii by optimizing the p-value
  Assert::IsProbability("DChiSquared::Radii", "alpha", alpha);
  std::vector<double> pars = {(double)_ndf, r, R};
  double n = _ndf;
  double max = (_ndf > 2) ? pow(2, -n/2)*pow((n-2)/exp(1), n/2-1)/tgamma(n/2) : 1./(n+1);
  TF1 fcont("", RChiSquaredPContent, 0, max, 3);
  fcont.SetParameters(&(pars[0]));
  l = fcont.GetX(alpha, 0, max);
  r = fcont.GetParameter(1);
  R = fcont.GetParameter(2);
}

double DChiSquared::Random() {

  // Faster to use ROOT to sample from it than any other technique, analytically challenging
  return _frand.GetRandom();
}

std::vector<double> DChiSquared::RandomVector() {

  // Generate a 1-vector of a Chi-squared random variable
  return {Random()};
}

bool DChiSquared::Initialize() {

  // Check that the number of degrees of freedom is at least 1
  Assert::IsGreater("DChiSquared::Initialize", "ndf", _ndf, size_t(1));

  // Get the dimension of the space, always 1
  _dim = 1;

  // Initialize the default name of the function and its title
  _name = "chisq"+std::to_string(_ndf);
  _title = "Chi-squared Distrubution ("+std::to_string(_ndf)+" DOF)";

  // Set up a default range around the mean (_ndf+/-5*sqrt(_ndf))
  _lower.push_back(std::max(0., _ndf-5*sqrt(_ndf)));
  _upper.push_back(_ndf+5*sqrt(_ndf));

  // Initialize the pseudorandom number generators
  time_t tt;
  _rdmzer = TRandom3(time(&tt));

  if ( gRandom )
     delete gRandom;
  gRandom = new TRandom3(time(&tt));

  _frand = TF1("", RChiSquared, _lower[0], _upper[0], 1);
  _frand.SetParameter(0, _ndf);

  return true;
}
