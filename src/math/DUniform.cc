#include "DUniform.hh"

// ROOT style function necessary to produce a 3D contour, must be at global level
double RUniformOffset(double *x, double *par) {

  // Get the parameters
  size_t n = par[0];

  // Return negative value if outside of the offset box, positive otherwise
  size_t i;
  for (i = 0; i < n; i++)
    if ( x[i] < par[1+i] || x[i] > par[1+n+i] )
	return -1.;

  return 1.;
}

DUniform::DUniform(const size_t n) :
  _l(n, -1.), _u(n, 1.), _vol(pow(2., n)) {

  try {
    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "DUniform::DUniform"));
  }
}

DUniform::DUniform(const double l,
		   const double u) :
  _l(1, l), _u(1, u), _vol(u-l) {

  try {
    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "DUniform::DUniform"));
  }
}

DUniform::DUniform(const std::vector<double>& l,
	     	   const std::vector<double>& u) :
  _l(l), _u(u), _vol(0.) {

  try {
    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "DUniform::DUniform"));
  }
}

DUniform::DUniform(const DUniform& df) {
  *this = df;
}

DUniform& DUniform::operator=(const DUniform& df) {
  if ( this == &df )
      return *this;

  _l = df._l;
  _u = df._u;
  _vol = df._vol;

  return *this;
}

DUniform::~DUniform () {}

double DUniform::CDF(const double& v) const {

  // CDF only defined for 1D
  Assert::IsEqual("DUniform::CDF", "Dimension", _dim, (size_t)1);
  if ( v <= _l[0] )
      return 0.;
  if ( v >= _u[0] )
      return 1.;

  return (v-_l[0])/(_u[0]-_l[0]);
}

double DUniform::CDFRadial(const double& R) const {

  // The triangular radial CDF follows a linear function
  Assert::IsGreater("DUniform::CDFRadial", "Radius", R, 0.);
  if ( R > 1. )
      return 1.;

  return pow(R, _dim);
}

std::vector<TLine*> DUniform::Contour1D(const double alpha) const {

  Assert::IsProbability("DUniform::Contour1D", "alpha", alpha);
  double mean = (_u[0]+_l[0])/2;
  double length = _u[0]-_l[0];
  double max = Level(alpha);
  TLine* line = new TLine(mean-alpha*length/2., max, mean+alpha*length/2., max);
  line->SetLineWidth(2);
  return {line};
}

std::vector<TObject*> DUniform::Contour2D(const double alpha) const {

  Assert::IsProbability("DUniform::Contour2D", "alpha", alpha);
  double salpha = sqrt(alpha);
  TBox* box = new TBox(0.5*(1+salpha)*_l[0]+0.5*(1-salpha)*_u[0],
		       0.5*(1+salpha)*_l[1]+0.5*(1-salpha)*_u[1],
		       0.5*(1-salpha)*_l[0]+0.5*(1+salpha)*_u[0],
		       0.5*(1-salpha)*_l[1]+0.5*(1+salpha)*_u[1]);
  box->SetFillStyle(0);
  box->SetLineWidth(2);
  return {box};
}

TF3* DUniform::Contour3D(const double alpha) const {

  // Find the box scaling factor, cubic root of alpha
  Assert::IsProbability("DUniform::Contour3D", "alpha", alpha);
  double C = pow(alpha, 1./3);

  // Define and return the offset function
  std::vector<double> pars = OffsetParameters(C);
  TF3* cube = new TF3("cube", RUniformOffset, _lower[0], _upper[0],
		      _lower[1], _upper[1], _lower[2], _upper[2], pars.size());
  cube->SetTitle(_title.c_str());
  cube->SetParameters(&(pars[0]));
  cube->SetFillColorAlpha(1, .25);
  cube->SetLineColor(1);
  return cube;
}

std::vector<TLine*> DUniform::ContourRadial(const double alpha) const {

  Assert::IsProbability("DUniform::ContourRadial", "alpha", alpha);
  double max = Level(alpha);
  double radius = Radius(alpha);
  TLine* line = new TLine(0, max, radius, max);
  line->SetLineWidth(2);
  return {line};
}

double DUniform::ContourVolume(const double alpha) const {

  // The contour volume is alpha times the total volume
  Assert::IsProbability("DUniform::ContourVolume", "alpha", alpha);
  return alpha*_vol;
}

double DUniform::Evaluate(const std::vector<double>& v) const {

  Assert::SameSize("DUniform::Evaluate", "Function space and argument", _lower, v);

  // Return 0 if outside of the box, 1/vol otherwise
  size_t i;
  for (i = 0; i < _dim; i++)
    if ( v[i] < _l[i] || v[i] > _u[i] )
	return 0.;

  return 1./_vol;
}

double DUniform::Level(const double alpha) const {

  Assert::IsProbability("DUniform::Level", "alpha", alpha);
  return 1./_vol;
}

double DUniform::Norm() const {

  // The distribution is always fully contained in this configuration
  return 1.;
}

double DUniform::Radial(const double R) const {

  Assert::IsGreater("DUniform::Radial", "Radius", R, 0.);

  // Return 0 if the supremum is larger than half the range
  if ( R > 1. )
      return 0.;

  return 1./_vol;
}

double DUniform::Radius(const double alpha) const {

  Assert::IsProbability("DGaus::Radius", "alpha", alpha);
  return pow(alpha, 1./_dim);
}

double DUniform::Random() {

  return _rdmzer.Uniform(_l[0], _u[0]);
}

std::vector<double> DUniform::RandomVector() {

  // Generate a vector of independant uniformly distributed variables
  std::vector<double> v(_dim);
  size_t i;
  for (i = 0; i < _dim; i++)
      v[i] = _rdmzer.Uniform(_l[i], _u[i]);

  return v;
}

bool DUniform::Initialize() {

  // Check the validity of the input
  Assert::IsNonZero("DUniform::Initialize", "Dimension", _l.size());
  Assert::SameSize("DUniform::Initialize", "Lower and upper boundaries", _l, _u);
  for (size_t i = 0; i < _l.size(); i++)
    if ( _l[i] >= _u[i] )
        throw(Exceptions::Exception(Exceptions::nonRecoverable,
	      "The boundaries of axis "+std::to_string(i)+" define a zero or negative interval",
	      "DUniform::Initialize"));

  // Get the dimension of the space
  _dim = _l.size();

  // Initialize the default name of the function and its title
  _name = "uni"+std::to_string(_dim);
  _title = std::to_string(_dim)+"D Uniform";

  // Compute the volume of the n-orthotope
  _vol = 1.;
  size_t i;
  for (i = 0; i < _dim; i++)
      _vol *= (_u[i]-_l[i]);

  // Set up a default range in each dimension to fully contain the distribution (+/-10%)
  for (i = 0; i < _dim; i++)
      _lower.push_back(_l[i]-.1*(_u[i]-_l[i]));
  for (i = 0; i < _dim; i++)
      _upper.push_back(_u[i]+.1*(_u[i]-_l[i]));

  // Initialize the pseudorandom number generator
  time_t tt;
  _rdmzer = TRandom3(time(&tt));

  return true;
}

std::vector<double> DUniform::OffsetParameters(const double C) const {

  std::vector<double> pars;
  pars.push_back(_dim);
  size_t i;
  for (i = 0; i < _dim; i++)
      pars.push_back(.5*(1.+C)*_l[i]+.5*(1.-C)*_u[i]);
  for (i = 0; i < _dim; i++)
      pars.push_back(.5*(1.-C)*_l[i]+.5*(1.+C)*_u[i]);

  return pars;
}
