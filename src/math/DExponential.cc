#include "DExponential.hh"

// ROOT style function necessary to produce a 3D contour, must be at global level
double RExponentialOffset(double *x, double *par) {

  // Get the parameters
  size_t n = par[1];
  double pref = par[2]*pow(.5, n);

  // Compute the inside of the exponential, sum of the lambda_i*|x_i-mu_i|
  double sfabs(0);
  size_t i;
  for (i = 0; i < n; i++)
      sfabs += par[3+n+i]*fabs(x[i]-par[3+i]);

  return pref*exp(-sfabs)-par[0];
}

DExponential::DExponential(const size_t n) :
  _mean(n, 0.), _lambda(n, 1.), _prod(0.) {

  try {
    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "DExponential::DExponential"));
  }
}

DExponential::DExponential(const double mean,
		   	   const double lambda) :
  _mean(1, mean), _lambda(1, lambda), _prod(0.) {

  try {
    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "DExponential::DExponential"));
  }
}

DExponential::DExponential(const std::vector<double>& mean,
	     	   	   const std::vector<double>& lambda) :
  _mean(mean), _lambda(lambda), _prod(0.) {

  try {
    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "DExponential::DExponential"));
  }
}

DExponential::DExponential(const DExponential& df) {
  *this = df;
}

DExponential& DExponential::operator=(const DExponential& df) {
  if ( this == &df )
      return *this;

  _mean = df._mean;
  _lambda = df._lambda;
  _prod = df._prod;

  return *this;
}

DExponential::~DExponential () {}

double DExponential::CDF(const double& v) const {

  // CDF only defined for 1D
  Assert::IsEqual("DExponential::CDF", "Dimension", _dim, (size_t)1);
  if ( v < _mean[0] )
      return .5*exp(_lambda[0]*(v-_mean[0]));
  return 1.-.5*exp(-_lambda[0]*(v-_mean[0]));
}

double DExponential::CDFRadial(const double& R) const {

  // The exponential radial CDF follows the lower incomplete Gamma function
  Assert::IsGreater("DExponential::CDFRadial", "Radius", R, 0.);
  return TMath::Gamma(_dim, R);
}

std::vector<TLine*> DExponential::Contour1D(const double alpha) const {

  Assert::IsProbability("DExponential::Contour1D", "alpha", alpha);
  double R = Radius(alpha)/_lambda[0];
  double level = Level(alpha);
  TLine* line = new TLine(_mean[0]-R, level, _mean[0]+R, level);
  line->SetLineWidth(2);
  return {line};
}

std::vector<TObject*> DExponential::Contour2D(const double alpha) const {

  Assert::IsProbability("DExponential::Contour2D", "alpha", alpha);
  double R = Radius(alpha);
  double R0(R/_lambda[0]), R1(R/_lambda[1]);
  std::vector<double> x({_mean[0], _mean[0]-R0, _mean[0], _mean[0]+R0, _mean[0]}),
		      y({_mean[1]-R1, _mean[1], _mean[1]+R1, _mean[1], _mean[1]-R1});
  TPolyLine* polygon = new TPolyLine(x.size(), &(x[0]), &(y[0]));
  polygon->SetLineWidth(2);
  return {polygon};
}

TF3* DExponential::Contour3D(const double alpha) const {

  // Find the level of the isosurface for the provided alpha
  Assert::IsProbability("DExponential::Contour3D", "alpha", alpha);
  double C = Level(alpha);

  // Define and return the offset function
  std::vector<double> pars = OffsetParameters(C);
  TF3* octahedron = new TF3("octahedron", RExponentialOffset, _lower[0], _upper[0],
		      	    _lower[1], _upper[1], _lower[2], _upper[2], pars.size());
  octahedron->SetTitle(_title.c_str());
  octahedron->SetParameters(&(pars[0]));
  octahedron->SetFillColorAlpha(1, .25);
  octahedron->SetLineColor(1);
  return octahedron;
}

std::vector<TLine*> DExponential::ContourRadial(const double alpha) const {

  Assert::IsProbability("DExponential::ContourRadial", "alpha", alpha);
  double R = Radius(alpha);
  double level = Radial(R);
  TLine* line = new TLine(0., level, R, level);
  line->SetLineWidth(2);
  return {line};
}

double DExponential::ContourVolume(const double alpha) const {

  // The contour volume is the volume of an n-orthoplexe of half axis R, scaled
  Assert::IsProbability("DExponential::ContourVolume", "alpha", alpha);
  return pow(2*Radius(alpha), _dim)/Math::Factorial(_dim)/_prod;
}

double DExponential::Evaluate(const std::vector<double>& v) const {

  Assert::SameSize("DExponential::Evaluate", "Function space and argument", _lower, v);

  // Compute the inside of the exponential, sum of the lambda_i*|x_i-mu_i|
  double sfabs(0);
  size_t i;
  for (i = 0; i < _dim; i++)
      sfabs += _lambda[i]*fabs(v[i]-_mean[i]);

  return pow(.5, _dim)*_prod*exp(-sfabs);
}

double DExponential::Level(const double alpha) const {

  Assert::IsProbability("DExponential::Level", "alpha", alpha);
  double R = Radius(alpha);
  return Radial(R);
}

double DExponential::Norm() const {

  double norm = 1.;
  size_t i;
  for (i = 0; i < _dim; i++)
      norm *= Norm1D(i);

  return norm;
}

double DExponential::Radial(const double R) const {

  Assert::IsGreater("DExponential::Radial", "Radius", R, 0.);
  return pow(.5, _dim)*_prod*exp(-R);
}

double DExponential::Radius(const double alpha) const {

  Assert::IsProbability("DExponential::Radius", "alpha", alpha);

  // Clear the easy case of n = 1
  if ( _dim == 1 )
      return -log(1.-alpha);

  // If not, need to numerically compute the radius by inverting the LI Gamma function
  TF1 fcont("", "TMath::Gamma([0], x)", 0, 100);
  fcont.SetParameter(0, _dim);
  return fcont.GetX(alpha);
}

double DExponential::Random() {

  // Randomize the sign, use ROOT to sample from a one-sided exponential
  int sign = (_rdmzer.Rndm() > .5) ? 1 : -1;
  return sign*_rdmzer.Exp(1./_lambda[0]);
}

std::vector<double> DExponential::RandomVector() {

  // Generate a vector of independant exponentially distributed variables,
  // offset by the mean and scale by the lambda vector
  std::vector<double> v(_dim);
  int sign;
  size_t i;
  for (i = 0; i < _dim; i++) {
    sign = (_rdmzer.Rndm() > .5) ? 1 : -1;
    v[i] = _mean[i] + sign*_rdmzer.Exp(1./_lambda[i]);
  }

  return v;
}

bool DExponential::Initialize() {

  // Check the validity of the input
  Assert::IsNonZero("DExponential::Initialize", "Dimension", _mean.size());
  Assert::SameSize("DExponential::Initialize", "Mean and lambda vectors", _mean, _lambda);
  for (size_t i = 0; i < _lambda.size(); i++)
      Assert::IsGreater("DExponential::Initialize", "lambda["+std::to_string(i)+"]",
			_lambda[i], 0., true);

  // Get the dimension of the space
  _dim = _mean.size();

  // Initialize the default name of the function and its title
  _name = "exp"+std::to_string(_dim);
  _title = std::to_string(_dim)+"D Exponential";

  // Compute the product of the scale factors
  _prod = 1.;
  size_t i;
  for(i = 0; i < _dim; i++)
      _prod *= _lambda[i];

  // Set up a default range around the mean (+/-5 taus)
  for (i = 0; i < _dim; i++)
      _lower.push_back(_mean[i]-5./_lambda[i]);
  for (i = 0; i < _dim; i++)
      _upper.push_back(_mean[i]+5./_lambda[i]);

  // Initialize the pseudorandom number generator
  time_t tt;
  _rdmzer = TRandom3(time(&tt));

  return true;
}

double DExponential::Norm1D(const size_t i) const {

  // Compute the normalised radii
  double Rl = _lambda[i]*fabs(_lower[i]-_mean[i]);
  double Ru = _lambda[i]*fabs(_upper[i]-_mean[i]);

  // Return the normalisation
  if ( _upper[i] < _mean[i] || _lower[i] > _mean[i]  )
      return .5*fabs(exp(-Ru) - exp(-Rl));

  return 1.-.5*(exp(-Rl) + exp(-Ru));
}

std::vector<double> DExponential::OffsetParameters(const double C) const {

  std::vector<double> pars;
  pars.push_back(C);
  pars.push_back(_dim);
  pars.push_back(_prod);
  size_t i;
  for (i = 0; i < _dim; i++)
      pars.push_back(_mean[i]);
  for (i = 0; i < _dim; i++)
      pars.push_back(_lambda[i]);

  return pars;
}
