#include "DTriangular.hh"

// ROOT style function necessary to numerically compute the radius of a given p-value
double RTriangularPContent(double *x, double *par) {

  if ( x[0] > 1. )
      return 1.;
  if ( par[0] == 1 )
      return 2*x[0]-x[0]*x[0];

  double n(par[0]);
  return pow(x[0], n)*((n+1) - n*x[0]);
}

// ROOT style function necessary to produce a 3D contour, must be at global level
double RTriangularOffset(double *x, double *par) {

  // Get the parameters
  size_t n = par[1];
  double pref = par[2]*Math::Factorial(n+1)*pow(.5, n);

  // Compute the inside of the exponential, sum of the lambda_i*|x_i-mu_i|
  double sfabs(0);
  size_t i;
  for (i = 0; i < n; i++)
      sfabs += par[3+n+i]*fabs(x[i]-par[3+i]);

  return pref*(1.-sfabs)-par[0];
}

DTriangular::DTriangular(const size_t n) :
  _mean(n, 0.), _lambda(n, 1.), _prod(0.) {

  try {
    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "DTriangular::DTriangular"));
  }
}

DTriangular::DTriangular(const double mean,
		   	 const double lambda) :
  _mean(1, mean), _lambda(1, lambda), _prod(0.) {

  try {
    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "DTriangular::DTriangular"));
  }
}

DTriangular::DTriangular(const std::vector<double>& mean,
	     	   	 const std::vector<double>& lambda) :
  _mean(mean), _lambda(lambda), _prod(0.) {

  try {
    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "DTriangular::DTriangular"));
  }
}

DTriangular::DTriangular(const DTriangular& df) {
  *this = df;
}

DTriangular& DTriangular::operator=(const DTriangular& df) {
  if ( this == &df )
      return *this;

  _mean = df._mean;
  _lambda = df._lambda;
  _prod = df._prod;

  return *this;
}

DTriangular::~DTriangular () {}

double DTriangular::CDF(const double& v) const {

  // CDF only defined for 1D
  Assert::IsEqual("DTriangular::CDF", "Dimension", _dim, (size_t)1);
  double z = _lambda[0]*(v-_mean[0]);
  if ( z <= -1. )
      return 0.;
  if ( z >= 1. )
      return 1.;
  if ( z <= 0 )
      return .5*pow(1+z, 2);

  return 1.-.5*pow(1-z, 2);
}

double DTriangular::CDFRadial(const double& R) const {

  // The triangular radial CDF follows a polynomial function
  Assert::IsGreater("DTriangular::CDFRadial", "Radius", R, 0.);
  if ( R > 1 )
      return 1.;
  return pow(R, _dim)*((_dim+1) - _dim*R);
}

std::vector<TLine*> DTriangular::Contour1D(const double alpha) const {

  Assert::IsProbability("DTriangular::Contour1D", "alpha", alpha);
  double R = Radius(alpha)/_lambda[0];
  double level = Level(alpha);
  TLine* line = new TLine(_mean[0]-R, level, _mean[0]+R, level);
  line->SetLineWidth(2);
  return {line};
}

std::vector<TObject*> DTriangular::Contour2D(const double alpha) const {

  Assert::IsProbability("DTriangular::Contour2D", "alpha", alpha);
  double R = Radius(alpha);
  double R0(R/_lambda[0]), R1(R/_lambda[1]);
  std::vector<double> x({_mean[0], _mean[0]-R0, _mean[0], _mean[0]+R0, _mean[0]}),
		      y({_mean[1]-R1, _mean[1], _mean[1]+R1, _mean[1], _mean[1]-R1});
  TPolyLine* polygon = new TPolyLine(x.size(), &(x[0]), &(y[0]));
  polygon->SetLineWidth(2);
  return {polygon};
}

TF3* DTriangular::Contour3D(const double alpha) const {

  // Find the level of the isosurface for the provided alpha
  Assert::IsProbability("DTriangular::Contour3D", "alpha", alpha);
  double C = Level(alpha);

  // Define and return the offset function
  std::vector<double> pars = OffsetParameters(C);
  TF3* nrhombus = new TF3("nrhombus", RTriangularOffset, _lower[0], _upper[0],
		      	  _lower[1], _upper[1], _lower[2], _upper[2], pars.size());
  nrhombus->SetTitle(_title.c_str());
  nrhombus->SetParameters(&(pars[0]));
  nrhombus->SetFillColorAlpha(1, .25);
  nrhombus->SetLineColor(1);
  return nrhombus;
}

std::vector<TLine*> DTriangular::ContourRadial(const double alpha) const {

  Assert::IsProbability("DTriangular::ContourRadial", "alpha", alpha);
  double R = Radius(alpha);
  double level = Radial(R);
  TLine* line = new TLine(0., level, R, level);
  line->SetLineWidth(2);
  return {line};
}

double DTriangular::ContourVolume(const double alpha) const {

  Assert::IsProbability("DTriangular::ContourVolume", "alpha", alpha);

  // The contour volume is the volume of an n-orthoplexe of half axis R, scaled
  return pow(2*Radius(alpha), _dim)/Math::Factorial(_dim)/_prod;
}

double DTriangular::Evaluate(const std::vector<double>& v) const {

  Assert::SameSize("DTriangular::Evaluate", "Function space and argument", _lower, v);

  // Compute the argument of the exponential, sum of the lambda_i*|x_i-mu_i|
  double sfabs(0);
  size_t i;
  for (i = 0; i < _dim; i++)
      sfabs += _lambda[i]*fabs(v[i]-_mean[i]);

  if ( sfabs > 1 )
      return 0;

  return Math::Factorial(_dim+1)*pow(.5, _dim)*_prod*(1.-sfabs);
}

double DTriangular::Level(const double alpha) const {

  Assert::IsProbability("DTriangular::Level", "alpha", alpha);
  double R = Radius(alpha);
  return Radial(R);
}

double DTriangular::Norm() const {

  // The distribution is always fully contained in this configuration
  return 1.;
}

double DTriangular::Radial(const double R) const {

  Assert::IsGreater("DTriangular::Radial", "Radius", R, 0.);
  if ( R > 1 )
      return 0.;

  return Math::Factorial(_dim+1)*pow(.5, _dim)*_prod*(1.-R);
}

double DTriangular::Radius(const double alpha) const {

  Assert::IsProbability("DTriangular::Radius", "alpha", alpha);

  // Clear the easy case of n = 1
  if ( _dim == 1 )
      return 1.-sqrt(1.-alpha);

  // If not, need to numerically compute the radius by inverting the P-value function
  TF1 fcont("", RTriangularPContent, 0, 1, 1);
  fcont.SetParameter(0, _dim);
  return fcont.GetX(alpha);
}

double DTriangular::Random() {

  // Simply find the radius that corresponds to the random variable
  // Randomize the sign, apply and return
  int sign = (_rdmzer.Rndm() > .5) ? 1 : -1;
  return _mean[0] + sign*Radius(_rdmzer.Rndm())/_lambda[0];
}

std::vector<double> DTriangular::RandomVector() {

  // One cannot simply generate independant variables as the condition that
  // sum_i |x_i| < 1 correlates the variables. For 1D and 2D, it is trivial
  std::vector<double> v(_dim);
  if ( _dim == 1 ) {
    return {Random()};
  } else if ( _dim == 2 ) {
    // First compute the L1 radius at which the random point lies
    double R = Radius(_rdmzer.Rndm());

    // The first RV, v0, is uniformly distributed between -R and R and v2 = +/-(R-v1)
    int sign = (_rdmzer.Rndm() > .5) ? 1 : -1;
    double ran = _rdmzer.Uniform(-R, R);
    v[0] = _mean[0] + ran/_lambda[0];
    v[1] = _mean[1] + sign*(R-fabs(ran))/_lambda[1];
    return v;
  }

  // For higher dimensions, simply pick n points from an n-exponential (L1 symmetric)
  double R = Radius(_rdmzer.Rndm());
  size_t i;
  int sign;
  double norm(0.);
  for (i = 0; i < _dim; i++) {
    sign = (_rdmzer.Rndm() > .5) ? 1 : -1;
    v[i] = sign*_rdmzer.Exp(1.);
    norm += fabs(v[i]);
  }

  for (i = 0; i < _dim; i++)
      v[i] *= R/norm;    
  return v;
}

bool DTriangular::Initialize() {

  // Check the the dimension is at least 1 and that lambdas are strictly greater than 0
  size_t i;
  Assert::IsNonZero("DTriangular::Initialize", "Dimension", _mean.size());
  Assert::SameSize("DTriangular::Initialize", "Mean and lambda vectors", _mean, _lambda);
  for (i = 0; i < _lambda.size(); i++)
      Assert::IsGreater("DTriangular::Initialize",
				"lambda["+std::to_string(i)+"]", _lambda[i], 0., true);  

  // Get the dimension of the space
  _dim = _mean.size();

  // Initialize the default name of the function and its title
  _name = "tri"+std::to_string(_dim);
  _title = std::to_string(_dim)+"D Triangular";

  // Compute the product of the scale factors
  _prod = 1.;
  for(i = 0; i < _dim; i++)
      _prod *= _lambda[i];

  // Set up a default range in each dimension to fully contain the distribution (+/-10%)
  for (i = 0; i < _dim; i++)
      _lower.push_back(_mean[i]-1./_lambda[i]-.1*2./_lambda[i]);
  for (i = 0; i < _dim; i++)
      _upper.push_back(_mean[i]+1./_lambda[i]+.1*2./_lambda[i]);

  // Initialize the pseudorandom number generator
  time_t tt;
  _rdmzer = TRandom3(time(&tt));

  return true;
}

std::vector<double> DTriangular::OffsetParameters(const double C) const {

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
