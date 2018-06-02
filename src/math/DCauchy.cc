#include "DCauchy.hh"

DCauchy::DCauchy() :
  _mean(0), _sigma(1.) {

  this->Initialize();
}

DCauchy::DCauchy(const double mu,
		 const double sigma) :
  _mean(mu), _sigma(sigma) {

  try {
    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "DCauchy::DCauchy"));
  }
}

DCauchy::DCauchy(const DCauchy& df) {
  *this = df;
}

DCauchy& DCauchy::operator=(const DCauchy& df) {
  if ( this == &df )
      return *this;

  _mean = df._mean;
  _sigma = df._sigma;

  return *this;
}

DCauchy::~DCauchy () {}

double DCauchy::CDF(const double& v) const {

  // CDF only defined for 1D
  Assert::IsEqual("DCauchy::CDF", "Dimension", _dim, (size_t)1);
  return .5+atan((v-_mean)/_sigma)/M_PI;
}

double DCauchy::CDFRadial(const double& R) const {

  // The Cauchy radial CDF combines both sides (factor 2)
  Assert::IsGreater("DCauchy::CDFRadial", "Radius", R, 0.);
  return 2*atan(R/_sigma)/M_PI;
}

std::vector<TLine*> DCauchy::Contour1D(const double alpha) const {

  Assert::IsProbability("DCauchy::Contour1D", "alpha", alpha);
  double R = _sigma*tan(M_PI*alpha/2);
  TLine* line = new TLine(-R, Evaluate({-R}), R, Evaluate({R}));
  line->SetLineWidth(2);
  return {line};
}

std::vector<TObject*> DCauchy::Contour2D(const double alpha) const {

  // 2D Cauchy distributions are not normalisable and as such do not exist, throw
  // If handled, it is possible to skip the ignore the empty output, will not crash
  throw(Exceptions::Exception(Exceptions::recoverable,
	"2D extension of the Cauchy distribution does not exist",
	"DCauchy::Contour2D"));
  return {new TObject()};
}

TF3* DCauchy::Contour3D(const double alpha) const {

  // 3D Cauchy distributions are not normalisable and as such do not exist, throw
  // If handled, it is possible to skip the ignore the empty output, will not crash
  throw(Exceptions::Exception(Exceptions::recoverable,
	"3D extension of the Cauchy distribution does not exist",
	"DCauchy::Contour3D"));
  return new TF3();
}

std::vector<TLine*> DCauchy::ContourRadial(const double alpha) const {

  Assert::IsProbability("DCauchy::ContourRadial", "alpha", alpha);
  double R = tan(M_PI*alpha/2);
  TLine* line = new TLine(0, Evaluate({R}), R, Evaluate({R}));
  line->SetLineWidth(2);
  return {line};
}

double DCauchy::ContourVolume(const double alpha) const {

  Assert::IsProbability("DCauchy::ContourVolume", "alpha", alpha);
  return 2*_sigma*tan(M_PI*alpha/2);
}

double DCauchy::Evaluate(const std::vector<double>& v) const {

  Assert::SameSize("DCauchy::Evaluate", "Function space and argument", _lower, v);
  return 1./(1.+pow((v[0]-_mean)/_sigma, 2))/(M_PI*_sigma);
}

double DCauchy::Level(const double alpha) const {

  Assert::IsProbability("DCauchy::Level", "alpha", alpha);
  return pow(cos(M_PI*alpha/2), 2)/(M_PI*_sigma);
}

double DCauchy::Norm() const {

  return (atan((_upper[0]-_mean)/_sigma)-atan((_lower[0]-_mean)/_sigma))/M_PI;
}

double DCauchy::Radial(const double R) const {

  Assert::IsGreater("DCauchy::Radial", "Radius", R, 0.);
  return 1./(1.+pow(R, 2))/(M_PI*_sigma);
}

double DCauchy::Random() {

  // Randomize the sign and simply return the position for a given p-value
  int sign = (_rdmzer.Rndm() > .5) ? 1 : -1;
  return _mean+sign*_sigma*tan(M_PI*_rdmzer.Rndm()/2);
}

std::vector<double> DCauchy::RandomVector() {

  // Generate a 1-vector of a Cauchy random variable
  return {Random()};
}

bool DCauchy::Initialize() {

  // Assert that sigma is strictly greater than 0
  Assert::IsGreater("DCauchy::Initialize", "sigma", _sigma, 0., true);

  // Get the dimension of the space, always 1
  _dim = 1;

  // Initialize the default name of the function and its title
  _name = "cauchy";
  _title = "Cauchy";

  // Set up a default range around the mean (+/-10 widths)
  _lower.push_back(_mean-10*_sigma);
  _upper.push_back(_mean+10*_sigma);

  // Initialize the pseudorandom number generator
  time_t tt;
  _rdmzer = TRandom3(time(&tt));

  return true;
}
