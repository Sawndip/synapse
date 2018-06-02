#include "PBATDE.hh"

PBATDE::PBATDE() :
  _dim(0), _J(0), _points(0), _tdes(0), _interp(false) {
}

PBATDE::PBATDE(const std::vector<std::vector<double>>& points,
	       const bool interp,
	       const size_t J) :
  _dim(0), _J(J), _points(points), _tdes(0), _interp(interp) {

  try {
    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "PBATDE::PBATDE"));
  }
}

PBATDE::PBATDE(const PBATDE& pbatde) {
  *this = pbatde;
}

PBATDE& PBATDE::operator=(const PBATDE& pbatde) {
  if ( this == &pbatde )
      return *this;

  _dim = pbatde._dim;
  _M = pbatde._M;
  _J = pbatde._J;
  _points = pbatde._points;
  _tdes = pbatde._tdes;
  _interp = pbatde._interp;

  return *this;
}

PBATDE::~PBATDE () {}

void PBATDE::Initialize() {

  // Check the input points, extract dim
  Assert::IsNotEmpty("PBATDE::Initialize", "Vector of input points", _points);
  Assert::IsNonZero("PBATDE::Initialize", "Dimension", _points[0].size());
  _dim = _points[0].size();

  // J is the number of samples in each of the bags. Use random resampling with
  // replacement to produce M samples of J points from a total of K points
  _M = _points.size()/_J;

  // Initialize one TDE per subsample
  std::vector<std::vector<double>> subsample(_J);
  size_t i, j;
  for (i = 0; i < _M; i++) {
    // Extract a subsample of J points from the total
//    subsample = Math::Resample(_points, _J, false);
    for (j = 0; j < _J; j++)
	subsample[j] = _points[i*_J+j];
    _tdes.push_back(TDE(subsample, _interp, false)); // No need for boundaries
  }
}

double PBATDE::Evaluate(const std::vector<double>& v) const {

  // Take the mean of the estimation of each TDE
  double value(0.);
  size_t i;
  for (i = 0; i < _M; i++)
      value += _tdes[i](v);

  return value/_M;

/*  size_t M = _M;
  double temp;
  for (i = 0; i < _M; i++) {
    temp = _tdes[i](v);
    if ( !temp ) {
      M -= 1;
      continue;
    }

    value += 1./temp;
  }

  if ( M )
      return M/value; // Yes it makes sense FD, do the math

  return 0.;*/
}

double PBATDE::Evaluate(const double& v) const {

  return Evaluate(std::vector<double>({v}));
}

double PBATDE::Evaluate(const double* v) const {

  std::vector<double> vv(v, v+_dim);
  return Evaluate(vv);
}
