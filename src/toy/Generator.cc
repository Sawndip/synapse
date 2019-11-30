#include "Generator.hh"

Generator::Generator() :
  _dim(), _mass(), _p(), _S() {
}

Generator::Generator(const size_t dim,
		     const double& mass) :
  _dim(dim), _mass(mass), _p(), _S() {

  try {
    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "Generator::Generator"));
  }
}

Generator::Generator(const Generator& gen) {
  *this = gen;
}

Generator& Generator::operator=(const Generator& gen) {
  if ( this == &gen )
      return *this;

  _dim = gen._dim;
  _mass = gen._mass;
  _p = gen._p;
  _S = gen._S;

  return *this;
}

Generator::~Generator () {}

void Generator::Initialize() {
  // Check that the dimension is acceptable
  if ( _dim != 2 && _dim != 4 )
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	    "Transverse phase space dimension not supported: "+std::to_string(_dim),
	    "Generator::Initialize"));

  // Check that the mass is positive
  if ( _mass <= 0. )
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	    "Particle mass is non physical: "+std::to_string(_mass),
	    "Generator::Initialize"));
}

void Generator::SetMatrixParametrisation(const double& p,
					 const double& eps,
		      	  		 const double& beta,
		      	  		 const double alpha) {

  // Initialize covariance, set the central momentum
  _S.Resize(_dim, _dim);
  _S.Identity();
  _p = p;

  // Generate a 2D correlation matrix
  if ( !alpha ) {
    // If alpha is 0., the beam is not correlated
    double sig2 = _mass*eps*beta/p;
    _S[0][0] = sig2;
    _S[1][1] = pow(_mass*eps, 2)/sig2;

  } else {
    // If alpha is non zero, produce correlated accordingly
    double gamma = (1.+alpha*alpha)/beta;
    _S[0][0] = beta*_mass*eps/p;
    _S[0][1] = -alpha*_mass*eps;
    _S[1][0] = _S[0][1];
    _S[1][1] = gamma*_mass*eps*p;
  }

  // If it is a 4D transverse phase space, add a block
  if ( _dim > 2 )
    for (size_t i = 0; i < 2; i++)
      for (size_t j = 0; j < 2; j++)
          _S[i+2][j+2] = _S[i][j];
}

Beam::BunchMap Generator::GaussianBunchMap(const size_t n) const {

  // Sample repetedly from a Gaussian
  Beam::BunchMap samples;
  std::vector<double> sample(_dim);

  std::vector<double> mean(_dim, 0.);
  DGaus gaus(mean, _S);
  for (size_t i = 0; i < n; i++) {
    sample = gaus.RandomVector();
    samples["x"].push_back(sample[0]);
    samples["px"].push_back(sample[1]);
    if ( _dim > 2 ) {
      samples["y"].push_back(sample[2]);
      samples["py"].push_back(sample[3]);
    }
    samples["pz"].push_back(_p);
  }

  return samples;
}

Beam::Bunch Generator::GaussianBunch(const size_t n) const {

  return Beam::Bunch(GaussianBunchMap(n), 0.);
}

Beam::BunchMap Generator::SpiralBunchMap(const size_t n,
					 const double& t) const {

  // Sample repetedly from a Gaussian
  Beam::BunchMap samples;
  std::vector<double> sample(_dim);

  std::vector<double> mean(_dim, 0.);
  DSpiral gaus(mean, _S, t);
  for (size_t i = 0; i < n; i++) {
    sample = gaus.RandomVector();
    samples["x"].push_back(sample[0]);
    samples["px"].push_back(sample[1]);
    if ( _dim > 2 ) {
      samples["y"].push_back(sample[2]);
      samples["py"].push_back(sample[3]);
    }
    samples["pz"].push_back(_p);
  }

  return samples;
}

Beam::Bunch Generator::SpiralBunch(const size_t n,
				   const double& t) const {

  return Beam::Bunch(SpiralBunchMap(n, t), 0.);
}
