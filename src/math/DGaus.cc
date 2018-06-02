#include "DGaus.hh"

// ROOT style function necessary to numerically compute 2 and 3D integrals (normalisation)
double RGaus(double *x, double *par) {

  // Clear the easy case of n = 1
  if ( par[0] == 1 )
      return exp(-pow(x[0]-par[2], 2)/(2*par[1]))/sqrt(2*M_PI*par[1]);

  // Extract the parameters, feed them to the function
  Matrix<double> icov(par[0], par[0]), X(par[0], 1), Xt;
  size_t i, j;
  for (i = 0; i < par[0]; i++) {
    X[i][0] = x[i]-par[2+i];
    for (j = 0; j < par[0]; j++)
        icov[i][j] = par[2+(int)par[0]+i*(int)par[0]+j];
  }
  Xt = X.Transpose();

  return exp(-(Xt*icov*X)[0][0]/2)/pow(2*M_PI, par[0]/2)/sqrt(par[1]);
}

// ROOT style function necessary to produce a 3D contour, must be at global level
double RGausOffset(double *x, double *par) {

  return RGaus(x, &par[1])-par[0];
}

DGaus::DGaus(const size_t n) :
  _mean(n, 0), _cov(n, n), _invcov(), _J(), _det(0.) {

  try {
    Assert::IsNonZero("DGaus::DGaus", "Dimension", n);
    _dim = n;
    _cov.Identity();

    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "DGaus::DGaus"));
  }
}

DGaus::DGaus(const double mu,
	     const double sigma) :
  _mean(1, mu), _cov(1, 1), _invcov(), _J(), _det(0.) {

  try {
    Assert::IsGreater("DGaus::DGaus", "sigma", sigma, 0., true);
    _dim = 1;
    _cov[0][0] = sigma*sigma;

    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "DGaus::DGaus"));
  }
}

DGaus::DGaus(const std::vector<double>& mu,
	     const std::vector<double>& sigma) :
  _mean(mu), _cov(sigma.size(), sigma.size()), _invcov(), _J(), _det(0.) {

  try {
    Assert::IsNonZero("DGaus::DGaus", "Dimension", mu.size());
    Assert::SameSize("DGaus::DGaus", "Mean and width vectors", mu, sigma);
    _dim = _mean.size();
    size_t i, j;
    for (i = 0; i < _dim; i++) {
      Assert::IsGreater("DGaus::DGaus", "sigma["+std::to_string(i)+"]", sigma[i], 0., true);
      for (j = 0; j < _dim; j++)
          _cov[i][j] = (i == j) ? sigma[i]*sigma[i] : 0.;
    }

    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "DGaus::DGaus"));
  }
}

DGaus::DGaus(const std::vector<double>& mu,
	     const Matrix<double>& cov) :
  _mean(mu), _cov(cov), _invcov(), _J(), _det(0.) {

  try {
    Assert::IsNonZero("DGaus::DGaus", "Dimension", cov.Nrows());
    Assert::SameSize("DGaus::DGaus", "Mean and covariance matrix", mu, cov[0]);
    Assert::IsSquare("DGaus::DGaus", "Covariance matrix", cov.std());
    _dim = mu.size();

    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "DGaus::DGaus"));
  }
}

DGaus::DGaus(const DGaus& df) {
  *this = df;
}

DGaus& DGaus::operator=(const DGaus& df) {
  if ( this == &df )
      return *this;

  _mean = df._mean;
  _cov = df._cov;
  _invcov = df._invcov;
  _J = df._J;
  _det = df._det;

  return *this;
}

DGaus::~DGaus () {}

double DGaus::CDF(const double& v) const {

  // CDF only defined for 1D
  Assert::IsEqual("DGaus::CDF", "Dimension", _dim, (size_t)1);
  return .5*(1.+erf((v-_mean[0])/(sqrt(2*_det))));
}

double DGaus::CDFRadial(const double& R) const {

  // The Gaussian radial CDF is the CDF of the chi squared distribution
  Assert::IsGreater("DGaus::CDFRadial", "Radius", R, 0.);
  return TMath::Gamma((double)_dim/2, R*R/2);
}

std::vector<TLine*> DGaus::Contour1D(const double alpha) const {

  Assert::IsProbability("DGaus::Contour1D", "alpha", alpha);
  double R = TMath::ErfInverse(alpha);
  double max = exp(-R*R)/sqrt(2*M_PI*_det);
  TLine* line = new TLine(_mean[0]-sqrt(2*_det)*R, max,
			  _mean[0]+sqrt(2*_det)*R, max);
  line->SetLineWidth(2);
  return {line};
}

std::vector<TObject*> DGaus::Contour2D(const double alpha) const {

  // Define the ellipse
  Assert::IsProbability("DGaus::Contour2D", "alpha", alpha);
  TEllipse *ell = new TEllipse();
  ell->SetFillColorAlpha(0, 0);
  ell->SetLineWidth(2);

  // First set the centre of the ellipse at the mean
  ell->SetX1(_mean[0]);
  ell->SetY1(_mean[1]);

  // Given the probability (p-value), find the corresponding value of the chi-squared distribution
  // The cumulative density distribution of a 2 DOF chi^2 distribution is expressed as
  // F(z) = gamma(1, z/2)/Gamma(1) = gamma(1, z/2) = (1-exp(-z/2)),
  // which allows us to yield z = -2*ln(1-p) as the corresponding value.
  double rad2 = -2*log(1.-alpha);

  // Compute the eigenvalues of the covariance matrix as they are related to the axes dimensions
  std::vector<double> lambdas;
  Matrix<double> V = _cov.EigenVectors(lambdas);
  ell->SetR1(sqrt(lambdas[0]*rad2));
  ell->SetR2(sqrt(lambdas[1]*rad2));

  // Compute the eigenvectors to yield the rotation angle of the ellipse
  ell->SetTheta(atan2(V[1][0], V[0][0])*360./(2*M_PI));
  ell->SetPhimin(0);
  ell->SetPhimax(360);

  return {ell};
}

TF3* DGaus::Contour3D(const double alpha) const {

  // Find the level of the isosurface for the provided alpha
  Assert::IsProbability("DGaus::Contour3D", "alpha", alpha);
  double C = Level(alpha);

  // Define and return the offset function
  std::vector<double> pars = OffsetParameters(C);
  TF3* ellipsoid = new TF3("ellipsoid", RGausOffset, _lower[0], _upper[0], _lower[1],
			   _upper[1], _lower[2], _upper[2], pars.size());
  ellipsoid->SetTitle(_title.c_str());
  ellipsoid->SetParameters(&(pars[0]));
  ellipsoid->SetFillColorAlpha(1, .25);
  ellipsoid->SetLineColor(1);
  return ellipsoid;
}

std::vector<TLine*> DGaus::ContourRadial(const double alpha) const {

  Assert::IsProbability("DGaus::ContourRadial", "alpha", alpha);
  double R = Radius(alpha);
  double level = Radial(R);
  TLine* line = new TLine(0, level, R, level);
  line->SetLineWidth(2);
  return {line};
}

double DGaus::ContourVolume(const double alpha) const {

  Assert::IsProbability("DGaus::ContourVolume", "alpha", alpha);

  // Clear the simple cases
  if ( _dim == 1 ) {
    return 2*sqrt(2*_det)*TMath::ErfInverse(alpha);
  } else if ( _dim == 2 ) {
    return M_PI*sqrt(_det)*(-2*log(1-alpha));
  }

  // In higher dimension, simply use the general definition for an n-ellipse
  double R = Radius(alpha);
  return pow(M_PI, (double)_dim/2)*pow(R, _dim)*sqrt(_det)/tgamma((double)_dim/2+1);
}

double DGaus::Evaluate(const std::vector<double>& v) const {

  Assert::SameSize("DGaus::Evaluate", "Function space and argument", _lower, v);

  // Clear the simple case of n = 1
  if ( _dim == 1 )
      return exp(-pow(v[0]-_mean[0], 2)/_det/2)/sqrt(2*M_PI*_det);

  // If higher dimension, need to define vectors
  Matrix<double> X(v), M(_mean), Xt;
  X = X-M;
  Xt = X.Transpose();
  return exp(-(Xt*_invcov*X)[0][0]/2)/pow(2*M_PI, (double)_dim/2)/sqrt(_det);
}

double DGaus::Level(const double alpha) const {

  Assert::IsProbability("DGaus::Level", "alpha", alpha);
  double R = Radius(alpha);
  return Radial(R);
}

double DGaus::Norm() const {

  // Clear the simple case of n = 1
  if ( _dim == 1 )
      return Norm1D(0);

  // For higher dimensions, it is a tricky analytical challenge to integrate a correlated
  // gaussian over a orthotopical interval. One needs to use a numerical method,
  // compute the integral over the interval, correct to have I = 1.
  // For n = 2 or 3, feed it to ROOT
  std::vector<double> pars = Parameters();
  if ( _dim == 2 ) {
    TF2 f("", RGaus, _lower[0], _upper[0], _lower[1], _upper[1], pars.size());
    f.SetParameters(&(pars[0]));
    return f.Integral(_lower[0], _upper[0], _lower[1], _upper[1]);
  } else if ( _dim == 3 ) {
    TF3 f("", RGaus, _lower[0], _upper[0], _lower[1], _upper[1], _lower[2], _upper[2], pars.size());
    f.SetParameters(&(pars[0]));
    return f.Integral(_lower[0], _upper[0], _lower[1], _upper[1], _lower[2], _upper[2]);
  }
  
  // For higher dimensions, no function available in ROOT, only support uncorrelated
  // Check that the matrix is diagonal, if not, throw recoverable
  size_t i, j;
  for (i = 0; i < _dim; i++)
    for (j = 0; j < _dim; j++)
      if ( i != j && _cov[i][j] )
	  throw(Exceptions::Exception(Exceptions::recoverable,
	        "The matrix is correlated, norm not supported",
		"DGaus::Norm"));

  // If uncorrelated, the norm is simply the product of the norms in each dimension
  double norm(1.);
  for (i = 0; i < _dim; i++)
      norm *= Norm1D(i);

  return norm;
}

double DGaus::Radial(const double R) const {

  Assert::IsGreater("DGaus::Radial", "Radius", R, 0.);
  return exp(-R*R/2)/pow(2*M_PI, (double)_dim/2.)/sqrt(_det);
}

double DGaus::Radius(const double alpha) const {

  Assert::IsProbability("DGaus::Radius", "alpha", alpha);

  // Use the algorithm in Algorithm AS 91 Appl. Statist. (1975) Vol.24, P.35
  return sqrt(TMath::ChisquareQuantile(alpha, _dim));
}

double DGaus::Random() {

  return _rdmzer.Gaus(_mean[0], sqrt(_cov[0][0]));
}

std::vector<double> DGaus::RandomVector() {

  // First generate a vector of uncorrelated normal variables
  Vector<double> x(_dim);
  size_t i;
  for (i = 0; i < _dim; i++)
      x[i] = _rdmzer.Gaus(0., 1);

  // Convert the coordinates in the reference of the n-sphere to coordinates in the
  // metric of the covariance matrix using the precomputed Jacobian
  x = _J*x;

  // Offset the vector by the mean and return
  for (i = 0; i < _dim; i++)
      x[i] += _mean[i];
  return x.std();
}

bool DGaus::Initialize() {

  // Initialize the default name of the function and its title
  _name = "gaus"+std::to_string(_dim);
  _title = std::to_string(_dim)+"D Gaussian";

  // Check that the covariance matrix is symmetric
  Assert::IsSymmetric("DGaus::Initialize", "Covariance matrix", _cov.std());

  // Compute the determinant of the covariance matrix, check that it is > 0
  _det = _cov.Determinant();
  Assert::IsGreater("DGaus::Initialize", "Covariance matrix determinant", _det, 0., true);

  // Compute the inverse covariance matrix
  _invcov = _cov.Inverse();

  // Compute the Jacobian
  std::vector<double> lambda(_dim);
  Matrix<double> U = _cov.EigenVectors(lambda);
  Matrix<double> Lsqrt(_dim, _dim);
  size_t i;
  for (i = 0; i < _dim; i++)
     Lsqrt[i][i] = sqrt(lambda[i]);

  _J = U*Lsqrt;

  // Set up a default range around the mean (+/-5 sigma)
  for (i = 0; i < _dim; i++)
      _lower.push_back(_mean[i]-5*sqrt(_cov[i][i]));
  for (i = 0; i < _dim; i++)
      _upper.push_back(_mean[i]+5*sqrt(_cov[i][i]));

  // Initialize the pseudorandom number generator
  time_t tt;
  _rdmzer = TRandom3(time(&tt));

  return true;
}

double DGaus::Norm1D(const size_t i) const {

  // Compute the central normal variables
  double Rl = (_lower[i]-_mean[i])/sqrt(_cov[i][i]);
  double Ru = (_upper[i]-_mean[i])/sqrt(_cov[i][i]);

  // Return the normalisation
  return .5*(erf(Ru/sqrt(2)) - erf(Rl/sqrt(2)));
}

std::vector<double> DGaus::Parameters() const {

  std::vector<double> pars;
  pars.push_back(_dim);

  pars.push_back(_det);

  size_t i, j;
  for (i = 0; i < _dim; i++)
      pars.push_back(_mean[i]);

  for (i = 0; i < _dim; i++)
    for (j = 0; j < _dim; j++)
	pars.push_back(_invcov[i][j]);

  return pars;
}

std::vector<double> DGaus::OffsetParameters(const double C) const {

  std::vector<double> pars = Parameters();
  pars.insert(pars.begin(), C);
  return pars;
}
