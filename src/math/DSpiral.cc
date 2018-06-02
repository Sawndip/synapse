#include "DSpiral.hh"

// ROOT style function necessary to numerically compute 2 and 3D integrals (normalisation)
double RSpiral(double *x, double *par) {

  // Clear the easy case of n = 1
  if ( par[0] == 1 )
      return exp(-pow(x[0]-par[2], 2)/(2*par[1]))/sqrt(2*M_PI*par[1]);

  // Extract the parameters
  Matrix<double> icov(par[0], par[0]), X(par[0], 1), Xt;
  size_t i, j;
  for (i = 0; i < par[0]; i++) {
    X[i][0] = x[i]-par[3+i];
    for (j = 0; j < par[0]; j++)
        icov[i][j] = par[3+(int)par[0]+i*(int)par[0]+j];
  }
  Xt = X.Transpose();

  // Find the squared radius of the point
  double R2 = (Xt*X)[0][0];

  // Define a 2D rotation matrix based on the radius
  double theta = par[2]*R2;
  Matrix<double> P(2, 2);
  P[0][0] = cos(theta);
  P[0][1] = -sin(theta);
  P[1][0] = -P[0][1];
  P[1][1] = P[0][0];

  // Rotate about each axis individually (d*(d-1)/2 dxd matrices)
  Matrix<double> T(par[0], par[0]), A(par[0], par[0]);
  T.Identity();
  for (size_t i = 0; i < par[0]; i++) {
    for (size_t j = i+1; j < par[0]; j++) {
      A.Identity();
      A[i][i] = P[0][0];
      A[i][j] = pow(-1, i+j+1)*P[0][1];
      A[j][i] = pow(-1, i+j+1)*P[1][0];
      A[j][j] = P[1][1];
      T *= A;
    }
  }
  Matrix<double> Tt = T.Transpose();

  // Apply the rotation matrix to the inverse covariance matrix
  R2 = (Xt*(T*icov*Tt)*X)[0][0];

  return exp(-R2/2)/pow(2*M_PI, par[0]/2)/sqrt(par[1]);
}

// ROOT style function necessary to produce a 3D contour, must be at global level
double RSpiralOffset(double *x, double *par) {

  return RSpiral(x, &par[1])-par[0];
}

DSpiral::DSpiral(const size_t n,
		 const double turn) :
  _mean(n, 0), _cov(n, n), _invcov(), _J(), _det(0.), _turn(turn) {

  try {
    Assert::IsNonZero("DSpiral::DSpiral", "Dimension", n);
    _dim = n;
    _cov.Identity();

    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "DSpiral::DSpiral"));
  }
}

DSpiral::DSpiral(const double mu,
	     	 const double sigma,
		 const double turn) :
  _mean(1, mu), _cov(1, 1), _invcov(), _J(), _det(0.), _turn(turn) {

  try {
    Assert::IsGreater("DSpiral::DSpiral", "sigma", sigma, 0., true);
    _dim = 1;
    _cov[0][0] = sigma*sigma;

    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "DSpiral::DSpiral"));
  }
}

DSpiral::DSpiral(const std::vector<double>& mu,
	     	 const std::vector<double>& sigma,
		 const double turn) :
  _mean(mu), _cov(sigma.size(), sigma.size()), _invcov(), _J(), _det(0.), _turn(turn) {

  try {
    Assert::IsNonZero("DSpiral::DSpiral", "Dimension", mu.size());
    Assert::SameSize("DSpiral::DSpiral", "Mean and width vectors", mu, sigma);
    _dim = _mean.size();
    size_t i, j;
    for (i = 0; i < _dim; i++) {
      Assert::IsGreater("DSpiral::DSpiral", "sigma["+std::to_string(i)+"]", sigma[i], 0., true);
      for (j = 0; j < _dim; j++)
          _cov[i][j] = (i == j) ? sigma[i]*sigma[i] : 0.;
    }

    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "DSpiral::DSpiral"));
  }
}

DSpiral::DSpiral(const std::vector<double>& mu,
	     	 const Matrix<double>& cov,
		 const double turn) :
  _mean(mu), _cov(cov), _invcov(), _J(), _det(0.), _turn(turn) {

  try {
    Assert::IsNonZero("DSpiral::DSpiral", "Dimension", cov.Nrows());
    Assert::SameSize("DSpiral::DSpiral", "Mean and covariance matrix", mu, cov[0]);
    Assert::IsSquare("DSpiral::DSpiral", "Covariance matrix", cov.std());
    _dim = mu.size();

    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "DSpiral::DSpiral"));
  }
}

DSpiral::DSpiral(const DSpiral& df) {
  *this = df;
}

DSpiral& DSpiral::operator=(const DSpiral& df) {
  if ( this == &df )
      return *this;

  _mean = df._mean;
  _cov = df._cov;
  _invcov = df._invcov;
  _J = df._J;
  _det = df._det;
  _turn = df._turn;

  return *this;
}

DSpiral::~DSpiral () {}

double DSpiral::CDF(const double& v) const {

  // CDF only defined for 1D
  Assert::IsEqual("DSpiral::CDF", "Dimension", _dim, (size_t)1);
  return .5*(1.+erf((v-_mean[0])/(sqrt(2*_det))));
}

double DSpiral::CDFRadial(const double& R) const {

  // The Spiral radial CDF is the CDF of the chi squared distribution
  Assert::IsGreater("DSpiral::CDFRadial", "Radius", R, 0.);
  return TMath::Gamma((double)_dim/2, R*R/2);
}

std::vector<TLine*> DSpiral::Contour1D(const double alpha) const {

  Assert::IsProbability("DSpiral::Contour1D", "alpha", alpha);
  double R = TMath::ErfInverse(alpha);
  double max = exp(-R*R)/sqrt(2*M_PI*_det);
  TLine* line = new TLine(_mean[0]-sqrt(2*_det)*R, max,
			  _mean[0]+sqrt(2*_det)*R, max);
  line->SetLineWidth(2);
  return {line};
}

std::vector<TObject*> DSpiral::Contour2D(const double alpha) const {

  // Initialize the contour
  double level = Level(alpha);
  std::vector<double> pars = Parameters();
  TF2* contour = new TF2("contour", RSpiral,
	_lower[0], _upper[0], _lower[1], _upper[1], pars.size());
  contour->SetNpx(250);
  contour->SetNpy(250);
  contour->SetParameters(&pars[0]);
  contour->SetContour(1);
  contour->SetContourLevel(0, level);
  contour->SetLineWidth(2);
  contour->SetLineColor(kBlack);
  return {contour};
}

TF3* DSpiral::Contour3D(const double alpha) const {

  // Find the level of the isosurface for the provided alpha
  Assert::IsProbability("DSpiral::Contour3D", "alpha", alpha);
  double C = Level(alpha);

  // Define and return the offset function
  std::vector<double> pars = OffsetParameters(C);
  TF3* spiraloid = new TF3("ellipsoid", RSpiralOffset, _lower[0], _upper[0], _lower[1],
			   _upper[1], _lower[2], _upper[2], pars.size());
  spiraloid->SetTitle(_title.c_str());
  spiraloid->SetParameters(&(pars[0]));
  spiraloid->SetFillColorAlpha(1, .25);
  spiraloid->SetLineColor(1);
  return spiraloid;
}

std::vector<TLine*> DSpiral::ContourRadial(const double alpha) const {

  Assert::IsProbability("DSpiral::ContourRadial", "alpha", alpha);
  double R = Radius(alpha);
  double level = Radial(R);
  TLine* line = new TLine(0, level, R, level);
  line->SetLineWidth(2);
  return {line};
}

double DSpiral::ContourVolume(const double alpha) const {

  Assert::IsProbability("DSpiral::ContourVolume", "alpha", alpha);

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

double DSpiral::Evaluate(const std::vector<double>& v) const {

  Assert::SameSize("DSpiral::Evaluate", "Function space and argument", _lower, v);

  // Clear the simple case of n = 1
  if ( _dim == 1 )
      return exp(-pow(v[0]-_mean[0], 2)/_det/2)/sqrt(2*M_PI*_det);

  // Find the squared radius of the point, get the twisting matrix and its transpose
  Matrix<double> X(v), M(_mean), Xt;
  X -= M;
  Xt = X.Transpose();
  double R2 = (Xt*X)[0][0];
  Matrix<double> T = TwistingMatrix(R2);
  Matrix<double> Tt = T.Transpose();

  // Apply the rotation matrix to the inverse covariance matrix
  R2 = (Xt*(T*_invcov*Tt)*X)[0][0];

  // If higher dimension, need to define vectors
  return exp(-R2/2)/pow(2*M_PI, (double)_dim/2)/sqrt(_det);
}

double DSpiral::Level(const double alpha) const {

  Assert::IsProbability("DSpiral::Level", "alpha", alpha);
  double R = Radius(alpha);
  return Radial(R);
}

double DSpiral::Norm() const {

  // Clear the simple case of n = 1
  if ( _dim == 1 )
      return Norm1D(0);

  // For higher dimensions, it is a tricky analytical challenge to integrate a correlated
  // spiral over a orthotopical interval. One needs to use a numerical method,
  // compute the integral over the interval, correct to have I = 1.
  // For n = 2 or 3, feed it to ROOT
  std::vector<double> pars = Parameters();
  if ( _dim == 2 ) {
    TF2 f("", RSpiral, _lower[0], _upper[0], _lower[1], _upper[1], pars.size());
    f.SetParameters(&(pars[0]));
    return f.Integral(_lower[0], _upper[0], _lower[1], _upper[1]);
  } else if ( _dim == 3 ) {
    TF3 f("", RSpiral, _lower[0], _upper[0], _lower[1], _upper[1], _lower[2], _upper[2], pars.size());
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
		"DSpiral::Norm"));

  // If uncorrelated, the norm is simply the product of the norms in each dimension
  double norm(1.);
  for (i = 0; i < _dim; i++)
      norm *= Norm1D(i);

  return norm;
}

double DSpiral::Radial(const double R) const {

  Assert::IsGreater("DSpiral::Radial", "Radius", R, 0.);
  return exp(-R*R/2)/pow(2*M_PI, (double)_dim/2.)/sqrt(_det);
}

double DSpiral::Radius(const double alpha) const {

  Assert::IsProbability("DSpiral::Radius", "alpha", alpha);

  // Use the algorithm in Algorithm AS 91 Appl. Statist. (1975) Vol.24, P.35
  return sqrt(TMath::ChisquareQuantile(alpha, _dim));
}

double DSpiral::Random() {

  return _rdmzer.Gaus(_mean[0], sqrt(_cov[0][0]));
}

std::vector<double> DSpiral::RandomVector() {

  // First generate a vector of uncorrelated normal variables
  Vector<double> x(_dim);
  size_t i;
  for (i = 0; i < _dim; i++)
      x[i] = _rdmzer.Gaus(0., 1);

  // Convert the coordinates in the reference of the n-sphere to coordinates in the
  // metric of the covariance matrix using the precomputed Jacobian
  x = _J*x;

  // Apply the twisting matrix to the point
  double R2 = x.dot(x);
  Matrix<double> T = TwistingMatrix(R2);
  x = T*x;

  // Offset the vector by the mean and return
  for (i = 0; i < _dim; i++)
      x[i] += _mean[i];
  return x.std();
}

bool DSpiral::Initialize() {

  // Initialize the default name of the function and its title
  _name = "spiral"+std::to_string(_dim);
  _title = std::to_string(_dim)+"D Spiral";

  // Check that the covariance matrix is symmetric
  Assert::IsSymmetric("DSpiral::Initialize", "Covariance matrix", _cov.std());

  // Compute the determinant of the covariance matrix, check that it is > 0
  _det = _cov.Determinant();
  Assert::IsGreater("DSpiral::Initialize", "Covariance matrix determinant", _det, 0., true);

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

double DSpiral::Norm1D(const size_t i) const {

  // Compute the central normal variables
  double Rl = (_lower[i]-_mean[i])/sqrt(_cov[i][i]);
  double Ru = (_upper[i]-_mean[i])/sqrt(_cov[i][i]);

  // Return the normalisation
  return .5*(erf(Ru/sqrt(2)) - erf(Rl/sqrt(2)));
}

std::vector<double> DSpiral::Parameters() const {

  std::vector<double> pars;
  pars.push_back(_dim);

  pars.push_back(_det);

  pars.push_back(_turn);

  size_t i, j;
  for (i = 0; i < _dim; i++)
      pars.push_back(_mean[i]);

  for (i = 0; i < _dim; i++)
    for (j = 0; j < _dim; j++)
	pars.push_back(_invcov[i][j]);

  return pars;
}

std::vector<double> DSpiral::OffsetParameters(const double C) const {

  std::vector<double> pars = Parameters();
  pars.insert(pars.begin(), C);
  return pars;
}

Matrix<double> DSpiral::TwistingMatrix(const double R2) const {

  // Define a 2D rotation matrix based on the radius
  double theta = _turn*R2;
  Matrix<double> P(2, 2);
  P[0][0] = cos(theta);
  P[0][1] = -sin(theta);
  P[1][0] = -P[0][1];
  P[1][1] = P[0][0];

  // Rotate about each axis individually (d*(d-1)/2 dxd matrices)
  Matrix<double> T(_dim, _dim), A(_dim, _dim);
  T.Identity();
  for (size_t i = 0; i < _dim; i++) {
    for (size_t j = i+1; j < _dim; j++) {
      A.Identity();
      A[i][i] = P[0][0];
      A[i][j] = pow(-1, i+j+1)*P[0][1];
      A[j][i] = pow(-1, i+j+1)*P[1][0];
      A[j][j] = P[1][1];
      T *= A;
    }
  }

  return T;
}
