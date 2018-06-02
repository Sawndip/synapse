#include "DMaxwell.hh"

// ROOT style function necessary to generate random numbers, it is M(r,:phi:)drd^(n-1)phi
double RMaxwellRadialDist(double *x, double *par) {

  return pow(x[0], par[0]+1)*exp(-x[0]*x[0]/2);
}

// ROOT style function necessary to numerically compute the optimal radii
double RMaxwellRadial(double *x, double *par) {

  return pow(x[0], 2)*exp(-x[0]*x[0]/2)/par[0]/pow(2*M_PI, par[0]/2.);
}

// ROOT style function necessary to numerically compute the radii for a given level
double RMaxwellPContent(double *x, double *par) {

  // Cannot compute analytically. Compute the radii by finding the right level
  TF1 f("", RMaxwellRadial, 0, 100, 1);
  f.SetParameter(0, par[0]);
  par[1] = f.GetX(x[0], 0, sqrt(2));
  par[2] = f.GetX(x[0], sqrt(2), 100);

  return TMath::Gamma(par[0]/2+1, par[2]*par[2]/2)-TMath::Gamma(par[0]/2+1, par[1]*par[1]/2);
}

// ROOT style function necessary to numerically compute 2 and 3D integrals (normalisation)
double RMaxwell(double *x, double *par) {

  // Clear the easy case of n = 1
  if ( par[0] == 1 )
      return pow(x[0]-par[2], 2)*exp(-pow(x[0]-par[2], 2)/(2*par[1]))/sqrt(2*M_PI*pow(par[1], 3));

  // Extract the parameters, feed them to the function
  Matrix<double> icov(par[0], par[0]), X(par[0], 1), Xt;
  size_t i, j;
  for (i = 0; i < par[0]; i++) {
    X[i][0] = x[i]-par[2+i];
    for (j = 0; j < par[0]; j++)
        icov[i][j] = par[2+(int)par[0]+i*(int)par[0]+j];
  }
  Xt = X.Transpose();
  double R2 = (Xt*icov*X)[0][0];

  return R2*exp(-R2/2)/par[0]/pow(2*M_PI, par[0]/2)/pow(par[1], 3./2);
}

// ROOT style function necessary to produce a 3D contour, must be at global level
double RMaxwellOffset(double *x, double *par) {

  return RMaxwell(x, &par[1])-par[0];
}

DMaxwell::DMaxwell(const size_t n) :
  _mean(n, 0), _cov(n, n), _invcov(), _J(), _det(0.), _frand() {

  try {
    Assert::IsNonZero("DMaxwell::DMaxwell", "Dimension", n);
    _dim = n;
    _cov.Identity();

    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "DMaxwell::DMaxwell"));
  }
}

DMaxwell::DMaxwell(const double mu,
	     	   const double sigma) :
  _mean(1, mu), _cov(1, 1), _invcov(), _J(), _det(0.), _frand() {

  try {
    Assert::IsGreater("DMaxwell::DMaxwell", "sigma", sigma, 0., true);
    _dim = 1;
    _cov[0][0] = sigma*sigma;

    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "DMaxwell::DMaxwell"));
  }
}

DMaxwell::DMaxwell(const std::vector<double>& mu,
	     	   const std::vector<double>& sigma) :
  _mean(mu), _cov(sigma.size(), sigma.size()), _invcov(), _J(), _det(0.), _frand() {

  try {
    Assert::IsNonZero("DMaxwell::DMaxwell", "Dimension", mu.size());
    Assert::SameSize("DMaxwell::DMaxwell", "Mean and width vectors", mu, sigma);
    _dim = _mean.size();
    size_t i, j;
    for (i = 0; i < _dim; i++) {
      Assert::IsGreater("DMaxwell::DMaxwell", "sigma["+std::to_string(i)+"]", sigma[i], 0., true);
      for (j = 0; j < _dim; j++)
          _cov[i][j] = (i == j) ? sigma[i]*sigma[i] : 0.;
    }

    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "DMaxwell::DMaxwell"));
  }
}

DMaxwell::DMaxwell(const std::vector<double>& mu,
	     	   const Matrix<double>& cov) :
  _mean(mu), _cov(cov), _invcov(), _J(), _det(0.), _frand() {

  try {
    Assert::IsNonZero("DMaxwell::DMaxwell", "Dimension", cov.Nrows());
    Assert::SameSize("DMaxwell::DMaxwell", "Mean and covariance matrix", mu, cov[0]);
    Assert::IsSquare("DMaxwell::DMaxwell", "Covariance matrix", cov.std());
    _dim = mu.size();

    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "DMaxwell::DMaxwell"));
  }
}

DMaxwell::DMaxwell(const DMaxwell& df) {
  *this = df;
}

DMaxwell& DMaxwell::operator=(const DMaxwell& df) {
  if ( this == &df )
      return *this;

  _mean = df._mean;
  _cov = df._cov;
  _invcov = df._invcov;
  _J = df._J;
  _det = df._det;
  _frand = df._frand;

  return *this;
}

DMaxwell::~DMaxwell () {}

double DMaxwell::CDF(const double& v) const {

  // CDF only defined for 1D
  Assert::IsEqual("DMaxwell::CDF", "Dimension", _dim, (size_t)1);
  double z = (v-_mean[0])/sqrt(2*_det);
  return .5*(1+erf(z)) - z*exp(-z*z)/sqrt(M_PI);
}

double DMaxwell::CDFRadial(const double& R) const {

  // The Maxwell radial CDF follows the lower incomplete Gamma function
  Assert::IsGreater("DMaxwell::CDFRadial", "Radius", R, 0.);
  return TMath::Gamma((double)_dim/2+1., R*R/2);
}

std::vector<TLine*> DMaxwell::Contour1D(const double alpha) const {

  Assert::IsProbability("DMaxwell::Contour1D", "alpha", alpha);
  double r, R, l;
  Radii(alpha, r, R, l);
  TLine* linel = new TLine(_mean[0]-R, l, _mean[0]-r, l);
  linel->SetLineWidth(2);
  TLine* liner = new TLine(_mean[0]+r, l, _mean[0]+R, l);
  liner->SetLineWidth(2);
  return {linel, liner};
}

std::vector<TObject*> DMaxwell::Contour2D(const double alpha) const {

  // Define two concentric ellipses
  Assert::IsProbability("DMaxwell::Contour2D", "alpha", alpha);
  std::vector<TEllipse*> ells = {new TEllipse(), new TEllipse()};
  for (TEllipse *ell : ells) {
    ell->SetFillColorAlpha(0, 0);
    ell->SetLineWidth(2);
  }

  // First set the centre of the ellipses at the mean
  for (TEllipse *ell : ells) {
    ell->SetX1(_mean[0]);
    ell->SetY1(_mean[1]);
  }

  // Given the probability (p-value), find the corresponding radii
  double r, R, l;
  Radii(alpha, r, R, l);

  // Compute the eigenvalues of the covariance matrix as they are related to the axes dimensions
  std::vector<double> lambdas;
  Matrix<double> V = _cov.EigenVectors(lambdas);
  ells[0]->SetR1(sqrt(lambdas[0])*r);
  ells[0]->SetR2(sqrt(lambdas[1])*r);
  ells[1]->SetR1(sqrt(lambdas[0])*R);
  ells[1]->SetR2(sqrt(lambdas[1])*R);

  // Compute the eigenvectors to yield the rotation angle of the ellipse
  for (TEllipse *ell : ells) {
    ell->SetTheta(atan2(V[1][0], V[0][0])*360./(2*M_PI));
    ell->SetPhimin(0);
    ell->SetPhimax(360);
  }

  return {ells[0], ells[1]};
}

TF3* DMaxwell::Contour3D(const double alpha) const {

  // Find the level of the isosurface for the provided alpha
  Assert::IsProbability("DMaxwell::Contour3D", "alpha", alpha);
  double C = Level(alpha);

  // Define and return the offset function
  std::vector<double> pars = OffsetParameters(C);
  TF3* shell = new TF3("shell", RMaxwellOffset, _lower[0], _upper[0], _lower[1],
			   _upper[1], _lower[2], _upper[2], pars.size());
  shell->SetTitle(_title.c_str());
  shell->SetParameters(&(pars[0]));
  shell->SetFillColorAlpha(1, .25);
  shell->SetLineColor(1);
  return shell;
}

std::vector<TLine*> DMaxwell::ContourRadial(const double alpha) const {

  Assert::IsProbability("DMaxwell::ContourRadial", "alpha", alpha);
  double r, R, l;
  Radii(alpha, r, R, l);
  TLine* line = new TLine(r, l, R, l);
  line->SetLineWidth(2);
  return {line};
}

double DMaxwell::ContourVolume(const double alpha) const {

  Assert::IsProbability("DMaxwell::ContourVolume", "alpha", alpha);

  // Cannot compute analytically. Compute the radii and by optimizing
  double r, R, l;
  Radii(alpha, r, R, l);

  // Return the difference between the two n-spheres
  if ( _dim == 1 )
      return 2*(R-r)*sqrt(_det);

  return pow(M_PI, (double)_dim/2.)*(pow(R, _dim)-pow(r, _dim))*sqrt(_det)
	 /tgamma((double)_dim/2.+1);
}

double DMaxwell::Evaluate(const std::vector<double>& v) const {

  Assert::SameSize("DMaxwell::Evaluate", "Function space and argument", _lower, v);

  // Clear the simple case of n = 1
  if ( _dim == 1 )
      return pow(v[0]-_mean[0], 2)*exp(-pow(v[0]-_mean[0], 2)/_det/2)/sqrt(2*M_PI*pow(_det, 3));

  // If higher dimension, need to define vectors
  Matrix<double> X(v), M(_mean), Xt;
  X = X-M;
  Xt = X.Transpose();
  double R2 = (Xt*_invcov*X)[0][0];
  return R2*exp(-R2/2)/pow(2*M_PI, (double)_dim/2)/sqrt(pow(_det, 3))/_dim;
}

double DMaxwell::Level(const double alpha) const {

  Assert::IsProbability("DMaxwell::Level", "alpha", alpha);
  double r, R, l;
  Radii(alpha, r, R, l);
  return l;
}

double DMaxwell::Norm() const {

  // Clear the simple case of n = 1
  if ( _dim == 1 )
      return Norm1D(0);

  // For higher dimensions, it is a tricky analytical challenge to integrate a correlated
  // Maxwell over a orthotopical interval. One needs to use a numerical method,
  // compute the integral over the interval, correct to have I = 1.

  // If n = 2 or 3, feed it to ROOT
  std::vector<double> pars = Parameters();
  if ( _dim == 2 ) {
    TF2 f("", RMaxwell, _lower[0], _upper[0], _lower[1], _upper[1], pars.size());
    f.SetParameters(&(pars[0]));
    return f.Integral(_lower[0], _upper[0], _lower[1], _upper[1]);
  } else if ( _dim == 3 ) {
    TF3 f("", RMaxwell, _lower[0], _upper[0], _lower[1], _upper[1], _lower[2], _upper[2], pars.size());
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
		"DMaxwell::Norm"));

  // If uncorrelated, the norm is simply the product of the norms in each dimension
  double norm(1.);
  for (i = 0; i < _dim; i++)
      norm *= Norm1D(i);

  return norm;
}

double DMaxwell::Radial(const double R) const {

  Assert::IsGreater("DMaxwell::Radial", "Radius", R, 0.);
  return R*R*exp(-R*R/2)/_dim/pow(2*M_PI, (double)_dim/2.)/sqrt(_det);
}

void DMaxwell::Radii(const double alpha, double& r, double& R, double& l) const {

  // Cannot compute analytically. Compute the radii by optimizing the p-value
  Assert::IsProbability("DMaxwell::Radii", "alpha", alpha);
  std::vector<double> pars = {(double)_dim, r, R};
  double max = 2./_dim/pow(2*M_PI, (double)_dim/2)/exp(1);
  TF1 fcont("", RMaxwellPContent, 0, max, 3);
  fcont.SetParameters(&(pars[0]));
  l = fcont.GetX(alpha, 0, max);
  r = fcont.GetParameter(1);
  R = fcont.GetParameter(2);
}

double DMaxwell::Random() {

  // Faster to use ROOT to sample from it than any other technique, analytically challenging
  double sign = (_rdmzer.Rndm() > .5) ? 1 : -1;
  return _mean[0] + _J[0][0]*sign*_frand.GetRandom();
}

std::vector<double> DMaxwell::RandomVector() {

  // For n < 3, simple methods exist
  _frand.GetRandom();	// Strange bug with ROOT...
  double R = _frand.GetRandom();
  Vector<double> x(_dim, 1);
  if ( _dim == 1 ) {
    // Simply randomize the sign
    double sign = (_rdmzer.Rndm() > .5) ? 1 : -1;
    return {_mean[0] + _J[0][0]*sign*R};

  } else if ( _dim == 2 ) {
    // Pick an angle uniformly and compute the corresponding position on a circle
    double th = 2*M_PI*(double)rand()/RAND_MAX;
    x[0] = R*cos(th);
    x[1] = R*sin(th);

  } else if ( _dim == 3 ) {
    // ROOT's method for sphere point picking
    _rdmzer.Sphere(x[0], x[1], x[2], R);

  } else {
    // For higher dimensions, simply pick n points from an n-gaussian (L2 symmetric)
    size_t i;
    for (i = 0; i < _dim; i++)
	x[i] = _rdmzer.Gaus();

    double norm = x.norm();
    for (i = 0; i < _dim; i++)
	x[i] *= R/norm;    
  }

  // Convert the coordinates in the reference of the n-sphere to coordinates in the
  // metric of the covariance matrix using the precomputed the Jacobian
  x = _J*x;

  // Offset the vector by the mean and return
  size_t i;
  for (i = 0; i < _dim; i++)
      x[i] += _mean[i];
  return x.std();
}

bool DMaxwell::Initialize() {

  // Initialize the default name of the function and its title
  _name = "max"+std::to_string(_dim);
  _title = std::to_string(_dim)+"D Maxwell";

  // Check that the covariance matrix is symmetric
  Assert::IsSymmetric("DMaxwell::Initialize", "Covariance matrix", _cov.std());

  // Compute the determinant of the covariance matrix, check that it is > 0
  _det = _cov.Determinant();
  Assert::IsGreater("DMaxwell::Initialize", "Covariance matrix determinant", _det, 0., true);

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

  // Initialize the pseudorandom number generators
  time_t tt;
  _rdmzer = TRandom3(time(&tt));

  if ( gRandom )
     delete gRandom;
  gRandom = new TRandom3(time(&tt));

  _frand = TF1("", RMaxwellRadialDist, 0., 10., 1);
  _frand.SetParameter(0, _dim);

  return true;
}

double DMaxwell::Norm1D(const size_t i) const {

  // Compute the central normal variables
  double Rl = fabs(_lower[i]-_mean[i])/sqrt(_cov[i][i]);
  double Ru = fabs(_upper[i]-_mean[i])/sqrt(_cov[i][i]);

  // Return the normalisation
  if ( _upper[i] < _mean[i] || _lower[i] > _mean[i]  )
      return .5*fabs(TMath::Gamma(3./2, Ru*Ru/2) - TMath::Gamma(3./2, Rl*Rl/2));

  return .5*(TMath::Gamma(3./2, Rl*Rl/2) + TMath::Gamma(3./2, Ru*Ru/2));
}

std::vector<double> DMaxwell::Parameters() const {

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

std::vector<double> DMaxwell::OffsetParameters(const double C) const {

  std::vector<double> pars = Parameters();
  pars.insert(pars.begin(), C);
  return pars;
}
