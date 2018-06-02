#include "DMultiGaus.hh"

// ROOT style function necessary to numerically compute 2 and 3D integrals (normalisation)
double RMultiGaus(double *x, double *par) {

  // Each of the peak a fraction 1/k of the total probability, loop and increment
  Matrix<double> icov(par[0], par[0]), X(par[0], 1), Xt;
  double val = 0.;
  size_t off;
  size_t i, j, k;
  for (i = 0; i < par[1]; i++) {

    // Set the counter offset to use for the current peak
    off = 2+i*(1+par[0]+par[0]*par[0]);
    
    // Clear the easy case of n = 1
    if ( par[0] == 1 ) {
      val += exp(-pow(x[0]-par[off+1], 2)/(2*par[off]))/sqrt(2*M_PI*par[off])/par[1];
      continue;
    }

    // Extract the parameters, feed them to the function
    for (j = 0; j < par[0]; j++) {
      X[j][0] = x[j]-par[off+1+j];
      for (k = 0; k < par[0]; k++)
          icov[j][k] = par[off+1+(int)par[0]+j*(int)par[0]+k];
    }
    Xt = X.Transpose();

    val += exp(-(Xt*icov*X)[0][0]/2)/pow(2*M_PI, par[0]/2)/sqrt(par[off])/par[1];
  }

  return val;
}

// ROOT style function necessary to produce a 3D contour, must be at global level
double RMultiGausOffset(double *x, double *par) {

  return RMultiGaus(x, &par[1])-par[0];
}

// ROOT style function necessary to numerically compute the optimal radii of each peak
double RMultiGausPContent(double *x, double *par) {

  // The input of this function is the level of the function (from 0 to max of tallest peak)
  size_t i;
  double alpha(0);
  for (i = 0; i < par[1]; i++) {

    // A peak contributes if its max if larger than the current level, radius 0 otherwise
    if ( x[0] > par[2+i] ) {
      par[2+(int)par[1]+i] = 0;
      continue;
    }

    // If the level if too close to zero, the radius is infinite and alpha is 1/k
    if ( x[0] < 1e-9 ) {
      alpha += 1./par[1];
      continue;
    }

    // The current radius can be simply found as
    par[2+(int)par[1]+i] = sqrt(-2*log(x[0]/par[2+i]));

    // Given the radius, increment the p-value
    alpha += TMath::Gamma(par[0]/2, par[2+(int)par[1]+i]*par[2+(int)par[1]+i]/2)/par[1];
  }

  // Return the p-value to optimize
  return alpha;
}

DMultiGaus::DMultiGaus(const size_t n,
		       const size_t k,
		       const double dist) :
  _k(k), _means(k), _covs(k), _invcovs(0), _Js(0), _dets(0) {

  try {
    // Assert that the dimension and number of peaks are non-zero and that they are not superimposed
    Assert::IsNonZero("DMultiGaus::DMultiGaus", "Dimension", n);
    Assert::IsNonZero("DMultiGaus::DMultiGaus", "Number of peaks", k);
    Assert::IsNonZero("DMultiGaus::DMultiGaus", "Distance between peaks", dist);

    // One can produce a maximum of n+1 mutually equidistant points in nD (standard n-simplex)
    // If k is smaller or equal to n+1, dispose them that way. For a higher dimension, it is the 
    // open problem of placing k points optimally on an n-sphere (TODO?), place periodically for now
    _dim = n;
    size_t i, j;
    if ( _k <= _dim+1 ) {
      // Get the vertices of the k points on the (k-1)-unit sphere
      _means = Math::SimplexVertices(k-1);

      // If k is smaller than n+1, fill the rest of the coordinates with 0.
      if ( k < _dim+1 )
        for (i = 0; i < k; i++)      
          for (j = k-1; j < _dim; j++)
	      _means[i].push_back(0.);

      // Scale the mean by the dist factor (distance between two points is sqrt(2.+2./n)
      for (i = 0; i < k; i++)      
        for (j = 0; j < _dim; j++)
	    _means[i][j] *= dist/sqrt(2.+2./n);

    } else {
      // Offset the means by dist in 1 dimension only
      for (i = 0; i < k; i++) {
        _means[i].resize(_dim);
        _means[i][0] = i*dist;
      }

    }

    // The covariance matrices are all identity matrices for normals
    for (i = 0; i < k; i++) {
      _covs[i].Resize(n, n);
      _covs[i].Identity();
    }

    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "DMultiGaus::DMultiGaus"));
  }
}

DMultiGaus::DMultiGaus(const std::vector<double>& mu,
	     	       const std::vector<double>& sigma) :
  _k(mu.size()), _means(mu.size()), _covs(mu.size()), _invcovs(0), _Js(0), _dets(0) {

  try {
    Assert::IsNonZero("DMultiGaus::DMultiGaus", "Number of peaks", mu.size());
    Assert::SameSize("DMultiGaus::DMultiGaus", "Mean and width vectors", mu, sigma);
    _dim = 1;
    size_t i;
    for (i = 0; i < _k; i++)
        _means[i].push_back(mu[i]);

    for (i = 0; i < _k; i++) {
      Assert::IsGreater("DMultiGaus::DMultiGaus", "sigma["+std::to_string((int)i)+"]", sigma[i], 0.);
      _covs[i].Resize(1, 1);
      _covs[i][0][0] = sigma[i]*sigma[i];
    }

    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "DMultiGaus::DMultiGaus"));
  }
}

DMultiGaus::DMultiGaus(const std::vector<std::vector<double>>& mu,
	     	       const std::vector<std::vector<double>>& sigma) :
  _k(mu.size()), _means(mu), _covs(mu.size()), _invcovs(0), _Js(0), _dets(0) {

  try {
    Assert::IsNonZero("DMultiGaus::DMultiGaus", "Number of peaks", mu.size());
    Assert::IsNonZero("DMultiGaus::DMultiGaus", "Dimension", mu[0].size());
    Assert::SameSize("DMultiGaus::DMultiGaus", "Mean and width vectors", mu, sigma);
    _dim = mu[0].size();
    size_t i, j;
    std::string id;
    for (i = 0; i < _dim; i++) {
      id = std::to_string((int)i);
      Assert::IsEqual("DMultiGaus::DMultiGaus",
		"mean["+id+"] size and dimension", mu[i].size(), _dim);
      Assert::IsEqual("DMultiGaus::DMultiGaus",
		"sigma["+id+"] size and dimension", sigma[i].size(), _dim);
      _covs[i].Resize(_dim, _dim);
      for (j = 0; j < _dim; j++) {
        Assert::IsGreater("DMultiGaus::DMultiGaus",
		"sigma["+id+"]["+std::to_string((int)j)+"]", sigma[i][j], 0.);
        _covs[i][j][j] = sigma[i][j]*sigma[i][j];
      }
    }

    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "DMultiGaus::DMultiGaus"));
  }
}

DMultiGaus::DMultiGaus(const std::vector<std::vector<double>>& mu,
	     	       const std::vector<Matrix<double>>& cov) :
  _k(mu.size()), _means(mu), _covs(cov), _invcovs(0), _Js(0), _dets(0) {

  try {
    Assert::IsNonZero("DMultiGaus::DMultiGaus", "Number of peaks", mu.size());
    Assert::IsNonZero("DMultiGaus::DMultiGaus", "Dimension", mu[0].size());
    Assert::IsEqual("DMultiGaus::DMultiGaus",
	"Mean and covariance vectors sizes", mu.size(), cov.size());

    _dim = mu[0].size();
    size_t i;
    std::string id;
    for (i = 0; i < _dim; i++) {
      id = std::to_string((int)i);
      Assert::IsEqual("DMultiGaus::DMultiGaus",
		"mean["+id+"] size and dimension", mu[i].size(), _dim);
      Assert::IsEqual("DMultiGaus::DMultiGaus",
		"cov["+id+"] size and dimension", mu[i].size(), _dim);
      Assert::IsSquare("DMultiGaus::DMultiGaus", "cov["+id+"]", cov[i].std());
    }

    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "DMultiGaus::DMultiGaus"));
  }
}

DMultiGaus::DMultiGaus(const DMultiGaus& df) {
  *this = df;
}

DMultiGaus& DMultiGaus::operator=(const DMultiGaus& df) {
  if ( this == &df )
      return *this;

  _k = df._k;
  _means = df._means;
  _covs = df._covs;
  _invcovs = df._invcovs;
  _Js = df._Js;
  _dets = df._dets;

  return *this;
}

DMultiGaus::~DMultiGaus () {}

double DMultiGaus::CDF(const double& v) const {

  // CDF only defined for 1D
  Assert::IsEqual("DMultiGaus::CDF", "Dimension", _dim, (size_t)1);

  // Add up the cumulative density functions of each of the peaks
  double cdf(0.);
  size_t i;
  for (i = 0; i < _k; i++)
      cdf += .5*(1.+erf((v-_means[i][0])/(sqrt(2*_dets[i]))))/_k;

  return cdf;
}

double DMultiGaus::CDFRadial(const double& R) const {

  throw(Exceptions::Exception(Exceptions::recoverable,
	"No radial symmetry for any definition of the p-norm of a vector",
	"DMultiGaus::CDFRadial"));
  return 0.;
}

bool DMultiGaus::AddPeak(const std::vector<double>& mu, const std::vector<double>& sigma) {

  Assert::SameSize("DMultiGaus::AddPeak", "Mean and width vectors", mu, sigma);
  Matrix<double> cov;
  size_t i;
  for (i = 0; i < _dim; i++)
      cov[i][i] = sigma[i]*sigma[i];

  return AddPeak(mu, cov);
}

bool DMultiGaus::AddPeak(const std::vector<double>& mu, const Matrix<double>& cov) {

  // Throw if incompatible with current set
  Assert::SameSize("DMultiGaus::AddPeak", "Peak dimension and function dimension", mu, _means[0]);

  // Add the mean and cov to the list, increment k
  _k++;
  _means.push_back(mu);
  _covs.push_back(cov);

  try {
    // Compute and add the deterimant to the list
    _dets.push_back(cov.Determinant());

    // Compute and add the inverse covariance matrix to the list
    _invcovs.push_back(cov.Inverse());
  
    // Compute and add the Jacobian to the list
    std::vector<double> lambda(_dim);
    Matrix<double> U, Lsqrt(_dim, _dim);
    U = cov.EigenVectors(lambda);
    size_t i;
    for (i = 0; i < _dim; i++)
       Lsqrt[i][i] = sqrt(lambda[i]);

    _Js.push_back(U*Lsqrt);
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not increment"+std::string(e.what()),
	  "DMultiGaus::AddPeak"));
  }

  // Extend the ranges if necessary
  size_t i;
  for (i = 0; i < _dim; i++) {
    if ( mu[i]-5*sqrt(cov[i][i]) < _lower[i] )
	_lower[i] = mu[i]-5*sqrt(cov[i][i]);
    if ( mu[i]+5*sqrt(cov[i][i]) > _upper[i] )
	_upper[i] = mu[i]+5*sqrt(cov[i][i]);
  }
  
  return true;
}

std::vector<TLine*> DMultiGaus::Contour1D(const double alpha) const {

  // Add a line for each of the peaks that have a non-zero radii, the level is always the same
  Assert::IsProbability("DMultiGaus::Contour1D", "alpha", alpha);
  std::vector<double> radii = Radii(alpha);
  double level = Level(alpha);
  std::vector<TLine*> lines;
  size_t i;
  for (i = 0; i < _k; i++) {
    if ( !radii[i] )
	continue;

    lines.push_back(new TLine(_means[i][0]-radii[i], level, _means[i][0]+radii[i], level));
    lines.back()->SetLineWidth(2);
  }

  return lines;
}

std::vector<TObject*> DMultiGaus::Contour2D(const double alpha) const {

  // Add an ellipse for each of the peaks that have a non-zero radii
  Assert::IsProbability("DMultiGaus::Contour2D", "alpha", alpha);
  std::vector<double> radii = Radii(alpha);
  std::vector<TObject*> ells;
  size_t i;
  for (i = 0; i < _k; i++) {
    if ( !radii[i] )
	continue;

    // Define the ellipse
    TEllipse *ell = new TEllipse();
    ell->SetFillColorAlpha(0, 0);
    ell->SetLineWidth(2);

    // First set the centre of the ellipse at the mean
    ell->SetX1(_means[i][0]);
    ell->SetY1(_means[i][1]);

    // Compute the eigenvalues of the covariance matrix as they are related to the axes dimensions
    std::vector<double> lambdas;
    Matrix<double> V = _covs[i].EigenVectors(lambdas);
    ell->SetR1(sqrt(lambdas[0]*radii[i]));
    ell->SetR2(sqrt(lambdas[1]*radii[i]));

    // Compute the eigenvectors to yield the rotation angle of the ellipse
    ell->SetTheta(atan2(V[1][0], V[0][0])*360./(2*M_PI));
    ell->SetPhimin(0);
    ell->SetPhimax(360);
    ells.push_back(ell);
  }

  return ells;
}

TF3* DMultiGaus::Contour3D(const double alpha) const {

  // Find the level of the isosurface for the provided alpha
  Assert::IsProbability("DMultiGaus::Contour3D", "alpha", alpha);
  double C = Level(alpha);

  // Define and return the offset function
  std::vector<double> pars = OffsetParameters(C);
  TF3* multiell = new TF3("multiell", RMultiGausOffset, _lower[0], _upper[0], _lower[1],
			  _upper[1], _lower[2], _upper[2], pars.size());
  multiell->SetTitle(_title.c_str());
  multiell->SetParameters(&(pars[0]));
  multiell->SetFillColorAlpha(1, .25);
  multiell->SetLineColor(1);
  return multiell;
}

std::vector<TLine*> DMultiGaus::ContourRadial(const double alpha) const {

  throw(Exceptions::Exception(Exceptions::recoverable,
	"No radial symmetry for any definition of the p-norm of a vector",
	"DMultiGaus::ContourRadial"));
  return {new TLine()};
}

double DMultiGaus::ContourVolume(const double alpha) const {

  // Get the radii and increment the volume accordingly for each peak
  Assert::IsProbability("DMultiGaus::ContourVolume", "alpha", alpha);
  std::vector<double> radii = Radii(alpha);
  double vol(0.);
  size_t i;
  for (i = 0; i < _k; i++) {
    // In lower dimensions, use simple formulas
    if ( _dim == 1 ) {
      vol += 2*radii[i]*sqrt(_dets[i]);
      continue;
    } else if ( _dim == 2 ) {
      vol += M_PI*radii[i]*radii[i]*sqrt(_dets[i]);
      continue;
    } else if ( _dim == 3 ) {
      vol += (4./3)*M_PI*radii[i]*radii[i]*radii[i]*sqrt(_dets[i]);
      continue;
    }

    // In higher dimension, simply use the general definition for an n-ellipse
    vol += pow(M_PI, (double)_dim/2)*pow(radii[i], _dim)
	   * sqrt(_dets[i])/tgamma((double)_dim/2+1);
  }

  return vol;
}

double DMultiGaus::Evaluate(const std::vector<double>& v) const {

  Assert::SameSize("DMultiGaus::Evaluate", "Function space and argument", _lower, v);

  // Each of the peak a fraction 1/k of the total probability, loop and increment
  double val(0.);
  size_t i;
  for (i = 0; i < _k; i++) {

    // Clear the simple case of n = 1
    if ( _dim == 1 ) {
      val += exp(-pow(v[0]-_means[i][0], 2)/_dets[i]/2)/sqrt(2*M_PI*_dets[i])/_k;
      continue;
    }

    // If higher dimension, need to define vectors
    Matrix<double> X(v), M(_means[i]), Xt;
    X = X-M;
    Xt = X.Transpose();
    val += exp(-(Xt*_invcovs[i]*X)[0][0]/2)/pow(2*M_PI, (double)_dim/2)/sqrt(_dets[i])/_k;
  }

  return val;
}

double DMultiGaus::Level(const double alpha) const {

  Assert::IsProbability("DMultiGaus::Level", "alpha", alpha);
  std::vector<double> radii = Radii(alpha);
  return Radial(0, radii[0]);
}

double DMultiGaus::Norm() const {

  // Clear the simple case of n = 1
  if ( _dim == 1 )
      return Norm1D(0);

  // For higher dimensions, it is a tricky analytical challenge to integrate a correlated
  // gaussian over a orthotopical interval. One needs to use a numerical method,
  // compute the integral over the interval, correct to have I = 1.

  // If n = 2 or 3, feed it to ROOT
  std::vector<double> pars = Parameters();
  if ( _dim == 2 ) {
    TF2 f("", RMultiGaus, _lower[0], _upper[0], _lower[1], _upper[1], pars.size());
    f.SetParameters(&(pars[0]));
    return f.Integral(_lower[0], _upper[0], _lower[1], _upper[1]);
  } else if ( _dim == 3 ) {
    TF3 f("", RMultiGaus, _lower[0], _upper[0], _lower[1], _upper[1], _lower[2], _upper[2], pars.size());
    f.SetParameters(&(pars[0]));
    return f.Integral(_lower[0], _upper[0], _lower[1], _upper[1], _lower[2], _upper[2]);
  }
  
  // For higher dimensions, no function available in ROOT, only support uncorrelated
  // Check that the matrices are all diagonal, if not, throw recoverable
  size_t i, j, k;
  for (k = 0; k < _k; k++)
    for (i = 0; i < _dim; i++)
      for (j = 0; j < _dim; j++)
        if ( i != j && _covs[k][i][j] )
	    throw(Exceptions::Exception(Exceptions::recoverable,
	          "Matrix "+std::to_string((int)k)+" is correlated, norm not supported",
		  "DMultiGaus::Norm"));

  // If uncorrelated, the norm is simply the product of the norms in each dimension
  double norm(1.);
  for (i = 0; i < _dim; i++)
      norm *= Norm1D(i);

  return norm;
}

double DMultiGaus::Radial(const double R) const {

  throw(Exceptions::Exception(Exceptions::recoverable,
	"No radial symmetry for any definition of the p-norm of a vector",
	"DMultiGaus::Radial"));
  return 0.;
}

double DMultiGaus::Radial(const size_t i,
			  const double R) const {

  Assert::IsGreater("DMultiGaus::Radial", "Radius", R, 0.);
  return exp(-R*R/2)/pow(2*M_PI, (double)_dim/2.)/sqrt(_dets[i])/_k;
}

std::vector<double> DMultiGaus::Radii(const double alpha) const {

  Assert::IsProbability("DMultiGaus::Radii", "alpha", alpha);

  // Start an array of parameters to feed to the minimizer
  std::vector<double> pars(2+2*_k);
  pars[0] = _dim;
  pars[1] = _k;

  // Fetch the maximums of each of the peaks, set them as parameters
  size_t i;
  double max(0);
  for (i = 0; i < _k; i++) {
    pars[2+i] = 1./pow(2*M_PI, (double)_dim/2)/sqrt(_dets[i])/_k;
    if ( pars[2+i] > max )
	max = pars[2+i];
  }

  // Use numerical method to find the optimal level so that the intersection of a plane at
  // that level with the MultiGaus function exactly contains a fraction alpha of the probability
  TF1 f("", RMultiGausPContent, 0, max, pars.size());
  f.SetParameters(&(pars[0]));
  f.GetX(alpha);

  // Fetch the values of the radii from the minimizer
  std::vector<double> radii(_k);
  for (i = 0; i < _k; i++)
      radii[i] = f.GetParameter(2+_k+i);

  return radii;
}

double DMultiGaus::Random() {

  // As each peak carries the same probability ,choose one randomly and then sample from it
  size_t rk = rand()%_k;
  return _rdmzer.Gaus(_means[rk][0], sqrt(_covs[rk][0][0]));
}

std::vector<double> DMultiGaus::RandomVector() {

  // First generate a vector of uncorrelated normal variables
  Vector<double> x(_dim);
  size_t i;
  for (i = 0; i < _dim; i++)
      x[i] = _rdmzer.Gaus(0., 1);

  // As each peak carries the same probability, choose one randomly and then sample from it
  // Convert the coordinates in the reference of the n-sphere to coordinates in the
  // metric of the covariance matrix using the precomputed Jacobian
  size_t rk = rand()%_k;
  x = _Js[rk]*x;

  // Return the answer in the form of a vector, offset by the mean
  std::vector<double> v(_dim);
  for (i = 0; i < _dim; i++)
      v[i] = _means[rk][i]+x[i];
  return v;
}

bool DMultiGaus::Initialize() {

  // Initialize the default name of the function and its title
  _name = std::to_string(_k)+"gaus"+std::to_string(_dim);
  _title = std::to_string(_dim)+"D "+std::to_string(_k)+"-modal Gaussian";

  // Compute the determinants of the covariance matrices
  size_t i, j;
  for (i = 0; i < _k; i++)
      _dets.push_back(_covs[i].Determinant());

  // Compute the inverse covariance matrices
  for (i = 0; i < _k; i++)
      _invcovs.push_back(_covs[i].Inverse());

  // Compute the Jacobians
  std::vector<double> lambda(_dim);
  Matrix<double> U, Lsqrt(_dim, _dim);
  for (i = 0; i < _k; i++) {
    U = _covs[i].EigenVectors(lambda);
    for (j = 0; j < _dim; j++)
       Lsqrt[j][j] = sqrt(lambda[j]);

    _Js.push_back(U*Lsqrt);
  }

  // Set up a default range. Compute mean +/- 5*sigma for each peak, find the lowest and
  // highest values, set them as the range limits
  std::vector<double> mins(_k), maxs(_k);
  for (j = 0; j < _dim; j++) {
    for (i = 0; i < _k; i++) {
      mins[i] = _means[i][j]-5*sqrt(_covs[i][j][j]);
      maxs[i] = _means[i][j]+5*sqrt(_covs[i][j][j]);
    }

    _lower.push_back(Math::Min(mins));
    _upper.push_back(Math::Max(maxs));
  }

  // Initialize the pseudorandom number generator
  time_t tt;
  _rdmzer = TRandom3(time(&tt));

  return true;
}

double DMultiGaus::Norm1D(const size_t i) const {

  // Each peak carries the same amount 1/k of total probability
  double norm(1.);
  double Rl, Ru, tempnorm;
  size_t j;
  for (j = 0; j < _k; j++) {
    // Compute the central normal variables
    Rl = (_lower[i]-_means[j][i])/sqrt(_covs[j][i][i]);
    Ru = (_upper[i]-_means[j][i])/sqrt(_covs[j][i][i]);

    // Get the norm of this particular peak, decrement
    tempnorm = .5*(erf(Ru/sqrt(2)) - erf(Rl/sqrt(2)));
    norm -= (1.-tempnorm)/_k;
  }

  return norm;
}

std::vector<double> DMultiGaus::Parameters() const {

  std::vector<double> pars;
  pars.push_back(_dim);
  pars.push_back(_k);

  size_t i, j, k;
  for (i = 0; i < _k; i++) {
    pars.push_back(_dets[i]);

    for (j = 0; j < _dim; j++)
        pars.push_back(_means[i][j]);

    for (j = 0; j < _dim; j++)
      for (k = 0; k < _dim; k++)
	  pars.push_back(_invcovs[i][j][k]);
  }

  return pars;
}

std::vector<double> DMultiGaus::OffsetParameters(const double C) const {

  std::vector<double> pars = Parameters();
  pars.insert(pars.begin(), C);
  return pars;
}
