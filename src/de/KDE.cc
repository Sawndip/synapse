#include "KDE.hh"

KDE::KDE() :
  _btype(""), _H(0, 0), _samples(0.), _kernel() {
}

KDE::KDE(std::vector<std::vector<double>> mat) :
  _btype("silverman"), _H(0, 0), _samples(mat), _kernel() {

  try {
    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "KDE::KDE"));
  }
}

KDE::KDE(const KDE& kde) {
  *this = kde;
}

KDE& KDE::operator=(const KDE& kde) {
  if ( this == &kde )
      return *this;

  _btype = kde._btype;
  _H = kde._H;
  _samples = kde._samples;
  _kernel = kde._kernel;

  return *this;
}

KDE::~KDE () {}

double KDE::Evaluate(const std::vector<double>& v) const {

  // Loop over the whole sample and increment the KDE
  size_t d = _samples.size();
  size_t n = _samples[0].size();

  std::vector<double> vdiff;
  double kde(0.);
  size_t i, j;
  for (i = 0; i < n; i++) {

    // Get the difference between the argument and the sample point
    vdiff.resize(d);
    for (j = 0; j < d; j++)
	vdiff[j] = v[j] - _samples[j][i];

    // Increment
    kde += _kernel(vdiff);
  }

  // Return the KDE with its factor
  return kde/n;
}

double KDE::Evaluate(const double& v) const {

  return Evaluate(std::vector<double>({v}));
}

double KDE::Evaluate(const double* v) const {

  std::vector<double> vv(v, v+_samples.size());
  return Evaluate(vv);
}

std::vector<double> KDE::FastFunction(std::vector<std::vector<double>> points,
				      double eps) {

  // Check the input array of points
  size_t i, j, k, l;
  Assert::IsEqual("KDE::FastFunction",
	"Test dimension and training dimension", points.size(), _samples.size());
  Assert::IsNonZero("KDE::FastFunction", "Number of test points", points[0].size());

  // First step is to identify the range in each axis, data and evaluation points
  size_t d = points.size();
  size_t n = _samples[0].size();
  size_t m = points[0].size();
  std::pair<double, double> mm_points, mm_samples;
  std::vector<std::pair<double, double>> minmax(d);
  for (i = 0; i < d; i++) {
    auto it_points = std::minmax_element(points[i].begin(), points[i].end());
    auto it_samples = std::minmax_element(_samples[i].begin(), _samples[i].end());
    mm_points = std::pair<double, double>(points[i][it_points.first-points[i].begin()], 
				          points[i][it_points.second-points[i].begin()]);
    mm_samples = std::pair<double, double>(_samples[i][it_samples.first-_samples[i].begin()], 
				           _samples[i][it_samples.second-_samples[i].begin()]);

    minmax[i] = std::pair<double, double>(std::min(mm_points.first, mm_samples.first),
					  std::max(mm_points.second, mm_samples.second));
  }

  // Rescale the data and the points to lie between 0 and 1
  std::vector<std::vector<double>> samples_scaled(d);
  for (i = 0; i < d; i++) {
    samples_scaled[i].resize(n);
    for (j = 0; j < n; j++)
        samples_scaled[i][j] = (_samples[i][j]-minmax[i].first)/(minmax[i].second-minmax[i].first);
  }

  std::vector<std::vector<double>> points_scaled(d);
  for (i = 0; i < d; i++) {
    points_scaled[i].resize(m);
    for (j = 0; j < m; j++)
        points_scaled[i][j] = (points[i][j]-minmax[i].first)/(minmax[i].second-minmax[i].first);
  }

  // Scale the bandwidths accordingly
  Matrix<double> Hscaled = _H;
  for (i = 0; i < d; i++)
	Hscaled[i][i] /= pow(minmax[i].second-minmax[i].first, 2);

  // Reinitialize the Kernel with the scaled bandwidth
  _kernel = DGaus(std::vector<double>(d, 0), Hscaled);

  // Provided the requested error, define eps'
  double _factor = 1./n/sqrt(Hscaled.Determinant());
  double q = _factor/pow(2*M_PI, (double)d/2);
  double epsprime = eps/(n*q);

  // Choose each interval "radius" to be half the bandwidth in each dimension.
  std::vector<double> rx(d);
  for (i = 0; i < d; i++)
      rx[i] = sqrt(Hscaled[i][i])/2;

  // Subdivide each axis of the space in L/h equal parts.
  std::vector<size_t> L(d), step(d); 	// Amount of points on each axis and stepping
  size_t N(1); 				// Total amount of points on the grid L_0x...xL_d
  int a;
  for (a = (int)d-1; a > -1; a--) {
    // Total number of subdivisions in this projection
    // Each of the interval centre is at (l+0.5)/L, l = 0, ..., L-1
    L[a] = 1./(2*rx[a]);
    N *= L[a];

    // The first axis is the MSB, last in the LSB
    if ( a < (int)d-1 ) {
      step[a] = L[a+1] * step[a+1];
    } else {
      step[a] = 1;
    }
  }

  // Create a grid, each point of which can host an array of indices (to refer to _samples)
  // Associate every point with the closest interval centre, store in arrays for each
  // possible inteval index combination, i.e. in a dxL_0x...xL_d matrix of vectors
  std::vector<std::vector<size_t>> c_array(N);
  size_t id(0);
  for (j = 0; j < n; j++) {

    // For each projection, find the order of the interval, increment id
    id = 0;
    for (i = 0; i < d; i++)
        id += (int)(samples_scaled[i][j]/(2*rx[i]))*step[i];

    // Add the point to the corresponding cell
    c_array[id].push_back(j);
  }

  // Define the cutoff distance
  std::vector<double> rex(d);
  for (i = 0; i < d; i++)
      rex[i] = rx[i]+4*rx[i]*sqrt(log(1./epsprime));

  // Need to chose the taylor order numerically, seems complicated (TODO)
  size_t p = 1;	// Arbitrary order of the taylor expension

  // Compute the prefactor B for each one of the cells
  std::vector<std::vector<double>> B(p);
  Hscaled.Print();
  std::vector<double> cl(d);
  std::vector<double> vdiff(d);
  size_t index, temp;
  double factor;
  for (k = 0; k < p; k++) {
    B[k].resize(N);
    for (l = 0; l < N; l++) {
      B[k][l] = 0.;
      // The vector position of this centre
      temp = l;
      for (i = 0; i < d; i++) {
	  index = temp/step[i];
	  temp -= index*step[i];
          cl[i] = (index+.5)*2*rx[i];
      }

      for (j = 0; j < c_array[l].size(); j++) {
	// Produce the difference vector
        factor = 1.;
	for (i = 0; i < d; i++) {
	  vdiff[i] = samples_scaled[i][c_array[l][j]] - cl[i];
          factor *= (samples_scaled[i][c_array[l][j]] - cl[i])/sqrt(Hscaled[i][i]);
        }
        
        B[k][l] += _factor*_kernel(vdiff)*pow(factor, k);
      }
    }
  }

  // Loop over all the measurement points, if they are close enough to a c_l, increment
  std::vector<double> results(m);  
  for (j = 0; j < m; j++)
      results[j] = 0;
  double radius;
  for (l = 0; l < N; l++) {

    // The vector position of this centre
    temp = l;
    for (i = 0; i < d; i++) {
	index = temp/step[i];
	temp -= index*step[i];
        cl[i] = (index+.5)*2*rx[i];
    }      

    // Loop over the measurement points
    for (j = 0; j < m; j++) {
      // Get the vector difference between the two
      factor = 0;
      for (i = 0; i < d; i++) {
	vdiff[i] = points_scaled[i][j] - cl[i];
        factor *= (points_scaled[i][j] - cl[i])/sqrt(Hscaled[i][i]);
      }

      // If the points are apart beyond the cutoff, proceed (ellipsoid cut)
      radius = 0;
      for (i = 0; i < d; i++)
	  radius += pow(vdiff[i]/rex[i], 2);
      if ( radius > 1)
	  continue;

      // For each taylor order and dimension, fill
      for (k = 0; k < p; k++) {
	results[j] += B[k][l]*pow(2*M_PI, (double)d/2)*_kernel(vdiff)*pow(factor, k);
      }
    }
  }

  return results;
}

double KDE::FunctionROOT(double *x, double *par) {

  // The parameter specifies the number of variables
  size_t d = par[0];
  std::vector<double> v(d);
  size_t i;
  for (i = 0; i < d; i++)
      v[i] = x[i];

  // Return the KDE
  return Evaluate(v);
}

void KDE::Initialize() {

  // Check that all the samples are the same size
  Assert::IsNonZero("KDE::Initialize", "Dimension", _samples.size());
  Assert::IsNotEmpty("KDE::Initialize", "Vector of input points", _samples[0]);

  // Get the standard deviations of each samples
  size_t d = _samples.size();
  std::vector<double> vsig(d);
  for (size_t i = 0; i < d; i++)
      vsig[i] = Math::RMS(_samples[i]);

  // Set the bandwidth matrix
  double n = _samples[0].size();
  if ( _btype == "silverman" ) {
    _H = SilvermanBandwidth(vsig, n);
  } else if ( _btype == "scott" ) {
    _H = ScottBandwidth(vsig, n);
  } else {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Unknown bandwidth type: "+_btype,
	  "KDE::Initialize"));
  }

  // Set the gaussian Kernel with the produced bandwidth
  _kernel = DGaus(std::vector<double>(d, 0), _H);
}

Matrix<double> KDE::SilvermanBandwidth(std::vector<double> vsig, size_t n) {

  // Get the dimensionality of the matrix
  size_t d = vsig.size();

  // Initialize and fill the matrix
  Matrix<double> H(d, d);
  size_t i, j;
  for (i = 0; i < d; i++) {
    for (j = 0; j < d; j++) {
      if ( i == j ) {
	H[i][j] = pow(pow(4./(d+2.), 1./(d+4))*pow((double)n, -1./(d+4))*vsig[i], 2);
      } else {
	H[i][j] = 0.;
      }
    }
  }

  return H;
}

Matrix<double> KDE::ScottBandwidth(std::vector<double> vsig, size_t n) {

  // Get the dimensionality of the matrix
  size_t d = vsig.size();

  // Initialize and fill the matrix
  Matrix<double> H(d, d);
  size_t i, j;
  for (i = 0; i < d; i++) {
    for (j = 0; j < d; j++) {
      if ( i == j ) {
	H[i][j] = pow(pow((double)n, -1./(d+4))*vsig[i], 2);
      } else {
	H[i][j] = 0.;
      }
    }
  }

  return H;
}
