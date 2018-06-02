#include "KNearestNeighbours.hh"

KNearestNeighbours::KNearestNeighbours() :
  _dim(0), _points(0), _k(0), _kdtree(NULL), _rotate(false), _metric(), _scale(1.) {
}

KNearestNeighbours::KNearestNeighbours(const std::vector<std::vector<double>>& points,
				       const size_t k,
				       const bool rotate) :
  _dim(0), _points(points), _k(k), _kdtree(NULL), _rotate(rotate), _metric(), _scale(1.) {

  try {
    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "KNearestNeighbours::KNearestNeighbours"));
  }
}

KNearestNeighbours::KNearestNeighbours(const KNearestNeighbours& knn) {
  *this = knn;
}

KNearestNeighbours& KNearestNeighbours::operator=(const KNearestNeighbours& knn) {
  if ( this == &knn )
      return *this;

  _dim = knn._dim;
  _points = knn._points;
  _k = knn._k;
  _rotate = knn._rotate;
  _metric.ResizeTo(knn._metric);
  _metric = knn._metric;
  _scale = knn._scale;

  // Special care has to be taken here as NanoFLANN does not have a copy constructor
  // Simply reinitialize the KD tree when a copy is made
  if ( _kdtree )
      delete _kdtree;
  if ( _dim && _points.size() ) {
    _kdtree = new KDTree(_dim, _points, 10);
    _kdtree->index->buildIndex();
  }

  return *this;
}

KNearestNeighbours::~KNearestNeighbours () {

  delete _kdtree;
}

void KNearestNeighbours::Initialize() {

  // Check that k is of a sensible value
  Assert::IsNonZero("KNearestNeighbours::Initialize", "Number of neighbours", _k);
  if ( _k > _points.size() )
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	    "Number of neighbours requested is greater than the number of input points",
	    "KNearestNeighbours::Initialize"));

  // Get the dimension
  _dim = _points[0].size();

  // If the cov flag is true, find the metric of the points and rotate them
  if ( _rotate ) {
    // Get the covariance matrix
    Matrix<double> covmat = Math::CovarianceMatrix(Matrix<double>(_points).Transpose());
    TMatrixD S(_dim, _dim);
    size_t i, j;
    for (i = 0; i < _dim; i++)
      for (j = 0; j < _dim; j++)
	  S[i][j] = covmat[i][j];

    // Find the scaling factor
    _scale = sqrt(S.Determinant());

    // Find the eigenvalues and eigenvectors of the covariance matrix
    TVectorD lambda;
    TMatrixD U = S.EigenVectors(lambda);
    U.T();

    // Find the metric to rotate the points back to the normal distribution
    TMatrixD Lsqrt(_dim, _dim);
    for (i = 0; i < _dim; i++)
       Lsqrt[i][i] = 1./sqrt(lambda[i]);

    _metric.ResizeTo(_dim, _dim);
    _metric = Lsqrt*U;

    // Rotate the points
    TMatrixD X(_dim, 1);
    for (i = 0; i < _points.size(); i++) {
      for (j = 0; j < _dim; j++)
          X[j][0] = _points[i][j];
      X = _metric*X;
      for (j = 0; j < _dim; j++)
	  _points[i][j] = X[j][0];
    }
  }

  // Initialize the k-d tree
  _kdtree = new KDTree(_dim, _points, 10);
  _kdtree->index->buildIndex();
}

double KNearestNeighbours::Evaluate(const std::vector<double>& v) const	{

  // Rotate the vector in the same metric of points if there is one
  std::vector<double> u = v;
  if ( _rotate ) {
    TMatrixD X(_dim, 1);
    for (size_t i = 0; i < _dim; i++)
	X[i][0] = u[i];
    X = _metric*X;
    for (size_t i = 0; i < _dim; i++)
	u[i] = X[i][0];
  }
  return Density(u)/_scale;
}

double KNearestNeighbours::Evaluate(const double& v) const {

  return Evaluate(std::vector<double>({v}));
}

double KNearestNeighbours::Evaluate(const double* v) const {

  std::vector<double> vv(v, v+_dim);
  return Evaluate(vv);
}

double KNearestNeighbours::Distance(const std::vector<double>& point) const {

  // Look into the indexed points of the k-d tree of the k closest points
  std::vector<size_t> indices(_k);
  std::vector<double> dists(_k);
  nanoflann::KNNResultSet<double> results(_k);
  results.init(&indices[0], &dists[0]);
  _kdtree->index->findNeighbors(results, &point[0], nanoflann::SearchParams(10));

  // Return the distance to the k^th point (NanoFLANN returns squared distance)
  return sqrt(dists[_k-1]);
}

double KNearestNeighbours::Volume(const std::vector<double>& point) const {

  // Compute the distance, return the volume of the hypersphere
  double dist = Distance(point);
  size_t n = point.size();
  return pow(sqrt(M_PI)*dist, n)/tgamma((double)n/2+1);
}

double KNearestNeighbours::Density(const std::vector<double>& point) const {

  // Compute the volume, return the density in the hypersphere
  double vol = Volume(point);
  return (double)_k/vol/_points.size();
}
