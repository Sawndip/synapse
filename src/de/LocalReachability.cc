#include "LocalReachability.hh"

LocalReachability::LocalReachability() :
  _dim(0), _points(0), _k(0), _kdtree(NULL) {
}

LocalReachability::LocalReachability(const std::vector<std::vector<double>>& points, const size_t k) :
  _dim(0), _points(points), _k(k), _kdtree(NULL) {

  try {
    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "LocalReachability::LocalReachability"));
  }
}

LocalReachability::LocalReachability(const LocalReachability& lr) {
  *this = lr;
}

LocalReachability& LocalReachability::operator=(const LocalReachability& lr) {
  if ( this == &lr )
      return *this;

  _dim = lr._dim;
  _points = lr._points;
  _k = lr._k;

  // Special care has to be taken here as NanoFLANN does not have a copy constructor
  if ( _kdtree )
      delete _kdtree;
  if ( _dim && _points.size() ) {
    _kdtree = new KDTree(_dim, _points, 10);
    _kdtree->index->buildIndex();
  }

  return *this;
}

LocalReachability::~LocalReachability () {

  delete _kdtree;
}

void LocalReachability::Initialize() {

  // Check that k is of a sensible value
  Assert::IsNonZero("LocalReachability::Initialize", "Number of neighbours", _k);
  if ( _k > _points.size() )
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	    "Number of neighbours requested is greater than the number of input points",
	    "LocalReachability::Initialize"));

  _dim = _points[0].size();

  // Initialize the k-d tree
  _kdtree = new KDTree(_dim, _points, 10);
  _kdtree->index->buildIndex();
}

double LocalReachability::Evaluate(const double& v) const {

  return Evaluate(std::vector<double>({v}));
}

double LocalReachability::Evaluate(const double* v) const {

  std::vector<double> vv(v, v+_dim);
  return Evaluate(vv);
}

double LocalReachability::Distance(const std::vector<double>& A, const std::vector<double>& B) const {

  double sdist(0.);
  size_t i;
  for (i = 0; i < _dim; i++)
      sdist += pow(B[i]-A[i], 2);

  return sqrt(sdist);
}

double LocalReachability::KDistance(const std::vector<double>& A) const {

  // Look into the indexed points of the k-d tree of the k closest points
  std::vector<size_t> indices(_k);
  std::vector<double> dists(_k);
  nanoflann::KNNResultSet<double> results(_k);
  results.init(&indices[0], &dists[0]);
  _kdtree->index->findNeighbors(results, &A[0], nanoflann::SearchParams(10));

  // Return the distance to the k^th point (NanoFLANN returns squared distance)
  return sqrt(dists[_k-1]);
}

double LocalReachability::RDistance(const std::vector<double>& A, const std::vector<double>& B) const {

  return std::max(Distance(A, B), KDistance(B));
}

double LocalReachability::RVolume(const double& radius) const {

  // Compute the distance, return the volume of the hypersphere
  size_t n = _points[0].size();
  return pow(sqrt(M_PI)*radius, n)/tgamma((double)n/2+1);
}

double LocalReachability::RDensity(const std::vector<double>& A) const {

  // Look into the indexed points of the k-d tree of the k closest points
  std::vector<size_t> indices(_k);
  std::vector<double> dists(_k);
  nanoflann::KNNResultSet<double> results(_k);
  results.init(&indices[0], &dists[0]);
  _kdtree->index->findNeighbors(results, &A[0], nanoflann::SearchParams(10));

  // The LRD is the inverse of the average reachability of A
  double rdist(0.);
  size_t i;
  for (i = 0; i < _k; i++)
      rdist += RDistance(A, _points[indices[i]])/_k;
  double vol = RVolume(rdist);
  return _k/vol/_points.size();
}
