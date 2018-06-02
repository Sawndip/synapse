#include "Hull.hh"

Hull::Hull() :
  _dim(0), _points(), _vol(0) {
}

Hull::Hull(const std::vector<std::vector<double>>& points) :
  _dim(0), _points(), _vol(0) {

  try {
    Initialize(points);
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "Hull::Hull"));
  }
}

Hull::Hull(const std::vector<double>& points, const size_t& dim) :
  _dim(0), _points(), _vol(0) {

  this->SetPointArray(points, dim);
}

Hull::Hull(const Hull& hull) {
  *this = hull;
}

Hull& Hull::operator=(const Hull& hull) {
  if ( this == &hull )
      return *this;

  _dim = hull._dim;
  _points = hull._points;
  _vol = hull._vol;

  return *this;
}

Hull::~Hull () {}

void Hull::Initialize(const std::vector<std::vector<double>>& points) {

  // Check that there are points to fit and there are at least n+1 points for dimension n
  Assert::IsNotEmpty("Hull::Initialize", "Vector of input points", points);
  Assert::IsNonZero("Hull::Initialize", "Dimension", points[0].size());

  // Record the dimension and store the points
  _dim = points[0].size();	// Dimension of the space
  _points = points;
}

void Hull::SetPointArray(const std::vector<double>& points, const size_t& dim) {

  _dim = dim;
  size_t N = points.size()/dim;
  _points.resize(N);
  size_t i, j;
  for (i = 0; i < N; i++) {
    _points[i].resize(0);
    for (j = 0; j < dim; j++)
        _points[i].push_back(points[i*dim+j]);
  }
}

Matrix<double> Hull::GetCWPoints() const {

  // Set a new set of vectors and sort them clockwise
  std::vector<Vector<double>> vectors(_points.size());
  size_t i, j;
  for (i = 0; i < _points.size(); i++)
      vectors[i] = Vector<double>(_points[i]);
  Math::SortCW(vectors);

  // Close the loop
  vectors.push_back(vectors.front());

  // Set a matrix with them and return
  Matrix<double> cwpoints(2, vectors.size());
  for (i = 0; i < 2; i++)
    for (j = 0; j < vectors.size(); j++)
	cwpoints[i][j] = vectors[j][i];

  return cwpoints;
}

TPolyLine* Hull::Polygon() const {

  // Throw if the wrong dimension is requested
  TPolyLine* polyline = NULL;
  Assert::IsEqual("Hull::Polygons", "Dimension", _dim, (size_t)2);

  Matrix<double> cmfv = GetCWPoints();
  polyline = new TPolyLine(cmfv.Ncols(), &(cmfv[0][0]), &(cmfv[1][0]));

  return polyline;
}
