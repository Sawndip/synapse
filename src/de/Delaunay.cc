#include "Delaunay.hh"

Delaunay::Delaunay() :
  _dim(), _vertices(), _facets(), _qhull(NULL) {
}

Delaunay::Delaunay(const std::vector<Vertex>& vertices) :
  _dim(), _vertices(vertices), _facets(), _qhull(NULL) {

  try {
    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "Delaunay::Delaunay"));
  }
}

Delaunay::Delaunay(const std::vector<std::vector<double>>& vertices) :
  _dim(), _vertices(0), _facets(), _qhull(NULL) {

  size_t i;
  for (i = 0; i < vertices.size(); i++)
      _vertices.push_back(Vertex(vertices[i]));

  try {
    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "Delaunay::Delaunay"));
  }
}

Delaunay::Delaunay(const Delaunay& del) {
  *this = del;
}

Delaunay& Delaunay::operator=(const Delaunay& del) {
  if ( this == &del )
      return *this;

  _dim = del._dim;
  _vertices = del._vertices;
  _facets = del._facets;

  if ( _qhull )
      delete _qhull;
  _qhull = del._qhull;

  return *this;
}

Delaunay::~Delaunay () {

//  delete _qhull;
}

bool Delaunay::Initialize() {

  // Check that there are vertices to fit and there are at least n+1 vertices for dimension n
  Assert::IsNotEmpty("Delaunay::Initialize", "Vector of input points", _vertices);
  Assert::IsNonZero("Delaunay::Initialize", "Dimension", _vertices[0].size());
  if ( _vertices[0].size() > _vertices.size()-1 )
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	    "Not enough points to triangulate, need at least dim+1",
	    "Delaunay::Initialize"));

  // Fill a vector with the right structure for the Qhull code
  size_t N = _vertices.size();  // Number of vertices
  _dim = _vertices[0].size();	// Dimension of the space
  std::vector<double> points;
  size_t i, j;
  for (i = 0; i < N; i++)
      for (j = 0; j < _dim; j++)
          points.push_back(_vertices[i][j]);

  // Feed the vector to Qhull and compute the Delaunay triangulation ("d")
  // Must use a pointer, the Qhull object owns the memory of its subcomponents, cannot be copied
  // TODO, the case of N=n+1 must be handled differently (not enough points for a Qhull!)
  _qhull = new orgQhull::Qhull("", _dim, N, &(points[0]), "d");

  // Store the facets for fast access. Each facet is an array of n+1 n-points (simplex)
  std::vector<Vertex> tempfacet(_dim+1);
  for (const orgQhull::QhullFacet& facet : _qhull->facetList().toStdVector()) {
    for (i = 0; i < _dim+1; i++)
	tempfacet[i] = _vertices[facet.vertices()[i].point().id()];

    _facets[facet.id()] = tempfacet;
  }

  return true;
}

const std::vector<Vertex>& Delaunay::GetBestFacet(const std::vector<double>& v, bool* isin) const {

  // First raise the point on the parabole as requested by Qhull
  std::vector<double> vr = v;
  qh_setdelaunay(_qhull->qh(), _dim+1, 1, &(vr[0]));

  // Feed the extended vector to the findbestfacet method, extract the set of vertices
  double bestdist;
  unsigned int in;
  facetT *bestf = qh_findbestfacet(_qhull->qh(), &vr[0], qh_NOupper, &bestdist, &in);

  // QHull only knows if the point is in the circumcircle of the Delaunay facet, need to check (TODO)
  *isin = in;
  /* *isin = false;
  std::vector<std::vector<double>> simplex(_dim+1);
  for (size_t i = 0; i < _dim+1; i++)
      simplex[i] = _facets.at(bestf->id)[i].GetCoordinates();
  if ( Math::IsInsideSimplex(v, simplex) )
      *isin = true;*/

  return _facets.at(bestf->id);
}

std::vector<TPolyLine*> Delaunay::Polygons() const {

  // Throw if the wrong dimension is requested
  std::vector<TPolyLine*> polylines;
  Assert::IsEqual("Delaunay::Polygons", "Dimension", _dim, (size_t)2);

  // Loop over the facets, define a polygon per facet
  std::vector<double> x, y;
  for (const std::pair<size_t, std::vector<Vertex>>& facet : _facets) {

    x.resize(0);
    y.resize(0);
    for (const Vertex& vertex : facet.second) {
      x.push_back(vertex[0]);
      y.push_back(vertex[1]);
    }
    x.push_back(x.front()); // Close the shape
    y.push_back(y.front()); // Close the shape

    polylines.push_back(new TPolyLine(4, &(x[0]), &(y[0])));
  }

  return polylines;
}
