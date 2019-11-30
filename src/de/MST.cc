#include "MST.hh"

MST::MST() :
  _dim(0), _points(), _edges(), _lengths(), _length(0.) {
}

MST::MST(const std::vector<std::vector<double>>& points) :
  _dim(0), _points(points), _edges(), _lengths(), _length(0.) {

  try {
    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "MST::MST"));
  }
}

MST::MST(const MST& mst) {
  *this = mst;
}

MST& MST::operator=(const MST& mst) {
  if ( this == &mst )
      return *this;

  _dim = mst._dim;
  _points = mst._points;
  _edges = mst._edges;
  _lengths = mst._lengths;
  _length = mst._length;

  return *this;
}

MST::~MST () {}

void MST::Initialize() {

  // Check that there are vertices to fit
  Assert::IsNotEmpty("MST::Initialize", "Vector of input points", _points);
  Assert::IsNonZero("MST::Initialize", "Dimension", _points[0].size());
  _dim = _points[0].size();
  size_t n = _points.size();

  // Build a dictionary of all the n*(n-1)/2 possible edges and their lengths
  std::vector<std::pair<EdgeIndex, double>> edges;
  size_t i, j;
  for (i = 0; i < n; i++)
    for (j = i+1; j < n; j++)
        edges.push_back({EdgeIndex({i, j}), Distance(i, j)});

  // Sort the edges by increasing order of size
  std::sort(edges.begin(), edges.end(),
	    [] (const std::pair<EdgeIndex, double>& a, const std::pair<EdgeIndex, double>& b) {
		      return a.second < b.second;
	    });

  // Here we use Kruskal's algorithm with the disjoint set data structure.
  // Start assembling a path using the smallest available edges.
  // If a new edge links two points of the same subset, skip it.
  std::vector<size_t> ids(n, 0);
  std::vector<std::vector<size_t>> subsets(1);

  size_t id(0), a, minid, maxid;
  EdgeIndex index;
  for (const std::pair<EdgeIndex, double>& edge : edges) {

    // Find the indices of the vertices joined by the current edge
    index = edge.first;
    i = index.first;
    j = index.second;

    if ( !ids[i] && !ids[j] ) {
      // Set a new subset index for the two vertices
      id++;
      ids[i] = id;
      ids[j] = id;
      subsets.push_back({i, j});

    } else if ( ids[i] != ids[j] ) {
      // Merge the subsets
      minid = std::min(ids[i], ids[j]);
      maxid = std::max(ids[i], ids[j]);
      a = (ids[i] == maxid) ? j : i;
      if ( minid != 0 ) {
	      for (const size_t& l : subsets[minid]) {
          ids[l] = maxid;
          subsets[maxid].push_back(l);
        }
	      subsets[minid].resize(0);
      } else {
        ids[a] = maxid;
        subsets[maxid].push_back(a);
      }

    } else {
      continue;

    }

    // Add the edge to the MST
    _edges.push_back(index);
    _lengths.push_back(edge.second);
    _length += edge.second;
  }
}

std::vector<TLine*> MST::Lines(bool ndc) const {

  // Throw if the wrong dimension is requested
  std::vector<TLine*> lines;
  if ( _dim != 2 ) {
    throw(Exceptions::Exception(Exceptions::recoverable,
	  "Cannot return MST for dimension other than two",
	  "MST::Lines"));
    return lines;
  }

  // If the NDC system is requested, find the bounding box
  Vector<double> lower, upper, range(_dim);
  size_t i, j, d;
  if ( ndc ) {
    Math::BoundingBox(_points, lower, upper, 0.1);
    for (i = 0; i < lower.size(); i++)
        range[i] = upper[i]-lower[i];
  }

  // Loop over the edges, define a line per edge
  TLine* line;
  std::vector<double> a, b;
  for (const EdgeIndex& index : _edges) {

    // Find the indices of the vertices joined by the current edge
    i = index.first;
    j = index.second;

    // Get the points, convert to NDC if requested
    a = _points[i];
    b = _points[j];
    if ( ndc ) {
      for (d = 0; d < _dim; d++) {
        a[d] = (a[d]-lower[d])/range[d];
        b[d] = (b[d]-lower[d])/range[d];
      }
    }

    // Initialize the line, add it to the list
    line = new TLine(a[0], a[1], b[0], b[1]);
    lines.push_back(line);
  }

  return lines;
}

std::vector<TPolyLine3D*> MST::Lines3D() const {

  // Throw if the wrong dimension is requested
  std::vector<TPolyLine3D*> lines;
  if ( _dim != 3 ) {
    throw(Exceptions::Exception(Exceptions::recoverable,
	  "Cannot return MST for dimension other than three",
	  "MST::Lines3D"));
    return lines;
  }

  // Loop over the edges, define a line per edge
  TPolyLine3D* line;
  std::vector<double> a, b, x, y, z;
  size_t i, j;
  for (const EdgeIndex& index : _edges) {

    // Find the indices of the vertices joined by the current edge
    i = index.first;
    j = index.second;

    // Get the points
    b = _points[j];
    a = _points[i];

    // Initialize the line, add it to the list
    line = new TPolyLine3D(2);
    line->SetLineWidth(2);
    line->SetPoint(0, a[0], a[1], a[2]);
    line->SetPoint(1, b[0], b[1], b[2]);
    lines.push_back(line);
  }

  return lines;
}

void MST::Draw(const bool points) const {

  if ( _dim == 2 ) {

    // Draw the points if requested
    if ( points ) {
      TGraph* graph = new TGraph(_points.size());
      graph->SetTitle("");
      graph->SetMarkerStyle(20);
      for (size_t i = 0; i < _points.size(); i++)
          graph->SetPoint(i, _points[i][0], _points[i][1]);
      graph->Draw("AP SAME");
    }

    // Draw the MST
    for (TLine* line : Lines(!points))
        line->Draw("SAME");

  } else if ( _dim == 3 ) {

    // Find the bounding box
    Vector<double> lower, upper;
    Math::BoundingBox(_points, lower, upper, .1);

    // Create a 3d view
    TView *view = TView::CreateView(1);
    view->SetRange(lower[0], lower[1], lower[2], upper[0], upper[1], upper[2]);

    // Draw the points if requested
    if ( points ) {
      TPolyMarker3D *markers = new TPolyMarker3D();
      markers->SetMarkerStyle(20);
      for (size_t i = 0; i < _points.size(); i++)
	        markers->SetNextPoint(_points[i][0], _points[i][1], _points[i][2]);
      markers->Draw();
    }

    // Draw the MST
    for (TPolyLine3D* line : Lines3D())
        line->Draw("SAME");

  } else {
    throw(Exceptions::Exception(Exceptions::recoverable,
	  "Cannot draw MST of dimension "+std::to_string(_dim),
	  "MST::Draw"));
  }
}

void MST::Paint(const std::string name,
		            const bool points,
		            const std::vector<std::string> exts) const {

  TCanvas *c = new TCanvas("c", "c", 900, 900);

  if ( _dim == 2 || _dim == 3 ) {
    Draw(points);
  } else {
    throw(Exceptions::Exception(Exceptions::recoverable,
	  "Cannot paint MST of dimension "+std::to_string(_dim),
	  "MST::Paint"));
  }

  for (const std::string& ext : exts)
      c->SaveAs(std::string(name+"_MST."+ext).c_str());

  delete c;
}

double MST::Distance(size_t i, size_t j) const {

  double dist(0.);
  for (size_t k = 0; k < _dim; k++)
      dist += (_points[j][k]-_points[i][k])*(_points[j][k]-_points[i][k]);

  return sqrt(dist);
}
