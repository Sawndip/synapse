#include "Voronoi.hh"

Voronoi::Voronoi() :
  _dim(0), _points(), _cells(), _bounded(false), _lower(), _upper(), _kdtree(NULL)  {
}

Voronoi::Voronoi(const std::vector<std::vector<double>>& points,
	  	 const bool bounded,
	  	 const double alpha,
	  	 const double eps) :
  _dim(), _points(points), _cells(), _bounded(bounded), _lower(), _upper(), _kdtree(NULL) {

  try {
    Initialize(alpha, eps);
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "Voronoi::Voronoi"));
  }
}

Voronoi::Voronoi(const Voronoi& vor) {
  *this = vor;
}

Voronoi& Voronoi::operator=(const Voronoi& vor) {
  if ( this == &vor )
      return *this;

  _dim = vor._dim;
  _points = vor._points;
  _cells = vor._cells;
  _bounded = vor._bounded;
  _lower = vor._lower;
  _upper = vor._upper;

  // Special care has to be taken here as NanoFLANN does not have a copy constructor
  if ( _kdtree )
      delete _kdtree;
  if ( _dim && _points.size() ) {
    _kdtree = new KDTree(_dim, _points, 10);
    _kdtree->index->buildIndex();
  }

  return *this;
}

Voronoi::~Voronoi () {

  delete _kdtree;
}

void Voronoi::Initialize(const double alpha, const double eps) {

  // Check that there are points to fit and there are at least n+1 points for dimension n
  Assert::IsNotEmpty("Voronoi::Initialize", "Vector of input points", _points);
  Assert::IsNonZero("Voronoi::Initialize", "Dimension", _points[0].size());
  size_t N = _points.size();
  _dim = _points[0].size();
  if ( !_bounded && _dim > N-1 )
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	    "Not enough points to tesselate, need at least dim+1",
	    "Voronoi::Initialize"));

  // Check the value of alpha
  Assert::IsProbability("Voronoi::Initialize", "PCVT factor", alpha, false, false);

  // Prepare a vector to harbour the voronoi cells, one per point
  _cells.resize(N);

  // If requested, find the boundaries of the n-orthotope that encompasses all points.
  // Must offset the boundaries to make sure that there are no points on the boundaries
  // which would result in zero-volume cells.
  size_t i, j, k, l;
  double tmin, tmax;
  double frac = .1;
  if ( _bounded ) {
    for (i = 0; i < _dim; i++) {
      // Find the extremums in this dimension
      tmin = DBL_MAX;
      tmax = -DBL_MAX;
      for (j = 0; j < N; j++) {
        if ( _points[j][i] < tmin )
	    tmin = _points[j][i];
        if ( _points[j][i] > tmax )
	    tmax = _points[j][i];
      }

      _lower.push_back(tmin-(tmax-tmin)*frac);
      _upper.push_back(tmax+(tmax-tmin)*frac);
    }
  }

  // Treat the 1D tessellation separately as it is not handled by QHull
  if ( _dim == 1 ) {
    Initialize1D(alpha, eps);
    return;
  }

  // Initialize the points vector to be fed to QHull and the default polytope centres
  std::vector<Vector<double>> centres(N);
  std::vector<double> points;
  for (j = 0; j < N; j++) {
    centres[j] = Vector<double>(_points[j]);
    for (k = 0; k < _dim; k++)
	points.push_back(_points[j][k]);
  }

  // Loop as many time as requested to compute the Penalised Centroidal Voronoi Tessellation (PCVT)
  // At each additional iteration, use the centroids of each cell as a new point for
  // the next iteration (Lloyd's algorithm). Use a sum of the distance between the points and the
  // cell centroids as well as the distance between current and orignal point as a penalty function
  orgQhull::Qhull *qhull(NULL), *qhull_temp(NULL);	// QHull owns the memory
  std::vector<size_t> edge_ids;				// IDs of the points on the hull
  size_t id;						// ID of the current point
  size_t npoints;					// Number of points on the polytope
  Vector<double> first(_dim), current(_dim);		// First and current point position
  double penalty(FLT_MAX), temp(DBL_MAX);
  size_t ite = 0;
  while ( fabs(penalty-temp) > eps ) {

    // Reinitialize the penatly function
    temp = penalty;
    penalty = 0.;

    // If a bounded Voronoi is requested, produce (2*_dim) more points for each of the N points.
    // We end up with (2*_dim+1) times the initial amount of points
    if ( _bounded ) {
      points.resize(_dim*N);
      for (j = 0; j < N; j++) {
        for (k = 0; k < 2*_dim; k++) {

	  // Copy the corresponding point, only one of the coordinates is reflected
          for (l = 0; l < _dim; l++)
	      points.push_back(points[j*_dim+l]);

	  // For odd numbers, add a reflexion point about the max, otherwise about the min
	  if ( k % 2 ) {
	    points[_dim*N+j*2*_dim*_dim+k*_dim+k/2] = 2*_upper[k/2]-points[j*_dim+k/2];
	  } else {
	    points[_dim*N+j*2*_dim*_dim+k*_dim+k/2] = 2*_lower[k/2]-points[j*_dim+k/2];
	  }
        }
      }
    }

    // Feed the vector to Qhull and compute the Voronoi Tessellation ("v")
    // Must use a pointer, the Qhull object owns the memory of its subcomponents, cannot be copied
    qhull_temp = new orgQhull::Qhull("", _dim, points.size()/_dim, &(points[0]), "v");

    // Compute a standard convex hull and store a list of the point IDs that lie on the hull
    // It is only necessary if the points are not bounded, only once as they don't change
    if ( !_bounded ) {
      orgQhull::Qhull qhull_edge("", _dim, N, &(points[0]), "");
      edge_ids.resize(0);
      for (const orgQhull::QhullVertex& vertex : qhull_edge.vertexList().toStdVector())
          edge_ids.push_back(vertex.point().id());
    }

    // Loop over the points, compute the centroids of their corresponding Voronoi cells
    for (const orgQhull::QhullVertex& vertex : qhull_temp->vertexList().toStdVector()) {

      // If the points are bounded, skip the construction points
      id = vertex.point().id();
      if ( _bounded && id >= N)
	  continue;

      // Test weather the vertex has an open voronoi cell or not, skip open cells
      centres[id] = Vector<double>(vertex.point().coordinates(), _dim);
      if ( !_bounded && std::find(edge_ids.begin(), edge_ids.end(), id) != edge_ids.end() )
	  continue;

      if ( alpha ) {
        // Get the position of this vertex and its original position
	first = Vector<double>(_points[id]);
	current = Vector<double>(vertex.point().coordinates(), _dim);

        // Loop over the facets adjacent to the vertex and compute their circumcentres barycentre
	npoints = vertex.neighborFacets().size();
        centres[id] = Vector<double>(_dim, 0);
        for (orgQhull::QhullFacet& facet : vertex.neighborFacets().toStdVector())
          for (k = 0; k < _dim; k++)
	      centres[id][k] += facet.voronoiVertex()[k]/npoints;

        // Increment the penalty function
        penalty += alpha*pow((centres[id]-current).norm(), 2);
        penalty += (1.-alpha)*pow((current-first).norm(), 2);
      } else {
	penalty = temp;
      }
    }

    // If the penalty function has decreased or the hull pointer hasn't been set, set it
    if ( !qhull ) {
       qhull = qhull_temp;
    } else if ( penalty < temp ) {
       delete qhull;
       qhull = qhull_temp;
    }

    // If no convergence can be reached, break
    ite ++;
    if ( ite > 1000 )
        break;

    // If the penalty function has decreased, update the points to the COMs
    for (j = 0; j < N; j++)
      for (k = 0; k < _dim; k++)
	  points[j*_dim+k] = alpha*centres[j][k]+(1.-alpha)*_points[j][k];
  }

  // Loop over the points, fill the containers with data
  std::vector<double> cellpoints;	// Array of in-line cell points
  Hull cell;				// Hull object for each Voronoi cell
  double vol;				// Volume of the Voronoi cell
  for (const orgQhull::QhullVertex& vertex : qhull->vertexList().toStdVector()) {

    // If the points are bounded, skip the construction points
    id = vertex.point().id();
    if ( _bounded && id >= N)
        continue;

    // Loop over the facets adjacent to the vertex and yield their circumcentres
    cellpoints.resize(0);
    for (orgQhull::QhullFacet& facet : vertex.neighborFacets().toStdVector())
      for (k = 0; k < _dim; k++)
          cellpoints.push_back(facet.voronoiVertex()[k]);

    // Fill the Voronoi cells, if the cell is on the edge, define the volume as -1 (infinite)
    // The volume of the cell is the convexe hull volume. If the cell is a simplex, bypass Qhull.
    vol = -1;
    if ( _bounded || std::find(edge_ids.begin(), edge_ids.end(), id) == edge_ids.end() ) {
      if ( cellpoints.size()/_dim == _dim+1 ) {
	vol = SimplexVolume(cellpoints);
      } else {
        orgQhull::Qhull subhull("", _dim, cellpoints.size()/_dim, &(cellpoints[0]), "");
        vol = subhull.volume();
      }
    }

    // Initialize the cell
    cell.SetPointArray(cellpoints, _dim);
    cell.SetVolume(vol);
    _cells[id] = cell;

    // update the points if necessary
    if ( alpha )
      for (k = 0; k < _dim; k++)
	  _points[id][k] = vertex.point().coordinates()[k];
 }

  // Clear the memory
  delete qhull;
//  this->Draw("test", true, false, true);
}

void Voronoi::Initialize1D(const double alpha, const double eps) {

  // Create an array that contains the coordinate of the point and retains its id
  std::vector<std::pair<size_t, double>> points;
  size_t j;
  for (j = 0; j < _points.size(); j++)
      points.push_back(std::pair<size_t, double>(j, _points[j][0]));

  // Sort the points from lowest to highest
  size_t N = points.size();
  std::sort(points.begin(), points.end(),
	    [] (const std::pair<size_t, double>& a, const std::pair<size_t, double>& b) {
		return a.second < b.second;
	    });

  // Optimize the penalty function if a PCVT was requested
  if ( alpha ) {
    std::vector<double> centres(N);
    double penalty(FLT_MAX), temp(DBL_MAX);
    size_t ite = 0;
    while ( fabs(penalty-temp)/temp > eps ) {

      temp = penalty;
      penalty = 0.;
      for (j = 0; j < N; j++) {

        // Compute the centroid of each voronoi cell, compare it with the previous and initial
        if ( !j ) {
          if ( _bounded ) {
	    centres[j] = (2*_lower[0]+points[0].second+points[1].second)/4;
	  } else {
	    centres[j] = points[0].second;
	  }
        } else if ( j == N-1 ) {
          if ( _bounded ) {
	    centres[j] = (2*_upper[0]+points[N-1].second+points[N-2].second)/4;
	  } else {
	    centres[j] = points[N-1].second;
	  }
        } else {
	  centres[j] = (points[j-1].second+2*points[j].second+points[j+1].second)/4;
        }

        // Increment the penalty function if necessary
        penalty += alpha*pow(points[j].second-centres[j], 2);
        penalty += (1.-alpha)*pow(points[j].second-_points[points[j].first][0], 2);
      }

    // If no convergence can be reached, break
    ite ++;
    if ( ite > 1000 )
        break;

    // Update the points to the COMs if necessary
    if ( alpha )
      for (j = 0; j < N; j++)
	  points[j].second = alpha*centres[j]+(1.-alpha)*_points[points[j].first][0];
    }
  }

  // Fill the cells and volumes with the (updated) points
  std::vector<std::vector<double>> cellpoints;
  Hull cell;
  double vol;
  for (j = 0; j < N; j++) {

    // Fill the limits of each cell and their volumes
    // The volume of a cell is V_i=(x_{i+1}-x_{i-1})/2
    cellpoints.resize(0);
    if ( !j ) {
      if ( !_bounded ) {
	vol = -1.;
      } else {
	cellpoints.push_back({_lower[0]});
	vol = (points[0].second+points[1].second)/2-_lower[0];
      }
      cellpoints.push_back({(points[0].second+points[1].second)/2});
    } else if ( j == N-1 ) {
      if ( !_bounded ) {
	vol = -1.;
      } else {
	cellpoints.push_back({_upper[0]});
	vol = _upper[0]-(points[N-1].second+points[N-2].second)/2;
      }
      cellpoints.push_back({(points[N-1].second+points[N-2].second)/2});
    } else {
      vol = (points[j+1].second-points[j-1].second)/2;
      cellpoints.push_back({(points[j].second+points[j-1].second)/2});
      cellpoints.push_back({(points[j].second+points[j+1].second)/2});
    }

    // Fill the cell and add it to the list
    cell.SetPointArray(cellpoints);
    cell.SetVolume(vol);
    _cells[points[j].first] = cell;

    // Update the points
    _points[points[j].first] = {points[j].second};
  }

  // Initialize the k-d tree
  _kdtree = new KDTree(_dim, _points, 10);
  _kdtree->index->buildIndex();
}

size_t Voronoi::GetBestCellID(const std::vector<double>& v) const {

  // Look into the indexed points of the k-d tree to find the closest point
  std::vector<size_t> indices(1);
  std::vector<double> dists(1);
  nanoflann::KNNResultSet<double> results(1);
  results.init(&indices[0], &dists[0]);
  _kdtree->index->findNeighbors(results, &v[0], nanoflann::SearchParams(10));
  return indices[0];
}

const Hull& Voronoi::GetBestCell(const std::vector<double>& v) const {

  // Fill the vertex vector with the right facet and its mapping
  return _cells[GetBestCellID(v)];
}

double Voronoi::SimplexVolume(const std::vector<double>& points) const {

  // Check that it is indeed a simplex in this space
  if ( points.size()/_dim != _dim+1 )
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	    "A d-simplex mst have exactly d+1 points",
	    "Voronoi::Initialize"));

  // Set the matrix of vector differences
  Matrix<double> M(_dim, _dim);
  size_t i, j;
  for (i = 0; i < _dim; i++)
    for (j = 0; j < _dim; j++)
	M[i][j] = points[(i+1)*_dim+j]-points[j];

  // Compute the factorial of the dimension and return
  double fact = 1.;
  for (i = 2; i < _dim; i++)
      fact *= i;
  return fabs(M.Determinant()/fact);
}

TPolyLine* Voronoi::BoundingBox() const {

  // Throw if the wrong dimension is requested
  TPolyLine* rectangle = NULL;
  Assert::IsEqual("Voronoi::BoundingBox", "Dimension", _dim, (size_t)2);

  // Create a regtangle with the limits of the bounding box
  std::vector<std::vector<int>> keys = {{0, 0}, {1, 0}, {1, 1}, {0, 1}, {0, 0}};
  std::vector<std::vector<double>> points(2);
  size_t i, j;
  for (i = 0; i < keys.size(); i++) {
    for (j = 0; j < _dim; j++) {
      if ( keys[i][j] ) {
	  points[j].push_back(_upper[j]);
      } else {
	  points[j].push_back(_lower[j]);
      }
    }
  }

  rectangle = new TPolyLine(keys.size(), &(points[0][0]), &(points[1][0]));
  return rectangle;
}

std::vector<TPolyLine*> Voronoi::Polygons() const {

  // Throw if the wrong dimension is requested
  std::vector<TPolyLine*> polylines;
  Assert::IsEqual("Voronoi::Polygons", "Dimension", _dim, (size_t)2);

  // Loop over the Voronoi cells, get a polygon per cell
  for (const Hull& cell : _cells)
    if ( _bounded || cell.GetVolume() > 0 )
      	polylines.push_back(cell.Polygon());

  return polylines;
}

TH2Poly* Voronoi::DensityProfile() const {

  // Throw if the wrong dimension is requested
  TH2Poly* profile(NULL);
  Assert::IsEqual("Voronoi::DensityProfile", "Dimension", _dim, (size_t)2);
  if ( _bounded ) {
    profile = new TH2Poly("", "", _lower[0], _upper[0], _lower[1], _upper[1]);
  } else {
    // If the tesselation is not bounded, find the minimum and maximum
    size_t i, j;
    double tmin, tmax;
    double frac = .1;
    std::vector<double> lower, upper;
    for (i = 0; i < _dim; i++) {
      // Find the extremums in this dimension
      tmin = DBL_MAX;
      tmax = -DBL_MAX;
      for (j = 0; j < _points.size(); j++) {
        if ( _points[j][i] < tmin )
	    tmin = _points[j][i];
        if ( _points[j][i] > tmax )
	    tmax = _points[j][i];
      }

      lower.push_back(tmin-(tmax-tmin)*frac);
      upper.push_back(tmax+(tmax-tmin)*frac);
    }
    profile = new TH2Poly("", "", lower[0], upper[0], lower[1], upper[1]);
  }

  // Loop over the Voronoi cells, add a bin for each one
  Matrix<double> cmfv;
  double vol;
  size_t bin;
  for (const Hull& cell : _cells) {
    vol = cell.GetVolume();
    if ( vol  > 0 ) {
      cmfv = cell.GetCWPoints();
      bin = profile->AddBin(cmfv.Ncols(), &(cmfv[0][0]), &(cmfv[1][0]));
      profile->SetBinContent(bin, 1./vol/_points.size());
    }
  }

  return profile;
}

void Voronoi::Draw(std::string name, bool color, bool points, bool axes) const {

  // Initialize the density profile and the canvas
  TH2Poly* profile = DensityProfile();
  if ( !profile )
      return;
//  profile->GetXaxis()->SetTitle("x [mm]");
//  profile->GetYaxis()->SetTitle("p_{x} [MeV/c]");
  gStyle->SetOptStat(0);
  TCanvas* c = new TCanvas("c", "c", 900, 900);

  gPad->SetRightMargin(.15);
  gPad->SetLogz();
  profile->GetXaxis()->SetTitle("x");
  profile->GetYaxis()->SetTitle("y");
  profile->GetZaxis()->SetTitle("#rho(x,y)");

  // If no axes are requested, remove the margins and get rid of the ticks
  if ( !axes ) {
    c->SetMargin(0., 0., 0., 0.);
    profile->GetXaxis()->SetTickLength(0);
    profile->GetXaxis()->SetLabelOffset(999);
    profile->GetYaxis()->SetTickLength(0);
    profile->GetYaxis()->SetLabelOffset(999);
    profile->GetZaxis()->SetTickLength(0);
    profile->GetZaxis()->SetLabelOffset(999);
  }

  // Depending on the color switch, draw the histogram with or without colors
  if ( color ) {
    profile->Draw("COLZ");
  } else {
    profile->Draw("L");
  }

  // If the points are requested, fill a TGraph and draw it
  TGraph* graph = new TGraph();
  if ( points ) {
    graph->SetMarkerStyle(20);
    size_t i;
    for (i = 0; i < _points.size(); i++)
	graph->SetPoint(i, _points[i][0], _points[i][1]);
    graph->Draw("PSAME");
  }

  c->SaveAs(std::string("VT_"+name+".pdf").c_str());
  delete c;
  delete profile;
  delete graph;
}
