#include "AlphaComplex.hh"

AlphaComplex::AlphaComplex() :
  _dim(0), _alpha(0.), _points(), _facets(), _circumradii(), _areas(), _vol(0.) {
}

AlphaComplex::AlphaComplex(const std::vector<std::vector<double>>& points,
		           const double alpha) :
  _dim(0), _alpha(alpha), _points(points), _facets(), _circumradii(), _areas(), _vol(0.) {

  try {
    Assert::IsGreater("AlphaComplex::AlphaComplex", "Parameter alpha", alpha, 0.);
    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "AlphaComplex::AlphaComplex"));
  }
}

AlphaComplex::AlphaComplex(const AlphaComplex& alc) {
  *this = alc;
}

AlphaComplex& AlphaComplex::operator=(const AlphaComplex& alc) {
  if ( this == &alc )
      return *this;

  _dim = alc._dim;
  _alpha = alc._alpha;
  _points = alc._points;
  _facets = alc._facets;
  _circumradii = alc._circumradii;
  _areas = alc._areas;
  _vol = alc._vol;

  return *this;
}

AlphaComplex::~AlphaComplex () {}

void AlphaComplex::Initialize() {

  // Check that there are vertices to fit and there are at least n+1 vertices for dimension n
  Assert::IsNotEmpty("AlphaComplex::Initialize", "Vector of input points", _points);
  Assert::IsNonZero("AlphaComplex::Initialize", "Dimension", _points[0].size());
  if ( _points[0].size() > _points.size()-1 )
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	    "Not enough points to triangulate, need at least dim+1",
	    "AlphaComplex::Initialize"));

  // Fill a vector with the right structure for the Qhull code
  size_t N = _points.size();  // Number of vertices
  _dim = _points[0].size();	// Dimension of the space
  std::vector<double> points;
  size_t i, j;
  for (i = 0; i < N; i++)
      for (j = 0; j < _dim; j++)
          points.push_back(_points[i][j]);

  // Feed the vector to Qhull and compute the Delaunay triangulation ("d")
  // Must use a pointer, the Qhull object owns the memory of its subcomponents, cannot be copied
  // TODO, the case of N=n+1 must be handled differently (not enough points for a Qhull!)
  orgQhull::Qhull qhull("", _dim, N, &(points[0]), "d");

  // Store the facets for fast access. Each facet is an array of n+1 n-points (simplex)
  std::vector<std::vector<double>> tempfacet(_dim+1);
  for (orgQhull::QhullFacet& facet : qhull.facetList().toStdVector()) {

    // Compute and store the circumradii of all the facets. This is done so that the fraction
    // alpha can be easily changed without having to recompute the Delaunay triangulation.
    // If alpha=0, this is equivalent to the convexe hull
    for (i = 0; i < _dim+1; i++)
	tempfacet[i] = _points[facet.vertices()[i].point().id()];

    _facets.push_back(tempfacet);
    _circumradii.push_back(Circumradius(facet));
    _areas.push_back(facet.facetArea());

    // Increment the volume of the complex if the circumradii is within range
    if ( _alpha )
      if ( _circumradii.back() > 1./_alpha )
	  continue;

    _vol += _areas.back();
  }
}

void AlphaComplex::SetAlpha(const double alpha) {

  // Reset alpha
  Assert::IsGreater("AlphaComplex::SetAlpha", "Parameter alpha", alpha, 0.);
  _alpha = alpha;

  // Recompute the volume
  _vol = 0;
  size_t i;
  for (i = 0; i < _facets.size(); i++) {
    if ( _alpha )
      if ( _circumradii[i] > 1./_alpha )
	  continue;
    _vol += _areas[i];
  }
}

std::vector<TPolyLine*> AlphaComplex::Polygons(const bool fill) const {

  // Throw if the wrong dimension is requested
  std::vector<TPolyLine*> polylines;
  if ( _dim != 2 ) {
    throw(Exceptions::Exception(Exceptions::recoverable,
	  "Cannot return alpha-complex polygons for dimensions other than 2",
	  "AlphaComplex::Polygons"));
    return polylines;
  }

  // Loop over the facets, define a polygon per facet
  std::vector<double> x, y;
  size_t i;
  for (i = 0; i < _facets.size(); i++) {
    // If the circumradius of this facet is too large, skip it
    if ( _alpha )
      if ( _circumradii[i] > 1./_alpha )
	  continue;

    // Construct the vertices of the TPolyLine
    x.resize(0);
    y.resize(0);
    for (const std::vector<double>& vertex : _facets[i]) {
      x.push_back(vertex[0]);
      y.push_back(vertex[1]);
    }
    x.push_back(x.front()); // Close the shape
    y.push_back(y.front()); // Close the shape

    polylines.push_back(new TPolyLine(4, &(x[0]), &(y[0])));

    // If fill is requested, fill the facets and remove the edges
    if ( fill ) {
      polylines.back()->SetLineWidth(0);
      polylines.back()->SetFillColorAlpha(1, .25);
    } else {
      polylines.back()->SetLineWidth(.1);
    }
  }

  return polylines;
}

TGeoVolume* AlphaComplex::Polyhedra() const {

  // Throw if the wrong dimension is requested
  if ( _dim != 3 ) {
    throw(Exceptions::Exception(Exceptions::recoverable,
	  "Cannot return alpha-complex polyhedra for dimensions other than 3",
	  "AlphaComplex::Polyhedra"));
    return new TGeoVolume();
  }

  // Find the boundaries of the vertices to get the range of the top volume box
  Vector<double> lower, upper;
  Math::BoundingBox(_points, lower, upper);
  Vector<double> pos = (lower+upper)/2;
  Vector<double> hl = (upper-lower)/2; // Half length of the interval

  // Set up the box that contains all of the simplices
  TGeoBBox *box = new TGeoBBox(hl[0], hl[1], hl[2], &pos[0]);
  TGeoVolume *top = new TGeoVolume("top", box);
  gGeoManager->SetTopVolume(top);

  // Loop over the facets and add a simplex for each of them
  TGeoVolume* temp(NULL);
  TGeoHMatrix* H(NULL);
  size_t i;
  for (i = 0; i < _facets.size(); i++) {
    // If the circumradius of this facet is too large, skip it
    if ( _alpha )
      if ( _circumradii[i] > 1./_alpha )
	  continue;

    // Add the facet to the list of volumes
    temp = Polyhedron(_facets[i], H);
    top->AddNode(temp, i, H);
  }

  return top;
}

void AlphaComplex::Draw(const bool fill,
			const bool points) const {

  if ( _dim == 2 ) {
    if ( points ) {
      TGraph* graph = new TGraph(_points.size());
      graph->SetTitle("");
      graph->SetMarkerStyle(20);
      graph->SetMarkerSize(.5);
      for (size_t i = 0; i < _points.size(); i++)
	  graph->SetPoint(i, _points[i][0], _points[i][1]);
      graph->Draw("AP SAME");
    } 

    for (TPolyLine* poly : Polygons(fill))
      if ( fill ) {
	  poly->Draw("F SAME");
      } else {
	  poly->Draw("SAME");
     }

  } else if ( _dim == 3 ) {
    gGeoManager = new TGeoManager("polyhedra", "3D Alpha-Complex");
    gGeoManager->SetVerboseLevel(0);
    gGeoManager->SetMaxVisNodes(1e5);
    if ( !fill ) {
        Polyhedra()->Draw("SAME");
        TView *view = gPad->GetView();
        view->ShowAxis();
	view->TopView();
    } else {
        Polyhedra()->Raytrace();
    }

    if ( points ) {
      TPolyMarker3D *marker = new TPolyMarker3D();
      marker->SetMarkerStyle(20);
      marker->SetMarkerSize(.5);
      for (size_t i = 0; i < _points.size(); i++)
	  marker->SetNextPoint(_points[i][0], _points[i][1], _points[i][2]);
      marker->Draw("SAME");
    }

  } else {
    throw(Exceptions::Exception(Exceptions::recoverable,
	  "Cannot draw alpha-complex in dimension "+std::to_string(_dim),
	  "AlphaComplex::Draw"));
  }
}

void AlphaComplex::Paint(const std::string name,
			 const bool fill,
			 const bool points,
			 const std::vector<std::string> exts) const {

  TCanvas *c = new TCanvas("c", "c", 900, 900);
  if ( _dim == 2 ) {
    if ( points ) {
      Draw(fill, points);
    } else {
      Vector<double> lower, upper;
      Math::BoundingBox(_points, lower, upper, .1);
      c->DrawFrame(lower[0], lower[1], upper[0], upper[1]);
      Draw(fill);
    }

  } else if ( _dim == 3 ) {
    Draw(fill, points);

  } else {
    throw(Exceptions::Exception(Exceptions::recoverable,
	  "Cannot paint alpha-complex in dimension "+std::to_string(_dim),
	  "AlphaComplex::Paint"));
  }

  for (const std::string& ext : exts)
      c->SaveAs(std::string(name+"_AC."+ext).c_str());

  delete c;
}

double AlphaComplex::Circumradius(orgQhull::QhullFacet& facet) const {

  // First need to form the Euclidean distance matrix (EDM)
  Matrix<double> EDM(_dim+1, _dim+1);
  size_t i, j;
  Vector<double> vi, vj;
  for (i = 0; i < _dim+1; i++) {
    for (j = 0; j < _dim+1; j++) {
      if ( j <= i ) 
	  EDM[i][j] = EDM[j][i];

      vi = Vector<double>(facet.vertices()[i].point().coordinates(), _dim);
      vj = Vector<double>(facet.vertices()[j].point().coordinates(), _dim);
      EDM[i][j] = (vj-vi).norm2();
    }
  }

  // Need to form the Cayley-Menger (n+2)*(n+2) matrix of the simplex
  Matrix<double> CMM(_dim+2, _dim+2);
  CMM[0][0] = 0.;
  for (i = 0; i < _dim+1; i++) {
    CMM[i+1][0] = 1.;
    CMM[0][i+1] = 1.;
  }
  for (i = 0; i < _dim+1; i++)
    for (j = 0; j < _dim+1; j++)
        CMM[i+1][j+1] = EDM[i][j];

  return sqrt(EDM.Determinant()/CMM.Determinant()/(-2.));
}

TGeoVolume* AlphaComplex::Polyhedron(const std::vector<std::vector<double>>& facet,
				     TGeoHMatrix*& H) const {

  // Convert the four simplex vertices to Vectors
  std::vector<Vector<double>> V(4);
  size_t i;
  for (i = 0; i < 4; i++)
      V[i] = Vector<double>(facet[i]);

  // To draw the simplex, one first needs to put oneself into the referential of the simplex.
  // Choose the unit-vector from p0 to p1 as ex', the normal as ez' and ey' as ez' cross ex' (RH).
  // The origin of the new reference frame is simply p0
  std::vector<double> pars = Math::PlaneParameters(V);
  Vector<double> v01 = V[1]-V[0];
  double n01 = v01.norm();

  Vector<double> ez(&pars[0], 3);
  Vector<double> ex = v01/n01;
  Vector<double> ey = ez.cross(ex);

  // Create a rotation matrix from the standard ex*ey*ez to the simplex basis
  Matrix<double> R(3, 3);
  R[0] = ex.std();
  R[1] = ey.std();
  R[2] = ez.std();

  // Rotate the vectors into the simplex basis: p' = R*(p-p_0). Sort the first 3 clockwise
  std::vector<Vector<double>> Vprime(3);
  for (i = 0; i < 3; i++)
      Vprime[i] = R*(V[i]-V[0]);
  Math::SortCW(Vprime);
  Vprime.push_back(R*(V[3]-V[0]));

  // The fourth point is the only one that is not coplanar, find the half-distance
  double dz = Vprime[3][2]/2;

  // Fill an array of points with the V', fill the base and the peak separately
  std::vector<double> base(8), peak(8);
  for (i = 0; i < 4; i++) {
    base[2*i] = (i == 3) ? Vprime[0][0] : Vprime[i][0];
    base[2*i+1] = (i == 3) ? Vprime[0][1] : Vprime[i][1];
    peak[2*i] = Vprime[3][0];
    peak[2*i+1] = Vprime[3][1];
  }

  // If dz is negative, the peak comes first, base otherwise
  std::vector<double> points;
  if ( dz < 0 ) {
    points.insert(points.end(), peak.begin(), peak.end());
    points.insert(points.end(), base.begin(), base.end());
  } else {
    points.insert(points.end(), base.begin(), base.end());
    points.insert(points.end(), peak.begin(), peak.end());
  }

  // Initialize the tetrahedron
  TGeoArb8 *tetra = new TGeoArb8(fabs(dz), &points[0]);

  // Set the volume in the shape of a tetrahedron and return
  TGeoVolume *simplex = new TGeoVolume("tetra", tetra);
  simplex->SetLineWidth(1);

  // Compute the rotation Euler angles of the inverse rotation matrix (bring them back)
  Matrix<double> Rinv = R.Transpose();
  double thx = atan2(Rinv[2][1], Rinv[2][2])*360/(2*M_PI);
  double thy = atan2(-Rinv[2][0], sqrt(pow(Rinv[2][1], 2)+pow(Rinv[2][2], 2)))*360/(2*M_PI);
  double thz = atan2(Rinv[1][0], Rinv[0][0])*360/(2*M_PI);

  // Set the rotation of the volume into the frame of the parent
  TGeoTranslation Tz(0., 0., dz);
  TGeoRotation Rxyz;
  Rxyz.RotateX(thx);
  Rxyz.RotateY(thy);
  Rxyz.RotateZ(thz);
  TGeoTranslation T(V[0][0], V[0][1], V[0][2]);
  TGeoHMatrix mat = T*Rxyz*Tz;
  H = new TGeoHMatrix(mat);

  return simplex;
}
