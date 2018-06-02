#include "DensityEstimator.hh"

DensityEstimator::DensityEstimator() :
  _dim(), _name("de"), _algo(""), _points(), _ob(), _knn(),
  _lrd(), _tde(), _pbatde(), _prob(), _ex(false), _norm(1.) {
}

DensityEstimator::DensityEstimator(const std::vector<std::vector<double>>& points,
			   	   const std::string algo,
			   	   const bool ex,
			   	   const double norm) :
  _dim(), _name("de"), _algo(algo), _points(points), _ob(), _knn(),
  _lrd(), _tde(), _pbatde(), _prob(), _ex(ex), _norm(norm) {

  try {
    _name += "_"+_algo;
    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "DensityEstimator::DensityEstimator"));
  }
}

DensityEstimator::DensityEstimator(const DensityEstimator& de) {
  *this = de;
}

DensityEstimator& DensityEstimator::operator=(const DensityEstimator& de) {
  if ( this == &de )
      return *this;

  _dim = de._dim;
  _name = de._name;
  _algo = de._algo;
  _points = de._points;
  _ob = de._ob;
  _knn = de._knn;
  _lrd = de._lrd;
  _tde = de._tde;
  _pbatde = de._pbatde;
  _prob = de._prob;
  _ex = de._ex;
  _norm = de._norm;

  return *this;
}

DensityEstimator::~DensityEstimator () {}

void DensityEstimator::Initialize() {

  // Extract the dimension of the points
  Assert::IsNotEmpty("DensityEstimator::Initialize", "Input vector of point", _points);
  Assert::IsNonZero("DensityEstimator::Initialize", "Point dimension", _points[0].size());
  _dim = _points[0].size();  

  // Perform the density estimation on the provided points
  if ( _algo == "bin" ) {
    // Bin out the data to extract a grid of density
    _ob = OptimalBinning(_points, "scott");

  } else if ( _algo.find("knn") < _algo.size() ) {
    // Get the number of neighbours, use the rule of thumb k=n^(4/(4+d))
//    size_t k = (size_t)pow((double)_points.size(), 4./(4.+_dim));
//    size_t k = 1.53*(size_t)pow((double)_points.size(), 0.6);
    size_t k = 1.584*pow((double)_points.size(), 0.448);
    if ( _algo != "knn" ) {
      _algo.erase(_algo.find("knn"), 3);
      k = atoi(_algo.c_str());
    }

    // Feed the points to the k nearest neighbours algorithm
    _algo = "knn";
    _knn = KNearestNeighbours(_points, k, true);

  } else if ( _algo == "lrd" ) {
    // Get the number of neighbours, use the rule of thumb k=sqrt(N)
    size_t k = (size_t)pow((double)_points.size(), 1./2);

    // Feed the points to the local reachability density algorithm
    _lrd = LocalReachability(_points, k);

  } else if ( _algo == "tde" || _algo.find("cvt") < _algo.size() ) {
    // Get the PCVT factor
    double alpha = 1.;
    if ( _algo == "tde" ) {
      alpha = 0.;
    } else if ( _algo != "cvt" ) {
      std::string factor = _algo;
      factor.erase(factor.find("pcvt"), 4);
      alpha = atof(factor.c_str());
    }

    // Feed the points to a tesselation density estimation algorithm
    _algo = "tde";
    _tde = TDE(_points, false, true, alpha);

  } else if ( _algo.find("pbatde") < _algo.size() ) {
    // Get the content of each subset to produce, use the rule of thumb J=n^(d/(4+d))
//    size_t J = (size_t)pow((double)_points.size(), (double)_dim/(4.+_dim));
//    size_t J = 2.11028*pow((double)_points.size(), 1./3);
    size_t J = 3.567*pow((double)_points.size(), 0.357);
    if ( _algo != "pbatde" ) {
      _algo.erase(_algo.find("pbatde"), 6);
      J = atoi(_algo.c_str());
    }

    // Feed the points to the PBA tesselation density estimation algorithm
    _algo = "pbatde";
    _pbatde = PBATDE(_points, false, J);

  } else {
    // Density estimator not recognized, throw
    throw(Exceptions::Exception(Exceptions::recoverable,
	  "Algorithm not recognized: "+_algo,
	  "DensityEstimator::Initialize"));
  }
}

double DensityEstimator::Evaluate(const std::vector<double>& v) const {

  if ( _algo == "bin" ) {
    // In this case, use the interpolator
    return _norm*_ob(v);

  } else if ( _algo == "knn" ) {
    // Find the kNN density around the requested point
    return _norm*_knn(v);

  } else if ( _algo == "lrd" ) {
    // Find the LRD around the requested point
    return _lrd(v);

  } else if ( _algo == "tde" ) {
    // In this case, leave it to the TDE
    return _tde(v);

  } else if ( _algo == "pbatde" ) {
    // In this case, leave it to the PBATDE
    return _pbatde(v);

  }

  throw(Exceptions::Exception(Exceptions::nonRecoverable,
	"Algorithm not recognized: "+_algo,
	"DensityEstimator::Evaluate"));  
  return 0.;
}

double DensityEstimator::Evaluate(const double& v) const {

  return Evaluate(std::vector<double>({v}));
}

double DensityEstimator::Evaluate(const double* v) const {

  std::vector<double> vv(v, v+_dim);
  return Evaluate(vv);
}

double DensityEstimator::ContourVolume(double alpha,
				       const std::string method,
				       const bool draw,
				       std::vector<double> lower,
				       std::vector<double> upper,
				       const size_t nsamples) {

  // Number of points on the grid in each dimension for the drawing grid
  Assert::IsProbability("DensityEstimator::ContourVolume", "alpha", alpha);
  std::vector<size_t> n;
  for (size_t i = 0; i < _dim; i++)
      n.push_back(pow(1e4, 1./_dim));

  // Compute a bounding box to evaluate the contour in if no bounding box is provided
  double boxvol = 1.;
  if ( (lower.size() != _dim || upper.size() != _dim) && method != "mc" ) {
    Vector<double> l, u;
    boxvol = Math::BoundingBox(_points, l, u, .1);
    lower = l.std();
    upper = u.std();
  } else if ( method != "mc" ) {
    for (size_t i = 0; i < _dim; i++)
	boxvol *= (upper[i]-lower[i]);
  }

  // Fill what's necessary depending on the estimator
  if ( method == "box" ) {
    if ( !_prob.IsInitialized() ) {
      if ( _algo == "bin" ) {
        _prob = ProbabilityContour(GetGrid(n));
      } else {
        _prob = ProbabilityContour(GetGrid(n, lower, upper));
      }
    }

  } else if ( method == "mc" ) {
    std::vector<Vertex> points;
    if ( !_prob.IsInitialized() ) {

      // Evaluate at the input points if not already done
      size_t i;
      for (i = 0; i < _points.size(); i++)
	  points.push_back(Vertex(_points[i], Evaluate(_points[i])));

      // Sort the points by increasing density, select a fraction alpha of the points
      std::sort(points.begin(), points.end(), [] (const Vertex& a, const Vertex& b)
		{ return a.GetValue() < b.GetValue(); });
    } else {
      points = _prob.GetVertexArray();
    }

    // For a given fraction alpha, construct a bounding box that encompasses
    // a fraction alpha of the points only, this increases resolution
    size_t id = (1.-alpha)*(points.size()-1);
    std::vector<std::vector<double>> sub;
    size_t i, j;
    for (i = id; i < points.size(); i++)
	sub.push_back(points[i].GetCoordinates());

    Vector<double> l, u;	
    boxvol = Math::BoundingBox(sub, l, u, 1e-2);
    lower = l.std();
    upper = u.std();

    // Produce MC points uniformly distributed inside the orthotope and evaluate
    time_t tt;
    TRandom3 rdmzer(time(&tt));
    Vertex v(lower.size());
    std::vector<Vertex> mcvertices(nsamples);
    for (i = 0; i < nsamples; i++) {
      for (j = 0; j < lower.size(); j++)
	  v[j] = rdmzer.Uniform(lower[j], upper[j]);
      v.SetValue(this->Evaluate(v.GetCoordinates()));
      mcvertices[i] = v;
    }

    if ( !_prob.IsInitialized() ) {
      // Produce a grid if the contour must be drawn
      if ( draw ) {
         _prob = ProbabilityContour(points, mcvertices, boxvol, GetGrid(n, lower, upper));
      } else {
         _prob = ProbabilityContour(points, mcvertices, boxvol);
      }
    } else {
      _prob.SetMCVertexArray(mcvertices);
      _prob.SetBoxVolume(boxvol);
    }

    // Correct alpha if the parent is not completely contained within the range
    alpha /= _norm;

  } else if ( method == "hull" ) {
    if ( !_prob.IsInitialized() ) {
      // Evaluate at the input points
      size_t i;
      std::vector<Vertex> points;
      for (i = 0; i < _points.size(); i++)
	  points.push_back(Vertex(_points[i], Evaluate(_points[i])));

      // Produce a grid if the contour must be drawn
      if ( draw ) {
         _prob = ProbabilityContour(points, GetGrid(n, lower, upper));
      } else {
         _prob = ProbabilityContour(points);
      }
    }

    // Correct alpha if the parent is not completely contained within the range
    alpha /= _norm;

  } else {
    throw(Exceptions::Exception(Exceptions::recoverable,
	  "Reconstruction method not supported: "+method,
	  "DensityEstimator::ContourVolume"));
    return 0.;
  }

  _prob.Evaluate(alpha);
  return _prob.GetVolume();
}

double DensityEstimator::ContourVolumeError() {


  if ( _prob.IsInitialized() ) {
    try {
      // Compute the statistical uncertainty inherent to the density estimator
//      double alpha = _prob.GetAlpha()/_norm;
//      double serror = sqrt((1.-alpha)/(alpha*_prob.GetVertexArray().size()))
//				*exp(TMath::ChisquareQuantile(alpha, _dim)/4)*_prob.GetVolume();
//      double serror = 0.; // TODO TODO 

      // Add the uncertainty inherent to the volume reconstruction
      return /*serror+*/_prob.GetContourVolumeError();
    } catch ( Exceptions::Exception& e ) {
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	    "Could not compute the uncertainty on the contour volume"+std::string(e.what()),
	    "DensityEstimator::ContourVolumeError"));
    }
  }

  throw(Exceptions::Exception(Exceptions::recoverable,
	"The probability contour was not initialized",
	"DensityEstimator::ContourVolumeError"));
  return 0.;
}

Grid DensityEstimator::GetGrid(const std::vector<size_t>& number,
	       		       std::vector<double> lower,
	       		       std::vector<double> upper) const {

  // If lower or upper are not specified, get them from the grid. If no grid, throw
  Grid grid;
  if ( !lower.size() || !upper.size() ) {
    if ( _algo != "bin" ) {
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	    "This estimator ("+_algo+") is not based on a grid, must provide boundaries",
	    "DensityEstimator::GetGrid"));  
      return grid;
    }
    if ( !lower.size() )
	lower = _ob.GetLowerBoundArray();
    if ( !upper.size() )
	upper = _ob.GetUpperBoundArray();
  } else if ( number.size() != lower.size() || lower.size() != upper.size() ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Boundaries dimensions incompatible with the number of axes",
	  "DensityEstimator::GetGrid"));   
    return grid; 
  }

  // Create a new grid, loop over the points and set their values
  grid = Grid(number, lower, upper);
  size_t i;
  for (i = 0; i < grid.size(); i++)
      grid[i].SetValue(this->Evaluate(grid[i].GetCoordinates()));

  return grid;
}

std::vector<TPolyLine*> DensityEstimator::Meshing() const {

  if ( _algo == "bin" ) {
    return _ob.Meshing();

  } else if ( _algo == "tde" ) {
    return _tde.Meshing();

  } 

  throw(Exceptions::Exception(Exceptions::nonRecoverable,
	"No meshing produced for this density estimator: "+_algo,
	"DensityEstimator::Meshing"));  
  return std::vector<TPolyLine*>();
}

std::vector<TPolyLine*> DensityEstimator::Tesselation() const {

  if ( _algo == "tde" )
      return _tde.Polygons();

  throw(Exceptions::Exception(Exceptions::nonRecoverable,
	"No tesselation produced for this density estimator: "+_algo,
	"DensityEstimator::Tesselation"));  
  return std::vector<TPolyLine*>();
}

TGraph* DensityEstimator::Graph(double xmin, double xmax, size_t idx, std::vector<double> x) const {

  // Make sure the the requested id is smaller than the dimension
  Assert::IsGreater("DensityEstimator::Graph", "Dimension", _dim, idx, true);

  // Make sure the vector is of the right dimension. If not provided, fill with 0
  if ( !x.size() )
      x = std::vector<double>(_dim, 0);
  Assert::IsEqual("DensityEstimator::Graph", "Dimension and vector size", _dim, x.size());
  
  // For a certain amount of points, initialize a TGraph
  const size_t npoints = 100;
  TGraph *graph = new TGraph(npoints);
  graph->SetTitle("Density estimation");
  graph->GetXaxis()->SetTitle(TString::Format("x_{%d}", (int)idx));

  // Fill a TGraph with the estimator values
  size_t i;
  for (i = 0; i < npoints; i++) {
    x[idx] = xmin+i*(xmax-xmin)/(npoints-1);
    graph->SetPoint(i, x[idx], Evaluate(x));
  }
    
  // Set the style and return
  graph->SetLineWidth(2);
  graph->SetLineColor(2);
  return graph;
}

TH2F* DensityEstimator::Graph2D(double xmin, double xmax, double ymin, double ymax,
				size_t idx, size_t idy, std::vector<double> x) const {

  // Make sure the the requested ids are smaller than the dimension
  Assert::IsGreater("DensityEstimator::Graph2D", "Dimension", _dim, idx, true);
  Assert::IsGreater("DensityEstimator::Graph2D", "Dimension", _dim, idy, true);

  // Make sure the vector is of the right dimension. If not provided, fill with 0
  if ( !x.size() )
      x = std::vector<double>(_dim, 0);
  Assert::IsEqual("DensityEstimator::Graph2D", "Dimension and vector size", _dim, x.size());
  
  // For a certain amount of points, compute the estimator and fill a TH2F
  const size_t npoints = 100;
  TH2F *graph = new TH2F(TString::Format("%s_%d%d", _name.c_str(), (int)idx, (int)idy),
			 "Density estimation", npoints, xmin,
			 xmax, npoints, ymin, ymax);
  graph->GetXaxis()->SetTitle(TString::Format("x_{%d}", (int)idx));
  graph->GetYaxis()->SetTitle(TString::Format("x_{%d}", (int)idy));
  size_t i, j;
  for (i = 0; i < (size_t)graph->GetNbinsX(); i++) {
    for (j = 0; j < (size_t)graph->GetNbinsY(); j++) {
      x[idx] = graph->GetXaxis()->GetBinCenter(i+1);
      x[idy] = graph->GetYaxis()->GetBinCenter(j+1);
      graph->SetBinContent(i+1, j+1, Evaluate(x));
    }
  }

  // Return
  return graph;
}

void DensityEstimator::Draw(double xmin, double xmax, double ymin, double ymax, 
			    const std::string opt, int idx, int idy, std::vector<double> x) const {

  // Sort dimension wise and produce the relevant estimators
  if ( _dim == 1 || idy < 0 ) {
    TGraph* graph = Graph(xmin, xmax, idx, x);
    graph->Draw(opt.c_str());
    return;
  }
  
  TH2F* graph = Graph2D(xmin, xmax, ymin, ymax, idx, idy, x);
  graph->Draw(opt.c_str());
  return;
}

void DensityEstimator::Print(double xmin, double xmax, double ymin, double ymax,
			     const std::string opt, int idx, int idy, std::vector<double> x) const {

  try {
    gStyle->SetOptStat(0);
    TCanvas *c = new TCanvas("c", "c", 1200, 800);
    Draw(xmin, xmax, ymin, ymax, opt, idx, idy, x);
    c->SaveAs(TString::Format("%s.pdf", _name.c_str()));
    delete c;
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not print the estimator"+std::string(e.what()),
	  "DensityEstimator::Print"));
  }
}
