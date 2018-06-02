#include "ProbabilityContour.hh"

ProbabilityContour::ProbabilityContour() :
  _algo(""), _grid(), _indices(0), _alpha(0), _volume(0), _level(0),
  _vertices(0), _mcvertices(0), _boxvol(0.) {
}

ProbabilityContour::ProbabilityContour(const std::vector<Vertex>& vertices,
				       const Grid grid,
				       const std::string algo) :
  _algo(algo), _grid(grid), _indices(0), _alpha(0), _volume(0), _level(0),
  _vertices(vertices), _mcvertices(0), _boxvol(0.) {

  // Sort the vertices in ascending order of density
  std::sort(_vertices.begin(), _vertices.end(), [] (const Vertex& a, const Vertex& b) {
		return a.GetValue() < b.GetValue();
	    });
}

ProbabilityContour::ProbabilityContour(const Grid& grid,
				       const std::string algo) :
  _algo(algo), _grid(grid), _indices(0), _alpha(0), _volume(0), _level(0),
  _vertices(grid.GetVertexArray()), _mcvertices(0), _boxvol(0.) {

  // Sort the vertices in ascending order of density
  std::sort(_vertices.begin(), _vertices.end(), [] (const Vertex& a, const Vertex& b) {
		return a.GetValue() < b.GetValue();
	    });
}

ProbabilityContour::ProbabilityContour(const std::vector<Vertex>& vertices,
		     		       const std::vector<Vertex>& mcvertices,
		     		       const double& boxvol,
				       const Grid grid) :
  _algo("mc"), _grid(grid), _indices(0), _alpha(0), _volume(0),
  _vertices(vertices), _mcvertices(mcvertices), _boxvol(boxvol) {

  // Sort the vertices and MC vertices in ascending order of density
  std::sort(_vertices.begin(), _vertices.end(), [] (const Vertex& a, const Vertex& b) {
		return a.GetValue() < b.GetValue();
	    });
  std::sort(_mcvertices.begin(), _mcvertices.end(), [] (const Vertex& a, const Vertex& b) {
		return a.GetValue() < b.GetValue();
	    });
}

ProbabilityContour::ProbabilityContour(const ProbabilityContour& cont) {
  *this = cont;
}

ProbabilityContour& ProbabilityContour::operator=(const ProbabilityContour& cont) {
  if ( this == &cont )
      return *this;

  _algo = cont._algo;
  _grid = cont._grid;
  _indices = cont._indices;
  _alpha = cont._alpha;
  _volume = cont._volume;
  _level = cont._level;
  _vertices = cont._vertices;
  _mcvertices = cont._mcvertices;
  _boxvol = cont._boxvol;

  return *this;
}

ProbabilityContour::~ProbabilityContour() {}

TH1F* ProbabilityContour::Contour1D(size_t idx) const {

  // Return if the dimension does not exist on the grid
  TH1F* contour = NULL;
  if ( idx > _grid.GetDimension() ) {
    throw(Exceptions::Exception(Exceptions::recoverable,
	  "The axis requested is not provided in this space",
	  "ProbabilityContour::Contour1D"));
    return contour;
  }

  // Initialize the histogram
  contour = new TH1F(TString::Format("contour_axis_id%d", (int)idx),
	    	     TString::Format("Optimal 2D contour on axes id %d", (int)idx),
	     	     _grid.GetNpoints(idx), _grid.GetLowerBound(idx), _grid.GetUpperBound(idx));
  contour->SetLineWidth(2);
  contour->SetLineColor(kGreen+2);
  contour->SetMarkerColor(kGreen+2);

  // Fill the histogram
  if ( _algo == "rectangle" || _algo == "trapezoid" ) {
    size_t i;
    for (i = 0; i < _indices.size(); i++) {
      contour->SetBinContent(_indices[i][idx]+1, _level);
      contour->SetBinError(_indices[i][idx]+1, 1e-9);
    }

  } else if ( _algo == "mc" || _algo == "hull" || _algo == "alphac" ) {
    size_t i;
    for (i = 0; i < _grid.GetNpoints(idx); i++)
      if ( _grid[i].GetValue() > _level ) {
	contour->SetBinContent(i+1, _level);
	contour->SetBinError(i+1, 1e-9);
      }

  } else {
    throw(Exceptions::Exception(Exceptions::recoverable,
	  "Algorithm not supported: "+_algo,
	  "ProbabilityContour::Contour1D"));
  }

  return contour;
}

TH2F* ProbabilityContour::Contour2D(size_t idx, size_t idy) const {

  // Return if the dimension does not exist or the list does not contain at least two points
  TH2F* contour = NULL;
  if ( idx > _grid.GetDimension() || idy > _grid.GetDimension() ) {
    throw(Exceptions::Exception(Exceptions::recoverable,
	  "The axes requested are not all provided in this space",
	  "ProbabilityContour::Contour2D"));
  }

  // Initialize the histogram
  contour = new TH2F(TString::Format("contour_axis_id%d-%d", (int)idx, (int)idy),
	    	     TString::Format("Optimal 2D contour on axes id %d-%d", (int)idx, (int)idy),
	     	     _grid.GetNpoints(idx), _grid.GetLowerBound(idx), _grid.GetUpperBound(idx),
	     	     _grid.GetNpoints(idy), _grid.GetLowerBound(idy), _grid.GetUpperBound(idy));
  contour->SetContour(1);
  contour->SetContourLevel(0, 1.);
  contour->SetLineWidth(2);
  contour->SetLineColor(kGreen+2);
  contour->SetFillStyle(1000);
  contour->SetFillColorAlpha(kGreen+2, 0.25);

  // Fill the histogram
  if ( _algo == "rectangle" || _algo == "trapezoid" ) {
    size_t i;
    for (i = 0; i < _indices.size(); i++)
        contour->SetBinContent(_indices[i][idx]+1, _indices[i][idy]+1, 1.);

  } else if ( _algo == "mc" || _algo == "hull" || _algo == "alphac" ) {
    size_t i, j;
    for (i = 0; i < _grid.GetNpoints(idx); i++)
      for (j = 0; j < _grid.GetNpoints(idy); j++)
        if ( _grid[std::vector<size_t>({i, j})].GetValue() > _level )
	    contour->SetBinContent(i+1, j+1, 1);
  } else {
    throw(Exceptions::Exception(Exceptions::recoverable,
	  "Algorithm not supported: "+_algo,
	  "ProbabilityContour::Contour2D"));
  }

  return contour;
}

TH3F* ProbabilityContour::Contour3D(size_t idx, size_t idy, size_t idz) const {

  // Return if the dimension does not exist or the list does not contain at least two points
  TH3F* contour = NULL;
  if ( idx > _grid.GetDimension() || idy > _grid.GetDimension() || idz > _grid.GetDimension() ) {
    throw(Exceptions::Exception(Exceptions::recoverable,
	  "The axes requested are not all provided in this space",
	  "ProbabilityContour::Contour3D"));
    return contour;
  }

  // Initialize the histogram
  contour = new TH3F(TString::Format("contour_axis_id%d-%d-%d", (int)idx, (int)idy, (int)idz),
	     	     TString::Format("Optimal 3D contour on axes id %d-%d-%d",
			     	     (int)idx, (int)idy, (int)idz),
	     	     _grid.GetNpoints(idx), _grid.GetLowerBound(idx), _grid.GetUpperBound(idx),
	     	     _grid.GetNpoints(idy), _grid.GetLowerBound(idy), _grid.GetUpperBound(idy),
	     	     _grid.GetNpoints(idz), _grid.GetLowerBound(idz), _grid.GetUpperBound(idz));
  contour->SetContour(1);
  contour->SetContourLevel(0, 1.);
  contour->SetLineWidth(2);
  contour->SetLineColor(kGreen+2);
  contour->SetFillColorAlpha(kGreen+2, 0.25);
  
  // Fill the histogram
  if ( _algo == "rectangle" || _algo == "trapezoid" ) {
    size_t i;
    for (i = 0; i < _indices.size(); i++)
        contour->SetBinContent(_indices[i][idx]+1, _indices[i][idy]+1, _indices[i][idz]+1, 1.);

  } else if ( _algo == "mc" || _algo == "hull" || _algo == "alphac" ) {
    size_t i, j, k;
    for (i = 0; i < _grid.GetNpoints(idx); i++)
      for (j = 0; j < _grid.GetNpoints(idy); j++)
        for (k = 0; k < _grid.GetNpoints(idy); k++)
        if ( _grid[std::vector<size_t>({i, j, k})].GetValue() > _level )
	    contour->SetBinContent(i+1, j+1, k+1, 1);
  } else {
    throw(Exceptions::Exception(Exceptions::recoverable,
	  "Algorithm not supported: "+_algo,
	  "ProbabilityContour::Contour3D"));
  }

  return contour;
}

void ProbabilityContour::Evaluate(double alpha, double eps) {

  _alpha = alpha;

  try {
    if ( _algo == "alphac" ) {
      AlphaComplexMethod(alpha, eps);
    } else if ( _algo == "hull" ) {
      HullMethod(alpha, eps);
    } else if ( _algo == "mc" ) {
      MCMethod(alpha, eps);
    } else if ( _algo == "rectangle" ) {
      RectangleMethod(alpha, eps);
    } else if ( _algo == "trapezoid" ) {
      TrapezoidMethod(alpha, eps);
    } else {
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	    "Algorithm not supported: "+_algo,
	    "ProbabilityContour::Evaluate"));
    }
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "The estimation using the "+_algo+" method failed"+std::string(e.what()),
	  "ProbabilityContour::Evaluate"));
  }
}

double ProbabilityContour::GetContourVolumeError() const {

  if ( _algo == "mc" ) {
    // The error on the estimated volume follows a binomial distribution
    double occ = _volume/_boxvol;
    return _volume*sqrt((1.-occ)/(occ*_mcvertices.size()));
  }

  throw(Exceptions::Exception(Exceptions::recoverable,
	"Algorithm not supported: "+_algo,
	"ProbabilityContour::GetContourVolumeError"));
  return 0.;
}

void ProbabilityContour::AlphaComplexMethod(double alpha, double eps) {

  // For a given fraction alpha, simply find the level of the (1-alpha)*N^th point (vertices sorted)
  size_t id = (1.-alpha)*(_vertices.size()-1);
  _level = _vertices[id].GetValue();

  // For a given level, one can estimate the expected average distance between points (1NN)
  double n = _vertices[0].size();
  double R = .98*pow(_level*_vertices.size(), -1./n);

  // Given the mean radius, use alpha=1/R as a parameter to build the alpha-complex
  std::vector<std::vector<double>> points;
  size_t i;
  for (i = id; i < _vertices.size(); i++)
      points.push_back(_vertices[i].GetCoordinates());

  AlphaComplex alphac(points, 1./R);
  _volume = alphac.GetVolume();
}

void ProbabilityContour::HullMethod(double alpha, double eps) {

  // For a given fraction alpha, simply find the level of the (1-alpha)*N^th point (vertices sorted)
  size_t id = (1.-alpha)*(_vertices.size()-1);
  _level = _vertices[id].GetValue();

  // Compute the convexe hull of the core alpha fraction of the points
  std::vector<double> points;
  size_t i, j;
  for (i = id; i < _vertices.size(); i++)
    for (j = 0; j < _vertices[0].size(); j++)
  	points.push_back(_vertices[i][j]);

  orgQhull::Qhull qhull("", _vertices[0].size(),
	points.size()/_vertices[0].size(), &(points[0]), "");
  _volume = qhull.volume();
}

void ProbabilityContour::MCMethod(double alpha, double eps) {

  // For a given fraction alpha, simply find the level of the (1-alpha)*N^th point (vertices sorted)
  size_t id = (1.-alpha)*(_vertices.size()-1);
  _level = _vertices[id].GetValue();

  // Compute the contour by finding the fraction of MC points higher than this level
  size_t i;
  double n(0.);
  for (i = 0; i < _mcvertices.size(); i++)
    if ( _mcvertices[i].GetValue() > _level )
	n++;

  _volume = (n/_mcvertices.size())*_boxvol;
}

void ProbabilityContour::RectangleMethod(double alpha, double eps) {

  // On a regular grid, a cell volume is always the same
  double cellvol = _grid.GetCellVolume();

  // Start from the top, increment volume until the precision decreases or
  // the requested integrated probability is exceeded (will not get closer)
  _volume = 0;
  _indices.resize(0);
  double prob(0);
  size_t i = _vertices.size()-1;
  while ( prob < alpha ) {
    // Increment the probability if it makes it better
    if ( fabs(prob+cellvol*_vertices[i].GetValue()-alpha) > fabs(prob-alpha) )
	break;
    prob += cellvol*_vertices[i].GetValue();

    // Increment the volume
    _volume += cellvol;

    // Store which of the elements were used
    _indices.push_back(_vertices[i].GetIndices());

    // Increment the index
    i--;
  }

  // The level of the contour is the lowest probabability attained
  _level = _grid[_indices.back()].GetValue();

  // Flag if the desired precision was not achieved
  if ( fabs(prob-alpha) > eps )
      Pitch::print(Pitch::warning, "ProbabilityContour::RectangleMethod", 
	"The required precision was not achieved. Spread: "+std::to_string(fabs(prob-alpha)));
}

void ProbabilityContour::TrapezoidMethod(double alpha, double eps) {

  // Throw as the method is currently not supported (TODO TODO TODO)
  throw(Exceptions::Exception(Exceptions::nonRecoverable,
	"Algorithm currently not supported",
	"ProbabilityContour::TrapezoidMethod"));

  // On a regular grid, a cell volume is always the same
  double cellvol = _grid.GetCellVolume();

  // In the traperoidal case, loop over the cells rather than vertices
  // Compute the central value of the function as the centre of each cell, then sort
  std::vector<Grid> cellarray = _grid.GetCellArray();
  std::vector<std::pair<std::vector<size_t>, double>> central;
  size_t i;
  double val;
  for (const Grid& cell : _grid.GetCellArray()) {
    val = 0.;
    for(const Vertex& vertex : cell.GetVertexArray())
	val += vertex.GetValue();
    val /= pow(2, _grid.GetDimension());
    central.push_back(std::pair<std::vector<size_t>, double>(cell[0].GetIndices(), val));
  }
  std::sort(central.begin(), central.end(),
	[] (const std::pair<std::vector<size_t>, double>& a,
	    const std::pair<std::vector<size_t>, double>& b) {
		return a.second < b.second;
	});

  // Start from the top, increment volume until the precision decreases or
  // the requested integrated probability is exceeded (will not get closer)
  _volume = 0;
  _indices.resize(0);
  double prob(0);
  i = central.size()-1;
  while ( prob < alpha ) {
    // Increment the probability if it makes it better
    if ( fabs(prob+cellvol*central[i].second-alpha) > fabs(prob-alpha) )
	break;
    prob += cellvol*central[i].second;

    // Increment the volume
    _volume += cellvol;

    // Store which of the elements were used
    _indices.push_back(central[i].first);

    // Increment the index
    i--;
  }
}
