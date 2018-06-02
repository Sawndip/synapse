#include "Histogram.hh"

Histogram::Histogram() :
  _grid(), _points(0), _contents(0), _outflows(0) {
}

Histogram::Histogram(std::vector<size_t> number,
	    	     const std::vector<double>& lower,
	    	     const std::vector<double>& upper) :
  _grid(), _points(0), _contents(0), _outflows(0) {

  try {
    // The number of points on the grid is N+1 the number of bins
    size_t i;
    for (i = 0; i < number.size(); i++)
        number[i]++;
    _grid = Grid(number, lower, upper);

    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "Histogram::Histogram"));
  }
}

Histogram::Histogram(const Histogram& hist) {
  *this = hist;
}

Histogram& Histogram::operator=(const Histogram& hist) {
  if ( this == &hist )
      return *this;

  _grid = hist._grid;
  _points = hist._points;
  _contents = hist._contents;
  _outflows = hist._outflows;

  return *this;
}

Histogram::~Histogram () {}

void Histogram::Initialize() {

  // Initialize the contents array and the over/underflow array
  size_t ntotal = _grid.GetVertexArray().size();
  _contents.resize(ntotal);
  std::fill(_contents.begin(), _contents.end(), 0.);

  _outflows.resize(pow(3, _grid.GetDimension()));
  std::fill(_outflows.begin(), _outflows.end(), 0.);
}

void Histogram::Fill(const std::vector<double>& point, const double w) {

  // Check that the data point is of the same dimension as the histogram
  size_t dim = _grid.GetDimension();
  Assert::IsEqual("Histogram::Fill", "Dimension and point size", dim, point.size());

  // Append the array
  _points.push_back(point);

  // If the point is outside, find the corresponding under/overflow bin (3^_dim-1 bins)
  size_t i;
  std::vector<size_t> trits;
  bool inside(true);
  for (i = 0; i < dim; i++) {
    if ( point[i] < _grid.GetLowerBound(i) ) {
      trits.push_back(0);
      inside = false;
    } else if ( point[i] > _grid.GetUpperBound(i) ) {
      trits.push_back(2);
      inside = false;
    } else {
      trits.push_back(1);
    }
  }

  // Look for the cell in which the point is contained, fill it
  size_t id;
  std::vector<size_t> indices;
  if ( !inside ) {
    // Get the outflow bin ID, the first axis is the most significant
    id = 0;
    for (i = 0; i < dim; i++)
	id += trits[i]*pow(3, dim-1-i);
    _outflows[id] += w;
  } else {
    indices = _grid.GetLowerIndices(point);
    _contents[_grid.ID(indices)] += w;
  }
}

void Histogram::FillN(const std::vector<std::vector<double>>& points,
		      const std::vector<double> w) {

  size_t i;
  if ( !w.size() ) {
    for (i = 0; i < points.size(); i++)
        this->Fill(points[i]);

    return;
  }

  for (i = 0; i < points.size(); i++)
      this->Fill(points[i], w[i]);
}

Bin Histogram::GetBin(const std::vector<size_t>& indices) const {

  Grid cell = _grid.GetCell(indices);
  return Bin(indices,
	     cell.GetLowerBoundArray(),
	     cell.GetUpperBoundArray(),
	     _contents[_grid.ID(indices)]);
}

std::vector<Bin> Histogram::GetBinArray() const {

  std::vector<size_t> indices;
  std::vector<Bin> bins;
  bool edge;
  size_t i, j;
  for (i = 0; i < _contents.size(); i++) {

    indices = _grid.Indices(i);
    edge = false;
    for (j = 0; j < _grid.GetDimension(); j++)
      if ( indices[j] == _grid.GetNpoints(j)-1 )
	edge = true;
    if ( !edge )
	bins.push_back(GetBin(_grid.Indices(i)));
  }

  return bins;
}

std::vector<double> Histogram::GetBinBarycentre(const std::vector<size_t>& indices) const {

  std::vector<double> barycentre;
  size_t i;
  for (i = 0; i < _grid.GetDimension(); i++)
      barycentre.push_back(GetBinCentre(indices[i], i));

  return barycentre;
}

double Histogram::GetBinCentre(const size_t index, const size_t i) const {

  // Check that the index is not out of bound for axis i
  size_t nbins = _grid.GetNpoints(i)-1;
  if ( index >= nbins  ) {
    throw(Exceptions::Exception(Exceptions::recoverable,
	  "Index out of bounds",
	  "Histogram::GetBinCentre"));
    return 0.;
  }

  // Return the centre of the bin
  double min = _grid.GetLowerBound(i);
  return min + (i+.5)*_grid.GetPeriod(i);
}

void Histogram::SetBinContent(const std::vector<size_t>& indices, const double value) {

  _contents[_grid.ID(indices)] = value;
}

const double& Histogram::GetBinContent(const std::vector<size_t>& indices) const {

  return _contents[_grid.ID(indices)];
}

double Histogram::GetBinDensity(const std::vector<size_t>& indices) const {

  return GetBinContent(indices)/(_points.size()*GetBinVolume());
}

std::vector<size_t> Histogram::GetBinID(const std::vector<double>& point) const {

  return _grid.GetLowerIndices(point);
}

double Histogram::GetBinVolume() const {

  double volume(1.);
  size_t i;
  for (i = 0; i < _grid.GetDimension(); i++)
      volume *= GetBinWidth(i);

  return volume;
}

Grid Histogram::GetInterpolationGrid() const {

  // An interpolation grid has one vertex per bin. The vertex is located
  // at its centre and has the local density as its mapping.
  std::vector<Vertex> points;
  for (const Bin& bin : this->GetBinArray() )
      points.push_back(Vertex(bin.GetBarycentre(), bin.GetContent()));
  
  return Grid(points);
}

size_t Histogram::GetNbins() const {

  size_t i;
  size_t N(1);
  for (i = 0; i < _grid.GetDimension(); i++)
      N *= (_grid.GetNpoints(i)-1);

  return N;
}

void Histogram::Scale(const double factor) {

  size_t i;
  for (i = 0; i < _contents.size(); i++)
      _contents[i] *= factor;
}

void Histogram::ScaleToDensity(const double norm) {

  size_t i;
  for (i = 0; i < _contents.size(); i++)
      _contents[i] = norm*GetBinDensity(_grid.Indices(i));
}

TObject* Histogram::ToROOT(const std::string name) const {

  size_t dim = _grid.GetDimension();
  if ( dim == 1 ) {
    TH1F* hist = new TH1F(name.c_str(), name.c_str(),
	_grid.GetNpoints(0)-1, _grid.GetLowerBound(0), _grid.GetUpperBound(0));

    size_t i;
    for (i = 0; i < _grid.GetNpoints(0)-1; i++)
      hist->SetBinContent(i+1, _contents[_grid.ID({i})]);

    return hist;
  } else if ( dim == 2 ) {
    TH2F* hist = new TH2F(name.c_str(), name.c_str(), 
	_grid.GetNpoints(0)-1, _grid.GetLowerBound(0), _grid.GetUpperBound(0),
	_grid.GetNpoints(1)-1, _grid.GetLowerBound(1), _grid.GetUpperBound(1));

    size_t i, j;
    for (i = 0; i < _grid.GetNpoints(0)-1; i++)
      for (j = 0; j < _grid.GetNpoints(1)-1; j++)
	hist->SetBinContent(i+1, j+1, _contents[_grid.ID({i, j})]);

    return hist;
  } else if ( dim == 3 ) {
    TH3F* hist = new TH3F(name.c_str(), name.c_str(), 
	_grid.GetNpoints(0)-1, _grid.GetLowerBound(0), _grid.GetUpperBound(0),
	_grid.GetNpoints(1)-1, _grid.GetLowerBound(1), _grid.GetUpperBound(1),
	_grid.GetNpoints(2)-1, _grid.GetLowerBound(2), _grid.GetUpperBound(2));

    size_t i, j, k;
    for (i = 0; i < _grid.GetNpoints(0)-1; i++)
      for (j = 0; j < _grid.GetNpoints(1)-1; j++)
        for (k = 0; k < _grid.GetNpoints(2)-1; k++)
	  hist->SetBinContent(i+1, j+1, k+1, _contents[_grid.ID({i, j, k})]);

    return hist;
  }

  throw(Exceptions::Exception(Exceptions::recoverable,
	"ROOT does not support histograms of dimension "+std::to_string(dim),
	"Histogram::ToROOT"));
  return (TObject*)NULL;
}
