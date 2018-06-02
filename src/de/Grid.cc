#include "Grid.hh"

Grid::Grid() :
  _dim(0), _vertices(0), _number(0), _step(0), _period(0), _lower(0), _upper(0) {
}

Grid::Grid(const std::vector<Vertex>& vertices) :
  _dim(0), _vertices(vertices), _number(0), _step(0), _period(0), _lower(0), _upper(0) {

  try {
    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize from vertices"+std::string(e.what()),
	  "Grid::Grid"));
  }
}

Grid::Grid(const Matrix<double>& mat, const std::vector<double>& vec) :
  _dim(0), _vertices(), _number(0), _step(0), _period(0), _lower(0), _upper(0) {

  try {
    InitializeFromMatrix(mat, vec);
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize from matrix"+std::string(e.what()),
	  "Grid::Grid"));
  }
}

Grid::Grid(double (*fcn)(const std::vector<double>& val),
       	   const std::vector<size_t>& number,
       	   const std::vector<double>& lower,
       	   const std::vector<double>& upper) :
  _dim(0), _vertices(0), _number(number), _step(0), _period(0), _lower(lower), _upper(upper) {

  try {
    InitializeFromFunction(fcn);
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize from function"+std::string(e.what()),
	  "Grid::Grid"));
  }
}

Grid::Grid(const std::vector<size_t>& number,
       	   const std::vector<double>& lower,
       	   const std::vector<double>& upper,
	   const std::string algo) :
  _dim(0), _vertices(0), _number(number), _step(0), _period(0), _lower(lower), _upper(upper) {

  try {
    InitializeRandom(algo);
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize from random"+std::string(e.what()),
	  "Grid::Grid"));
  }
}

Grid::Grid(const Grid& grid) {
  *this = grid;
}

Grid& Grid::operator=(const Grid& grid) {
  if ( this == &grid )
      return *this;

  _dim = grid._dim;
  _vertices = grid._vertices;
  _number = grid._number;
  _step = grid._step;
  _period = grid._period;
  _lower = grid._lower;
  _upper = grid._upper;

  return *this;
}

Grid::~Grid () {}

void Grid::Initialize() {

  // First check that we have enough vertices to build a grid
  Assert::IsNotEmpty("Grid::Initialize", "Vector of input vertices", _vertices);
  Assert::IsNonZero("Grid::Initialize", "Dimension", _vertices[0].size());
  size_t N = _vertices.size();
  _dim = _vertices[0].size();
  if ( N < pow(2, _dim) )
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	    "Not enough vertices to build a grid, need at least 2^d",
	    "Grid::Initialize"));

  // Figure out the lower and upper boundaries if not provided
  size_t i;
  if ( !_lower.size() && !_upper.size() ) {
    for (i = 0; i < _dim; i++) {
      _lower.push_back(_vertices[0][i]);
      _upper.push_back(_vertices[N-1][i]);
    }
  }

  // Check that the boundaries are ordered properly and define a non-zero interval
  for (i = 0; i < _dim; i++)
      Assert::IsGreater("Grid::Initialize",
	"The interval on axis "+std::to_string(i), _upper[i]-_lower[i], 0., true);

  // Find the number of points, stepping and periodicity in each of the projection
  _number.resize(_dim);
  _step.resize(_dim);
  _period.resize(_dim);
  int j;
  for (j = _dim-1; j >= 0; j--) {

    // Get the stepping of the axis
    _step[j] = 1.;
    for (i = j+1; i < _dim; i++)
	_step[j] *= _number[i];

    // Get the number of points on this axis according to the stepping
    for (i = 0; i < N+1; i += _step[j]) {
      if ( i == N ) {
	_number[j] = N/_step[j];
	break;
      } else if ( i && _vertices[i][j] == _lower[j] ) {
	_number[j] = i/_step[j];
	break;
      }
    }

    // Get the period between two steps
    _period[j] = (_upper[j]-_lower[j])/(_number[j]-1);
  }
}

void Grid::InitializeFromMatrix(const Matrix<double>& mat, const std::vector<double>& vec) {

  // Check that the vertex matrix has the right dimensions
  Assert::IsNonZero("Grid::InitializeFromMatrix", "Matrix dimension", mat.Nrows());  
  _dim = mat.Nrows();
  if ( mat.Ncols() < pow(2, _dim) )
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	    "Not enough vertices to build a grid, need at least 2^d",
	    "Grid::InitializeFromMatrix"));

  // Check that we have a value for each of the N=N_1*...*N_n vertices
  Assert::SameSize("Grid::InitializeFromMatrix", "Grid point vector and value vector", mat[0], vec);

  // Create a vertex vector corresponding to the input
  size_t i;
  for (i = 0; i < mat.Ncols(); i++)
      _vertices.push_back(Vertex(mat.Column(i), vec[i]));

  // Initialize the grid
  Initialize();
}

void Grid::InitializeFromFunction(double (*fcn)(const std::vector<double>& val)) {

  // Check that we have at least two vertices per dimension
  Assert::IsNonZero("Grid::InitializeFromFunction", "Dimension", _number.size());
  _dim = _number.size();
  size_t i, j;
  for (i = 0; i < _dim; i++)
    if ( _number[i] < 2 )
        throw(Exceptions::Exception(Exceptions::nonRecoverable,
	      "The number of points must be two or more per dimension",
	      "Grid::InitializeFromFunction"));

  // Check that the boundaries are ordered properly and define a non-zero interval
  for (i = 0; i < _dim; i++)
      Assert::IsGreater("Grid::InitializeFromFunction",
	"The interval on axis "+std::to_string(i), _upper[i]-_lower[i], 0., true);

  // The total number of vertices is the product of every number of points N_1*...*N_n
  size_t N = 1;
  for (i = 0; i < _dim; i++)
      N *= _number[i];
  _vertices.resize(N);

  // Get the stepping and periodicity of each of the grid axes
  _step.resize(_dim);
  _period.resize(_dim);
  for (i = 0; i < _dim; i++) {

    // Get the stepping of the axis
    _step[i] = 1;
    for (j = i+1; j < _dim; j++)
	_step[i] *= _number[j];

    // Get the period between two steps
    _period[i] = (_upper[i]-_lower[i])/(_number[i]-1);
  }

  // Fill the grid of vertices, read the indices using the stepping
  double value;
  std::vector<size_t> indices(_dim);
  std::vector<double> coords(_dim);
  for (j = 0; j < N; j++) {

    // Read the indices for each dimension separately and get the corresponding position
    indices = Indices(j);
    for (i = 0; i < _dim; i++)
        coords[i] = _lower[i]+indices[i]*_period[i];

    // Evaluate the function in the coords
    value = fcn(coords);

    // Set the vertex
    _vertices[j] = Vertex(coords, value, indices);
  }
}

void Grid::InitializeRandom(const std::string& algo) {

  // Check that we have at least two vertices per dimension
  Assert::IsNonZero("Grid::InitializeRandom", "Dimension", _number.size());
  _dim = _number.size();
  size_t i, j;
  for (i = 0; i < _dim; i++)
    if ( _number[i] < 2 )
        throw(Exceptions::Exception(Exceptions::nonRecoverable,
	      "The number of points must be two or more per dimension",
	      "Grid::InitializeRandom"));

  // Check that the boundaries are ordered properly and define a non-zero interval
  for (i = 0; i < _dim; i++)
      Assert::IsGreater("Grid::InitializeRandom",
	"The interval on axis "+std::to_string(i), _upper[i]-_lower[i], 0., true);

  // The total number of vertices is the product of every number of points N_1*...*N_n
  size_t N = 1;
  for (i = 0; i < _dim; i++)
      N *= _number[i];
  _vertices.resize(N);

  // Get the stepping and periodicity of each of the grid axes
  _step.resize(_dim);
  _period.resize(_dim);
  for (i = 0; i < _dim; i++) {

    // Get the stepping of the axis
    _step[i] = 1;
    for (j = i+1; j < _dim; j++)
	_step[i] *= _number[j];

    // Get the period between two steps
    _period[i] = (_upper[i]-_lower[i])/(_number[i]-1);
  }

  // Fill the grid of vertices, read the indices using the stepping
  std::vector<size_t> indices(_dim);
  std::vector<double> coords(_dim);
  for (j = 0; j < N; j++) {

    // Read the indices for each dimension separately and get the corresponding position
    indices = Indices(j);
    for (i = 0; i < _dim; i++)
        coords[i] = _lower[i]+indices[i]*_period[i];

    // Set the vertex
    if ( algo == "flat" ) {
      _vertices[j] = Vertex(coords, 1., indices);
    } else if ( algo == "uniform" ) {
      _vertices[j] = Vertex(coords, (double)rand()/RAND_MAX, indices);
    } else {
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	    "Grid generation mapping algorithm not recognized: "+algo,
	    "Grid::InitializeRandom"));
    }
  }
}

Grid Grid::GetBestCell(const std::vector<double>& v) const {

  std::vector<size_t> minids = this->GetLowerIndices(v);
  return this->GetCell(minids);
}

std::vector<size_t> Grid::GetLowerIndices(const std::vector<double>& v) const {

  // For each axis, find the ID of the point that immidiately preceeds v
  std::vector<size_t> minids(_dim);
  size_t i ;
  for (i = 0; i < _dim; i++) {

    // If v[i] is beyond the grid, extrapolate from closest instead
    if ( v[i] <= _lower[i] ) {
      minids[i] = 0;
      continue;
    } else if ( v[i] >= _upper[i]-_period[i] ) {
      minids[i] = _number[i]-2;
      continue;
    }

    minids[i] = (size_t)std::floor((v[i]-_lower[i])/_period[i]);
  }

  return minids;
}

Grid Grid::GetCell(const std::vector<size_t>& indices) const {

  // Set a grid of 2^n vertices that constitute the cell
  size_t N = pow(2, _dim);		// Number of vertices in an n-orthotope
  std::vector<Vertex> vertices;
  std::vector<size_t> ids(_dim);
  size_t i, j, k;
  for (i = 0; i < N; i++) {

    // Choose the right vertex on the grid
    // Read the jth bit of i as the indication for incrementation
    for (j = 0; j < _dim; j++) {
      k = (i & ( 1 << j )) >> j;
      ids[_dim-1-j] = indices[_dim-1-j] + k;
    }

    vertices.push_back(_vertices[ID(ids)]);
  }

  return Grid(vertices);
}

std::vector<Grid> Grid::GetCellArray() const {

  // Loop over all the possible indices, skip the upper edge, increment the vector
  std::vector<Grid> cells;
  size_t i, j;
  std::vector<size_t> indices;
  bool edge;
  for (i = 0; i < _vertices.size(); i++) {
    // Check that it is not on the upper edge
    edge = false;
    indices = Indices(i);
    for (j = 0; j <_dim; j++)
      if ( indices[j] == _number[j]-1 )
	edge = true;

    // If not on the uppe edge, simply increment the vector
    if ( !edge )
        cells.push_back(GetCell(indices));
  }

  return cells;
}

double Grid::GetCellVolume() const {

  double cellvol(1.);
  size_t i;
  for (i = 0; i < _dim; i++)
      cellvol *= _period[i];
  return cellvol;
}

double Grid::GetPosition(const size_t i, const size_t id) const {
  return _lower[i]+id*_period[i];
}

std::vector<double> Grid::GetPosition(const std::vector<size_t>& vi,
				      const std::vector<size_t>& ids) const {

  std::vector<double> pos(vi.size());
  size_t i;
  for (i = 0; i < vi.size(); i++)
	pos[i] = _lower[vi[i]]+ids[i]*_period[vi[i]];

  return pos;
}

double Grid::GetVolume() const {

  double vol(1.);
  size_t i;
  for (i = 0; i < _dim; i++)
      vol *= (_upper[i]-_lower[i]);
  return vol;
}

size_t Grid::ID(const std::vector<size_t>& indices) const {

  size_t i, id(0);
  for (i = 0; i < indices.size(); i++)
      id += indices[i]*_step[i];

  return id;
}

std::vector<size_t> Grid::Indices(size_t id) const {

  // Read the indices for each dimension separately
  size_t i;
  std::vector<size_t> indices(_dim);
  for (i = 0; i < _dim; i++) {
    indices[i] = id/_step[i];
    id -= indices[i]*_step[i];
  }

  return indices;
}

bool Grid::IsInside(const std::vector<double>& v) const {

  size_t i;
  for (i = 0; i < _dim; i++)
    if ( v[i] < _lower[i] || v[i] > _upper[i] )
	return false;

  return true;
}

std::vector<TPolyLine*> Grid::Polygons() const {

  // Throw if the wrong dimension is requested
  std::vector<TPolyLine*> polylines;
  Assert::IsEqual("Grid::Polygons", "Dimension", _dim, (size_t)2);

  // Loop over the cells, define a polygon per cell
  size_t i, j;
  std::vector<std::vector<size_t>> aclockw = {{0, 0}, {1, 0}, {1, 1}, {0, 1}, {0, 0}};
  std::vector<double> x, y;
  for (i = 0; i < _number[0]-1; i++) {
    for (j = 0; j < _number[1]-1; j++) {
      x.resize(0);
      y.resize(0);
      for (const std::vector<size_t> id : aclockw ) {
	x.push_back(_lower[0]+(i+id[0])*_period[0]);
	y.push_back(_lower[1]+(j+id[1])*_period[1]);
      }

      polylines.push_back(new TPolyLine(5, &(x[0]), &(y[0])));
    }
  }

  return polylines;
}
