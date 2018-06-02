#include "LinearInterpolator.hh"

LinearInterpolator::LinearInterpolator() :
  _dim(0), _params(0) {
}

LinearInterpolator::LinearInterpolator(const Grid& grid) :
  _dim(grid.GetDimension()), _params(0) {

  try {
    Initialize(grid);
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::recoverable,
	  "Could not initialize this grid of vertices"+std::string(e.what()),
	  "LinearInterpolator::LinearInterpolator"));
  }
}

LinearInterpolator::LinearInterpolator(const LinearInterpolator& lint) {
  *this = lint;
}

LinearInterpolator& LinearInterpolator::operator=(const LinearInterpolator& lint) {
  if ( this == &lint )
      return *this;

  _dim = lint._dim;
  _params = lint._params;

  return *this;
}

LinearInterpolator::~LinearInterpolator () {}

void LinearInterpolator::Initialize(const Grid& grid) {

  // Check that the grid has the right dimensions
  Assert::IsNonZero("LinearInterpolator::Initialize", "Dimension", _dim);
  size_t i, j;
  for (i = 0; i < _dim; i++)
    if ( grid.GetNpoints(i) != 2 )
        throw(Exceptions::Exception(Exceptions::nonRecoverable,
	      "The "+std::to_string(i)+"th axis of the grid does not cointain"
	      "exactly two points. Not an n-orthotope, cannot initialize.",
	      "LinearInterpolator::Initialize"));

  // Compute the base matrix
  size_t N = grid.size();
  Matrix<double> base(N, N);

  // Loop over the vertices (N of them)
  for (i = 0; i < N; i++) {
    // Loop over all the possible combinations
    for (j = 0; j < N; j++) {

      // Set the ij matrix element, ith vertex and jth combination of its coordinates
      // (j & ( 1 << k )) >> k corresponds to the k^th bit of j
      double element = 1.;
      for (size_t k = 0; k < _dim; k++)
	  element *= pow(grid[i][k], (j & ( 1 << k )) >> k);

      base[i][j] = element;
    }
  }

  // Solve in the parameters by inverting the base matrix
  std::vector<double> v(N);  // Column vector that contains the values of the function
  for (i = 0; i < N; i++)
      v[i] = grid[i].GetValue();

  _params = base.Solve(v);
}

double LinearInterpolator::Evaluate(const std::vector<double>& v) const {

  // Need to be of the same dimension as the rest
  Assert::IsEqual("LinearInterpolator::Evaluate", "Dimension and vector size", _dim, v.size());

  // Loop over all the possible combinations, increment function
  size_t N = size_t(pow(2, _dim));
  double function(0.);
  for (size_t i = 0; i < N; i++) {

    // Set the ith combination of the coordinates
    // (i & ( 1 << j )) >> j corresponds to the j^th bit of i
    double element = 1.;
    for (size_t j = 0; j < _dim; j++)
	element *= pow(v[j], (i & ( 1 << j )) >> j);

    function += element*_params[i];
  }

  // Return the interpolated function in the n-point
  return function;
}
