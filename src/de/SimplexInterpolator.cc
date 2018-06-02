#include "SimplexInterpolator.hh"

SimplexInterpolator::SimplexInterpolator() :
  _params(0) {
}

SimplexInterpolator::SimplexInterpolator(const std::vector<Vertex>& vertices) :
  _params(0) {

  try {
    Initialize(vertices);
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::recoverable,
	  "Could not initialize this array of vertices"+std::string(e.what()),
	  "SimplexInterpolator::SimplexInterpolator"));
  }
}

SimplexInterpolator::SimplexInterpolator(const SimplexInterpolator& sint) {
  *this = sint;
}

SimplexInterpolator& SimplexInterpolator::operator=(const SimplexInterpolator& sint) {
  if ( this == &sint )
      return *this;

  _params = sint._params;

  return *this;
}

SimplexInterpolator::~SimplexInterpolator () {}

void SimplexInterpolator::Initialize(const std::vector<Vertex>& vertices) {

  // Check that the grid has the right dimensions
  size_t N = vertices.size();
  size_t i, j;
  Assert::IsNotEmpty("SimplexInterpolator::Initialize", "Vector of vertices", vertices);
  for (i = 0; i < N; i++)
    if ( vertices[i].size() != N-1 )
        throw(Exceptions::Exception(Exceptions::nonRecoverable,
	      "The "+std::to_string(i)+"th vertex of the array is not"
	      "of dimension N-1. Not an simplex, cannot initialize",
	      "SimplexInterpolator::Initialize"));

  // Compute the base matrix, need to find the n-hyperplane in the (n+1) space
  Matrix<double> base(N, N);

  // Loop over the vertices (N of them), fill the first column with 1, others with coordinates
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
	base[i][j] = !j ? 1. : vertices[i][j-1];

  // If the base matrix is singular, the points are collinear and do not form an n-simplex
  if ( !base.Determinant() )
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	    "The simplex vertices are affinely dependant, cannot intialize",
	    "SimplexInterpolator::Initialize"));

  // Solve in the parameters
  std::vector<double> v(N);  // Vector that contains the values of the function
  for (i = 0; i < N; i++)
      v[i] = vertices[i].GetValue();

  _params = base.Solve(v);
}

double SimplexInterpolator::Evaluate(const std::vector<double>& v) const {

  // Need to be of the same dimension as the rest
  Assert::IsEqual("SimplexInterpolator::Evaluate", 
	"Dimension and vector size", _params.size()-1, v.size());

  // Loop over all the dimensions and increment the vector product _params*v
  double function(0.);
  size_t i;
  for (i = 0; i <_params.size(); i++)
      function += !i ? _params[i] : v[i-1]*_params[i];

  // Return the interpolated function in the n-point
  return function;
}
