#include "Interpolator.hh"

Interpolator::Interpolator() :
  _dim(0), _algo(""), _grid(), _del(), _ex(false) {
}

Interpolator::Interpolator(const Grid& grid,
			   const bool ex,
			   const std::string algo) :
  _dim(grid.GetDimension()), _algo(algo), _grid(grid), _del(), _ex(ex) {

}

Interpolator::Interpolator(const std::vector<Vertex>& vertices,
			   const bool ex,
			   const std::string algo) :
  _dim(vertices[0].size()), _algo(algo), _grid(), _del(vertices), _ex(ex) {
}

Interpolator::Interpolator(const Interpolator& inter) {
  *this = inter;
}

Interpolator& Interpolator::operator=(const Interpolator& inter) {
  if ( this == &inter )
      return *this;

  _dim = inter._dim;
  _algo = inter._algo;
  _grid = inter._grid;
  _del = inter._del;
  _ex = inter._ex;

  return *this;
}

Interpolator::~Interpolator () {}

double Interpolator::Evaluate(const std::vector<double>& v) const {

  try {
    if ( _algo == "linear" ) {
      // Return 0 if the extrapolation is not requested and the point is outside of the grid
      if ( !_ex && !_grid.IsInside(v) )
	  return 0.;

      // Get the cell surrounding the point
      Grid cell = _grid.GetBestCell(v);

      // Feed it to the linear interpolation algorithm
      LinearInterpolator linint(cell);
      return linint(v);

    } else if ( _algo == "simplex" ) {
      // Find the facet in which the point lies or the closest simplex if extrapolation
      bool isin;
      std::vector<Vertex> vertices = _del.GetBestFacet(v, &isin);

      // Return 0 if the extrapolation is not requested and the point is outside of the hull
      if ( !_ex && !isin )
	  return 0.;

      // Feed the facet to the simplex interpolation algorithm
      SimplexInterpolator simpint(vertices);
      return simpint(v);

    }
    throw(Exceptions::Exception(Exceptions::recoverable,
	  "Algorithm is not supported: "+_algo,
	  "Interpolator::Evaluate"));

  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::recoverable,
	  "Could not evalute in the requested point"+std::string(e.what()),
	  "Interpolator::Evaluate"));
  }

  return 0.;
}

double Interpolator::Evaluate(const double& v) const {

  return Evaluate(std::vector<double>({v}));
}


double Interpolator::Evaluate(const double* v) const {

  std::vector<double> vv(v, v+_dim);
  return Evaluate(vv);
}

std::vector<TPolyLine*> Interpolator::Meshing() const {

  try {
    if ( _algo == "linear" ) {
      return _grid.Polygons();

    } else if ( _algo == "simplex" ) {
      return _del.Polygons();

    }
    throw(Exceptions::Exception(Exceptions::recoverable,
	  "No meshing produced for this class of interpolator: "+_algo,
	  "Interpolator::Meshing"));

  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::recoverable,
	  "Could not return a meshing"+std::string(e.what()),
	  "Interpolator::Meshing"));
  }

  return std::vector<TPolyLine*>();  
}
