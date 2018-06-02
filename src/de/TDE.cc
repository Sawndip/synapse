#include "TDE.hh"

TDE::TDE() :
  _vor(), _int(), _weights(), _interp() {
}

TDE::TDE(std::vector<std::vector<double>> vertices,
	 const bool interp,
	 const bool bounded,
	 const double alpha,
	 const double eps) :
  _vor(), _int(), _weights(), _interp(interp) {

  try {
    // Check for duplicates in the point sample, remove duplicates and retain weight (TODO)
    size_t weight;
    size_t i, j;
    for (i = 0; i < vertices.size(); i++) {
      weight = 1;
      for (j = i+1; j < vertices.size(); j++)
        if ( vertices[i] == vertices[j] ) {
	  weight++;
	  vertices.erase(vertices.begin()+j);
        }
      _weights.push_back(weight);
    }

    _vor = Voronoi(vertices, bounded, alpha, eps);
    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "TDE::TDE"));
  }
}

TDE::TDE(const TDE& tde) {
  *this = tde;
}

TDE& TDE::operator=(const TDE& tde) {
  if ( this == &tde )
      return *this;

  _vor = tde._vor;
  _int = tde._int;
  _weights = tde._weights;
  _interp = tde._interp;

  return *this;
}

TDE::~TDE () {}

void TDE::Initialize() {

  // If the interpolation is needed, estimate the density at each of the vertices,
  // feed the array to a simplexial interpolator
  if ( _interp ) {
    std::vector<Vertex> vertices;
    size_t n = _vor.GetNpoints();
    size_t i;
    for (i = 0; i < n; i++) {
      vertices.push_back(Vertex(_vor.GetPoint(i)));
      if ( _vor.GetCell(i).GetVolume() > 0 ) {
	  vertices[i].SetValue(_weights[i]/_vor.GetCell(i).GetVolume()/n);
      } else {
	  vertices[i].SetValue(0.);
      }
    }

    _int = Interpolator(vertices);
  }
}

double TDE::Evaluate(const std::vector<double>& v) const {

  // If interpolated, leave it to the interpolator
  if ( _interp )
	return _int(v);

  // If not, find the closest vertex and return its Voronoi density
  size_t bestid = _vor.GetBestCellID(v);
  if ( _vor.GetCell(bestid).GetVolume() > 0 )
      return _weights[bestid]/_vor.GetCell(bestid).GetVolume()/_vor.GetNpoints();

  return 0.;
}

double TDE::Evaluate(const double& v) const {

  return Evaluate(std::vector<double>({v}));
}

double TDE::Evaluate(const double* v) const {

  std::vector<double> vv(v, v+_vor.GetDimension());
  return Evaluate(vv);
}

std::vector<TPolyLine*> TDE::Meshing() const {

  // Throw if there is no interpolation on this estimate
  std::vector<TPolyLine*> meshing;
  if ( _interp )
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	    "No interpolation mesh produced for this TDE",
	    "TDE::Meshing"));

  return _int.Meshing();
}
