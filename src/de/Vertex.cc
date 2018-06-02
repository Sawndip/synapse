#include "Vertex.hh"

Vertex::Vertex() :
  _coords(0), _value(0.), _indices(0) {
}

Vertex::Vertex(const std::vector<double>& coords,
	       const double& value,
	       const std::vector<size_t>& indices) :
  _coords(coords), _value(value), _indices(indices) {
}

Vertex::Vertex(const size_t n) :
  _coords(std::vector<double>(n)), _value(0.), _indices(0) {
}


Vertex::Vertex(const Vertex& vertex) {
  *this = vertex;
}

Vertex& Vertex::operator=(const Vertex& vertex) {
  if ( this == &vertex )
      return *this;

  _coords = vertex._coords;
  _value = vertex._value;
  _indices = vertex._indices;

  return *this;
}

Vertex::~Vertex () {}

double Vertex::GetNorm() const {

  double sum(0);
  for (const double x : _coords )
      sum += x*x;

  return sqrt(sum);
}
