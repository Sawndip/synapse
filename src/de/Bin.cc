#include "Bin.hh"

Bin::Bin() :
  _dim(0), _indices(0), _lower(0), _upper(0), _content(0.) {
}

Bin::Bin(const std::vector<size_t>& indices,
      	 const std::vector<double>& lower,
      	 const std::vector<double>& upper,
      	 const double& content) :
  _dim(indices.size()), _indices(indices), _lower(lower), _upper(upper), _content(content) {
}

Bin::Bin(const Bin& bin) {
  *this = bin;
}

Bin& Bin::operator=(const Bin& bin) {
  if ( this == &bin )
      return *this;

  _dim = bin._dim;
  _indices = bin._indices;
  _lower = bin._lower;
  _upper = bin._upper;
  _content = bin._content;

  return *this;
}

Bin::~Bin () {}

std::vector<double> Bin::GetBarycentre() const {

  std::vector<double> barycentre;
  size_t i;
  for (i = 0; i < _dim; i++)
      barycentre.push_back(GetCentre(i));

  return barycentre;
}

double Bin::GetCentre(const size_t i) const {

  return (_lower[i]+_upper[i])/2.;
}

double Bin::GetVolume() const {

  double vol(1.);
  size_t i;
  for (i = 0; i < _dim; i++)
      vol *= GetWidth(i);

  return vol;
}

double Bin::GetWidth(const size_t i) const {

  return (_upper[i]-_lower[i]);
}
