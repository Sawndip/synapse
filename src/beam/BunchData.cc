#include "BunchData.hh"

namespace Beam {

BunchData::BunchData() :
  _trans(), _etrans(), _S(), _opt(), _eopt(), _T(), _pz(), _epz(), _axes(), _keys(), _ids() {
}

BunchData::BunchData(const BunchMap& samples,
  	     	     const BunchMap& errors) :
  _trans(), _etrans(), _S(), _opt(), _eopt(), _T(), _pz(), _epz(), _axes(), _keys(), _ids() {

  try {
    Initialize(samples, errors);
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "BunchData::BunchData"));
  }
}

BunchData::BunchData(const BunchData& bunchdata) {
  *this = bunchdata;
}

BunchData& BunchData::operator=(const BunchData& bunchdata) {
  if ( this == &bunchdata )
      return *this;

  _trans = bunchdata._trans;
  _etrans = bunchdata._etrans;
  _S = bunchdata._S;
  _opt = bunchdata._opt;
  _eopt = bunchdata._eopt;
  _T = bunchdata._T;
  _pz = bunchdata._pz;
  _epz = bunchdata._epz;
  _axes = bunchdata._axes;
  _keys = bunchdata._keys;
  _ids = bunchdata._ids;

  return *this;
}

BunchData::~BunchData () {}

const std::vector<double>& BunchData::Samples(const std::string var) const {

  if ( var == "pz" )
      return _pz;

  return _trans[_ids.at(var)];
}

const std::vector<double>& BunchData::Errors(const std::string var) const {

  if ( var == "pz" )
      return _epz;

  return _etrans[_ids.at(var)];
}

const double& BunchData::S(const std::string vara, const std::string varb) const {

  return _S[_ids.at(vara)][_ids.at(varb)];
}

const std::vector<double>& BunchData::Opticals(const std::string var) const {

  return _opt[_ids.at(var)];
}

const std::vector<double>& BunchData::OpticalErrors(const std::string var) const {

  return _eopt[_ids.at(var)];
}

const double& BunchData::T(const std::string vara, const std::string varb) const {

  return _T[_ids.at(vara)][_ids.at(varb)];
}

bool BunchData::Contains(const std::string& var) const {

  if ( var == "pz" )
      return true;

  return _ids.find(var) != _ids.end();
}

void BunchData::AssertContains(const std::string& var) const {

  if ( !Contains(var) )
      throw(Exceptions::Exception(Exceptions::recoverable,
	    "Requested variable not provided with samples: "+var,
	    "BunchData::AssertContains"));
}

void BunchData::Initialize(const BunchMap& samples,
			   const BunchMap& errors) {

  // First of all, extract the transverse phase space
  for (const std::pair<std::string, std::vector<double>>& row : samples) {
    if ( row.first == "pz" ) {
      _pz = row.second;
      continue;
    }
    _ids[row.first] = _keys.size();
    _keys.push_back(row.first);
  }
  size_t dim = _keys.size();

  // If no sample, empty samples or inconsistent samples, abort
  Assert::IsNonZero("BunchData::Initialize", "Dimension", dim);
  Assert::IsNotEmpty("BunchData::Initialize", "Input particle sample", samples.at(_keys[0]));
  size_t N = samples.at(_keys[0]).size();
  for (const std::string& key : _keys) {
    if ( samples.at(key).size() != N ) {
      std::ostringstream message;
      message << "samples[" << _keys[0] << "] and samples[" << key << "]";
      Assert::SameSize("BunchData::Initialize", message.str(), samples.at(_keys[0]), samples.at(key));
      return;
    }
  }

  // Get the axes from the list of keys
  if ( dim == 2 ) {
    std::string axis = _ids.find("x") != _ids.end() ? "x" : "y";
    _axes.push_back(axis);
  } else if ( dim == 4 ) {
    _axes = {"x", "y"};
  } else {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Number of dimensions not currently supported: "+std::to_string(dim),
	  "BunchData::Initialize"));
    return;
  }

  // Initialize the sample matrix, exclude the longitudinal momentum
  size_t i;
  _trans = Matrix<double>(dim, N);
  for (i = 0; i < dim; i++)
      _trans[i] = samples.at(_keys[i]);

  // Fill the matrix of errors. If the measurement errors on the phase-space vectors
  // are not provided, focus on the statistical errors only
  if ( errors.size() ) {
    size_t N = errors.at(_keys[0]).size();
    _etrans = Matrix<double>(dim, N);
    for (i = 0; i < _keys.size(); i++)
        _etrans[i] = errors.at(_keys[i]);
    if ( errors.find("pz") != errors.end() )
        _epz = errors.at("pz");
  }

  // Initialize the phase space covariance matrix
  _S = Math::CovarianceMatrix(_trans);

  // Initialize the optical sample
  _opt = _trans;
  for (const std::string& axis : _axes)
    for (i = 0; i < N; i++) {
      Assert::IsNonZero("BunchData::Initialize", "pz["+std::to_string(i)+"]", _pz[i]);
      _opt[_ids["p"+axis]][i] /= _pz[i];
    }

  // Initialize the optical sample errors
  _eopt = _etrans;
  for (const std::string& axis : _axes)
    for (i = 0; i < _eopt.Ncols(); i++)
      _eopt[_ids["p"+axis]][i] =
	sqrt(pow(_pz[i]*_etrans[_ids["p"+axis]][i], 2) + 
	     pow(_trans[_ids["p"+axis]][i]*_epz[i], 2))/pow(_pz[i], 2);

  // Initialize the trace space covariance matrix
  _T = Math::CovarianceMatrix(_opt);
}

} // namespace Beam
