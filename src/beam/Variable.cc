#include "Variable.hh"

namespace Beam {

Variable::Variable() :
  _val(0), _merr(0), _serr(0) {
}

Variable::Variable(const double& val,
	   	   const double merr,
	   	   const double serr) :
  _val(val), _merr(merr), _serr(serr) {
}

Variable::Variable(const Variable& var) {
  *this = var;
}

Variable& Variable::operator=(const Variable& var) {
  if ( this == &var )
      return *this;

  _val = var._val;
  _merr = var._merr;
  _serr = var._serr;

  return *this;
}

Variable::~Variable () {}

Variable::operator double() const {

  return _val;
}

std::ostream& operator<<(std::ostream& os, const Variable& var) {

  os << var.GetValue() << " +/- " << var.GetSError() << " (stat) +/- " 
				  << var.GetMError() << " (meas)";
  return os;
}
} // namespace Beam
