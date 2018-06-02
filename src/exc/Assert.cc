#include "Assert.hh"

namespace Assert {

void IsWithinBounds(const std::string& location,
		    const size_t n,
		    const size_t i) {

  if ( i >= n )
      throw(Exceptions::Exception(Exceptions::nonRecoverable, "Index out of bounds", location));
}

void IsProbability(const std::string& location,
		   const std::string& what,
		   const double& p,
		   const bool nonzero,
		   const bool noncert) {

  if ( !nonzero && !noncert && (p < 0. || p > 1.) ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  what+" does not represent a probability value",
	  location));
  } else if ( nonzero && !noncert && (p <= 0. || p > 1.) ) {
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	    what+" does not represent a non-zero probability value",
	    location));
  } else if ( !nonzero && noncert && (p < 0. || p >= 1.) ) {
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	    what+" does not represent a non-certain probability value",
	    location));
  } else if ( nonzero && noncert && (p <= 0. || p >= 1.) ) {
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	    what+" does not represent a non-extreme probability value",
	    location));
  }
}

} // namespace Assert
