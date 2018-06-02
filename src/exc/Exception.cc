#ifndef NO_STACKTRACE
#include <stdlib.h>
#include <execinfo.h>
#endif

#include "Pitch.hh"
#include "Exception.hh"

namespace Exceptions {

const size_t   Exception::_maxStackSize = 100;
bool Exception::_willDoStackTrace = true;

Exception::Exception() throw() {
  Pitch::mout();
}

Exception::Exception(exceptionLevel level,
		     std::string errorMessage,
		     std::string location) throw()
 : exception(), _message(errorMessage), _location(location),
   _stacktrace(""), _what(), _level(level) {

  // Make sure we initialise Pitch otherwise can get segv
  Pitch::mout();

  // Use std::vector as dynamic array
  SetWhat("\n["+_location+"] "+_message);
  if (_willDoStackTrace)
      _stacktrace = MakeStackTrace(2);
}

void Exception::Print() {
  Pitch::mout(_level) << _message << "\n";
  Pitch::mout(Pitch::debug) << "Error at " << _location << "\n";
  if (_willDoStackTrace && _stacktrace != "")
    Pitch::mout(Pitch::debug) << "Stack trace\n" << GetStackTrace() << "\n";
}

void Exception::SetWillDoStackTrace(bool willDoStackTrace) {
    _willDoStackTrace = willDoStackTrace;
}

bool Exception::GetWillDoStackTrace() {
    return _willDoStackTrace;
}

#ifndef NO_STACKTRACE
std::string Exception::MakeStackTrace(size_t skipTrace) {

  //
  size_t stackSize;
  void * stackAddress[_maxStackSize];
  char **stackNames;

  std::stringstream sstr;
  stackSize  = backtrace(stackAddress, _maxStackSize);
  stackNames = backtrace_symbols(stackAddress, stackSize);

  for (size_t i = skipTrace; i < stackSize; i++)
      sstr << stackNames[i] << std::endl;
  free(stackNames);
  return sstr.str();
}
#endif

#ifdef NO_STACKTRACE
std::string Exception::MakeStackTrace(size_t skipTrace) {return "";}
#endif
}
