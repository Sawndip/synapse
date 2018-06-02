#ifndef EXCEPTION_HH
#define EXCEPTION_HH

// C++ includes
#include <cstdio>
#include <cstring>
#include <exception>
#include <string>
#include <vector>
#include <ostream>
#include <sstream>

// Other includes
#include "ExceptionLevel.hh"

/** @brief Environment of custom exceptions to be thrown by the code
 */
namespace Exceptions {

/** @brief Global exception object for all non-standard exceptions thrown by the code.
 *
 * Exception has a severity (exceptionLevel), an error message and a
 * location. Also a function to return the stack trace (gcc only).
 * 
 * If not using gcc, define NO_STACKTRACE at compile time
 * 
 * Note that the throw() directives are there to explicitly declare that
 * exceptions cannot be thrown in the constructors and destructors
 */
class Exception : public std::exception {
 public:

  /** @brief Default constructor, initilizes everything to 0 */
  Exception() throw();

  /** @brief Normal constructor with error level, error message and location
   *
   *  @param	level		Severity of the exeption
   *  @param	errorMessage	Error message to throw
   *  @param	location	Location of the exception in the code
   */
  Exception(exceptionLevel level,
	    std::string errorMessage,
	    std::string location) throw();

  /** @brief Destructor */
  ~Exception() throw() {}

  /** @brief Return char buffer of  message+" at "+location, memory remains owned by Exception
   *
   */
  const char* what() const throw() 		{ return &_what[0]; }

  /** @brief Print the error message and location
   *
   *  	     Prints to error message to Pitch::mout(Pitch::error) and Location to
   *  	     Pitch::mout(Pitch::debug)
   */
  void Print();

  /** @brief Get the severity of the exception */
  inline exceptionLevel GetErrorLevel() const 	{ return _level; }

  /** @brief Get the error message (as defined by constructor) */
  inline std::string GetMessage() const 	{ return _message; }

  /** @brief Set the error message */
  inline void SetMessage(std::string new_message);

  /** @brief Get the location (as defined by Exception constructor) of the error */
  std::string GetLocation() const 		{ return _location; }

  /** @brief  Return the stack trace if it was stored */
  std::string GetStackTrace() const		{ return _stacktrace; }

  /** @brief Gcc-specific code to recover the stack trace as a string.
   *
   *  	     Will skip traces below skipTrace (so, for example, if we want to know
   *  	     where the Exception was thrown we might set skipTrace to 2, to skip
   *  	     MakeStackTrace and Exception constructor).
   */
  static std::string MakeStackTrace(size_t skipTrace);

  /** @brief Set that it makes a stack trace on exception */
  static void SetWillDoStackTrace(bool willDoStackTrace);

  /** @brief Return whether makes a stack trace on exception */
  static bool GetWillDoStackTrace();

 private:
  /** @brief Sets the reason for the exception */
  inline void SetWhat(std::string what_str);

  static const size_t 	_maxStackSize;		///< Maximum memory size of the stack
  static bool 		_willDoStackTrace;	///< Returns the stack trace if true
  std::string       	_message;		///< Message to be thrown
  std::string       	_location;		///< Name of the function that threw
  std::string       	_stacktrace;		///< Stack trace as an std::string
  std::vector<char> 	_what;			///< Reason for throwing
  exceptionLevel    	_level;			///< Severity of the exception
};

void Exception::SetMessage(std::string new_message) {
    _message = new_message;
    SetWhat("\n["+_location+"] "+_message);
}

void Exception::SetWhat(std::string what_str) {
    _what = std::vector<char>(what_str.size()+1);
    snprintf(&_what[0], what_str.size()+1, "%s", what_str.c_str());
}

} // namespace Exceptions

#endif
