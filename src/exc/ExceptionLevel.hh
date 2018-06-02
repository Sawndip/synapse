#ifndef EXCEPTIONLEVEL_HH
#define EXCEPTIONLEVEL_HH

/** @brief Environment of custom exceptions to be thrown by the code
 */
namespace Exceptions {

  /** @brief ExceptionLevel enumerates the severity of the exception.
   *
   *  	     An enumeration allows to distinguish between different error levels.
   *   	      - nonRecoverable means the internal state of the programme is no longer
   *            well-defined i.e. some memory problem.
   *  	      - recoverable means that in principle the code could keep on running, although
   *    	most of the time this results in end of run
   */
  enum exceptionLevel {recoverable, nonRecoverable};
}

#endif
