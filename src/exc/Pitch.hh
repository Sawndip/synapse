#ifndef PITCH_HH
#define PITCH_HH

// C++ includes
#include <ostream>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <sstream>

// Other includes
#include "Exception.hh"
#include "ExceptionLevel.hh"

/** @brief Redirect output to std::out, std::cerr or file depending on severity
 *
 * 	   Usage:
 * 	   Pitch::mout(Pitch::debug) << "Blah blah" << std::endl;
 *
 * 	   Pitch is used to redirect output at run time to some output file or
 * 	   standard outputs. Redirection is controlled by errorLevel:
 *   	     debug - debugging information. Deliver information related to debugging
 *           	     errors in the code; verbose developer level output
 *   	     info - run control information. Deliver information related to run
 *          	    information; aimed at users to give summary information about the
 *          	    run (idea is to give a constant, but limited, stream of data to
 *          	    info)
 *   	     warning - information to warn the user of suspected errors in running
 *   	     error - information to warn the user of certain errors in running
 *   	     fatal - final information to tell the user if the code is crashing
 *
 * 	   We also redirect output according to exceptionLevel
 *   	     recoverable - redirects to error
 *   	     nonRecoverable - redirects to fatal
 *
 * 	   Pitch is a singleton class - so initialisation is done whenever a function
 * 	   is called and the thing never destructs
 */
class Pitch {
 public:
  /** @brief Error level that determines where the output goes */
  enum errorLevel {debug, info, warning, error, fatal, log};

  /** @brief Shorthand that returns output to Pitch::mout(Pitch::debug) */
  static std::ostream & mout();

  /** @brief Returns a reference to the ostream corresponding to the errorLevel */
  static std::ostream & mout(errorLevel level);

  /** @brief Returns a reference to the ostream corresponding to the exceptionLevel */
  static std::ostream & mout(Exceptions::exceptionLevel level);

  /** @brief Set the ostream for a given error level (Exceptions::Exception) */
  static void setAnOutput(errorLevel level, std::ostream& out);

  /** @brief Set the ostream for all items below "verboseLevel" to /dev/null
   *  or the log file (if the logLevel is set to 2)
   *  If verboseLevel is less than or equal to:
   *   - Pitch::debug then mout(Pitch::debug) redirects to std::cout
   *   - Pitch::info then mout(Pitch::info) redirects to std::clog
   *   - Pitch::warning then mout(Pitch::warning) redirects to std::cerr
   *   - Pitch::error then mout(Pitch::error) redirects to std::cerr
   *   - Pitch::fatal then mout(Pitch::fatal) redirects to std::cerr
   *  Note that the redirection is independent of setStandardOutputs status
   */
  static void setOutputs(int verboseLevel, int logLevel);

  /** @brief Turn on/off std::cout, std::cerr, std::clog
   * Redirect the standard outputs to /dev/null or log file depending on verboseLevel and logLevel:
   * - if verboseLevel > int(Pitch::debug), set std::cout to /dev/null or log file
   * - if verboseLevel > int(Pitch::info), set std::clog to /dev/null or log file
   * - if verboseLevel > int(Pitch::warning), set std::cerr to /dev/null or log file
   */
  static void setStandardOutputs(int verboseLevel);

  /** @brief Set up the log file stream
   *  - if logLevel = 0 then Pitch::log will go to /dev/null
   *  - if logLevel = 1 file opened and Pitch::log will go there
   *  - if logLevel = 2 file opened, Pitch::log and any other streams not going to screen go there
   */
  static void setLog(int logLevel);

  /** @brief Close the log file and free the memory */
  static void closeLog();

  /** @brief Return the log file name */
  static std::string getLogName() 		{ return logname; }

  /** @brief Set the log file name */
  static void setLogName(std::string aLogname)	{ logname = aLogname; }

  /** @brief Activate std::cout
   *  	    If isActive is false, redirects std::cout to /dev/null; if isActive is
   *  	    true, redirects it to std::cout
   */
  static void activateCout(bool isActive);

  /** @brief Return false if std::cout points to voidout; else true */
  static bool coutIsActive();

  /** @brief Activate std::cerr
   *  	    If isActive is false, redirects std::cerr to /dev/null; if isActive is
   *  	    true, redirects it to std::cerr
   */
  static void activateCerr(bool isActive);

  /** @brief Return false if std::cerr points to voidout; else true */
  static bool cerrIsActive();

  /** @brief Activate std::clog
   *  	    If isActive is false, redirects std::clog to /dev/null; if isActive is
   *  	    true, redirects it to std::clog
   */
  static void activateClog(bool isActive);

  /** @brief Return false if std::clog points to voidout; else true */
  static bool clogIsActive();

  /** @brief Return ostream pointing to /dev/null (irrespective of verbose level) */
  static std::ostream& nullOut();

  /** @brief Return ostream pointing to standard output (irrespective of verb level) */
  static std::ostream& coutOut();

  /** @brief Return ostream pointing to standard log (irrespective of verbose level) */
  static std::ostream& clogOut();

  /** @brief Return ostream pointing to standard error (irrespective of verb level) */
  static std::ostream& cerrOut();

  /** @brief Return the level of the error in plain english */
  static std::string levelStr(errorLevel level);

  /** @brief Prints a message with a certain severity level 
   *  @param	level	Level of severity of the message
   *  @param	message	Message to print out along with the erreor
   */
  static void print(errorLevel level, std::string message);

  /** @brief Prints a message in a given function with a certain severity level 
   *  @param	level	Level of severity of the message
   *  @param	message	Message to print out along with the erreor
   *  @param	func	Function in which the error was encountered
   */
  static void print(errorLevel level, std::string message, std::string func);

 private:
  /** @brief Constructor is called by any call to mout. Defines std::map output*/
  Pitch();

  /** @brief Ensures class is always instantiated whenever needed.
   *  @returns Pointer to the singleton class instance
   */
  static Pitch* getInstance();

  /** @brief Setup outputs */
  static void initialiseOutputs();

  /** @brief Set Exception */
  static void setException();

  /* Static variables */
  static std::map<errorLevel,
		std::ostream*>	output;			///< Links an ostream with each error level
  static std::ostream* 		voidout;		///< Dummy ostream to send stuff to /dev/null
  static std::ostream* 		stdout;			///< ostream that points to std::out
  static std::ostream* 		stdlog;			///< ostream that points to std::log
  static std::ostream* 		stderr; 		///< ostream that points to std::err
  static std::ofstream* 	logfile;		///< ostream that points to a log file
  static const errorLevel 	default_error_level;	///< The default error level
  static const int 		default_log_level;	///< The default log level
  static std::string 		logname;		///< The log file name
  static Pitch* 		instance;		///< Pointer to the singleton instance
};

#endif
