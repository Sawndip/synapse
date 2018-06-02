#ifndef GLOBALS_HH
#define GLOBALS_HH

// Cpp includes
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <cstdlib>
#include <regex>
#include <getopt.h>

// ROOT includes
#include "TSystem.h"
#include "TStyle.h"
#include "TVector3.h"
#include "TColor.h"

// Other includes
#include "Exception.hh"
#include "Pitch.hh"
#include "GeometryHandler.hh"

/** @brief Represents a single global variable and defines the methods to read it
 */
class Global {
 public:

  /** @brief Default constructor, initializes everything to 0 */
  Global() : _val(""), _des("")				{}

  /** @brief Proper constructor */
  Global(const std::string& val): _val(val), _des("")	{}

  /** @brief Copy constructor */
  Global(const Global& global)			{ *this = global; }

  /** @brief Equality operator */
  Global& operator=(const Global& global)	{ _val = global._val; 
						  _des = global._des;
						  return *this; }
  /** @brief Destructor */
  ~Global()					{}

  /** @brief Overloads the ostream operator */
  friend std::ostream& operator<<(std::ostream& os, const Global& global)
						{ os << global._val; return os; }

  /** @brief Returns the variable as a bool */
  bool AsBool()	const				{ return (bool)strtol(&_val[0], NULL, 10); }

  /** @brief Returns the variable as an int */
  int AsInt() const				{ return (int)strtol(&_val[0], NULL, 10); }

  /** @brief Returns the variable as a long int */
  long int AsLongInt() const			{ return strtol(&_val[0], NULL, 10); }

  /** @brief Returns the variable as an int in the requested base */
  long int AsBase(int base) const		{ return strtol(&_val[0], NULL, base); }

  /** @brief Returns the variable as a float */
  float AsFloat() const				{ return strtof(&_val[0], NULL); }

  /** @brief Returns the variable as a double */
  double AsDouble() const			{ return strtod(&_val[0], NULL); }

  /** @brief Sets the decription of the varibale */
  void SetDescription(const std::string& des)	{ _des = des; }

  /** @brief Returns the decription of the varibale */
  const std::string& Description() const	{ return _des; }

  /** @brief Returns the variable as an std::string */
  const std::string& AsString()	const		{ return _val; };

  /** @brief Overload the unary ! operator */
  bool operator !() const			{ return !AsBool(); }

  /** @brief Overload the bool cast operator */
  operator bool() const				{ return AsBool(); }

  /** @brief Overload the int cast operator */
  operator int() const				{ return AsInt(); }

  /** @brief Overload the long int cast operator */
  operator long int() const			{ return AsLongInt(); }

  /** @brief Overload the size_t cast operator */
  operator size_t() const			{ return (size_t)AsLongInt(); }

  /** @brief Overload the float cast operator */
  operator float() const			{ return AsFloat(); }

  /** @brief Overload the double cast operator */
  operator double() const			{ return AsDouble(); }

  /** @brief Overload the std::string cast operator */
  operator std::string() const			{ return _val; }

 private:

  std::string 		_val;		///< Global variable stored as an std::string
  std::string		_des;		///< Description of the global variable
};

/** @brief Imports, stores and manages the global variables.
 *
 *	   It imports and stores the global reconstruction variables in a dictionary.
 *	   The variables can then be exported as any type requested by the code.
 *
 *	   This is a singleton class. It is intialized once the first time it is
 *	   instentiated and is called by the other objects as is.
 */
class Globals {
 public:

  /** @brief Returns the single instance of the class, initialize if not yet done
   *
   *  @param	argc	Number of command line arguments
   *  @param	argv	List of command line arguments
   */
  static Globals& GetInstance(int argc = 0,
	  		      char** argv = NULL);

  /** @brief Returns the list of data files */
  const std::vector<std::string>& GetDataFiles() const	{ return _data_files; }

  /** @brief Returns the geometry handelr */
  const GeometryHandler& GetGeometryHandler() const	{ return _geoh; }

  /** @brief Returns the requested global variable */
  const Global& operator()(const std::string& var) const;

  /** @brief Returns the requested global variable */
  const Global& operator[](const std::string& var) const;

  /** @brief Print the list of datacards and their value */
  void Print(Pitch::errorLevel level=Pitch::info) const;

 private:

  /** @brief Command line construct, takes extra arguments from the command line 
   *
   *  @param	argc	Number of command line arguments
   *  @param	argv	List of command line arguments
   */
  Globals(int argc,
	  char** argv);

  /** @brief Copy constructor, not implemented for singleton */
  Globals(const Globals& globals);

  /** @brief Equality operator, not implemented for singleton */
  Globals& operator=(const Globals& globals);

  /** @brief Destructor */
  ~Globals();

  /** @brief Loads the global constants from a datacards file */
  void LoadFile(std::string filename);

  /** @brief Load the extra arguments from the command line arguments
   *
   *  @param	argc	Number of command line arguments
   *  @param	argv	List of command line arguments
   */
  void LoadCommandLine(int argc,
		       char** argv);

  /** @brief Loads the geometry handler */
  void LoadGeometry();

  /** @brief Sets the default style of the ROOT plots */
  void SetDefaultStyle();

  std::map<std::string, Global>		_globals;	///< Dictionary of global variables
  std::vector<std::string>		_data_files;	///< List of runs (remaining arguments)
  GeometryHandler			_geoh;		///< Geometry handler
};

#endif
