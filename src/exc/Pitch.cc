#include "Pitch.hh"

// Initialise the static member variables
std::ostream* Pitch::stdout    = NULL;
std::ostream* Pitch::stdlog    = NULL;
std::ostream* Pitch::stderr    = NULL;
std::ostream* Pitch::voidout   = NULL;
std::ofstream* Pitch::logfile  = NULL;
Pitch* Pitch::instance        = NULL;

const Pitch::errorLevel Pitch::default_error_level = Pitch::warning;
const int Pitch::default_log_level = 0;
std::string Pitch::logname = "emittance.log";

std::map<Pitch::errorLevel, std::ostream*> Pitch::output;

Pitch::Pitch() {
  // Do nothing
}

std::ostream & Pitch::mout() {
  getInstance();
  return Pitch::mout(debug);
}

std::ostream & Pitch::mout(errorLevel level) {
  getInstance();
  return (*output[level]);
}

std::ostream & Pitch::mout(Exceptions::exceptionLevel level) {
  getInstance();
  // Eventually a map or somesuch I suspect
  if (level == Exceptions::recoverable)
    return (*output[Pitch::error]);
  else
    return (*output[Pitch::fatal]);
}

void  Pitch::setAnOutput(errorLevel level, std::ostream& out) {
  getInstance();
  output[level] = &out;
}

Pitch * Pitch::getInstance() {
  if (instance  == NULL) {
      instance = new Pitch();
      initialiseOutputs();
      setOutputs(default_error_level, default_log_level);
      setStandardOutputs(default_error_level);
  }
  return instance;
}

void Pitch::setOutputs(int verboseLevel, int logLevel) {
  getInstance();
  setLog(logLevel);
  std::ostream* out[5] = {stdout, stdlog, stderr, stderr, stderr};
  for (int i = 0; i <= fatal; ++i) {
    if (i >= verboseLevel)
      output[Pitch::errorLevel(i)] = out[i];
    else if (logLevel > 1)
      output[Pitch::errorLevel(i)] = logfile;
    else
      output[Pitch::errorLevel(i)] = voidout;
  }
}

void Pitch::setLog(int logLevel) {
  getInstance();
  // At present, any log level above 0 just goes to a standard file
  if (logLevel) {
    if (!logfile) {
      logfile = new std::ofstream(logname.c_str());
    } else if (!logfile->is_open()) {
      logfile->open(logname.c_str());
    }
    output[Pitch::log] = logfile;
  } else {
    output[Pitch::log] = voidout;
  }
}

void Pitch::closeLog() {
  if (logfile) {
    if (logfile->is_open())
      logfile->close();
    delete logfile;
    logfile = NULL;
  }
}

void Pitch::setStandardOutputs(int verboseLevel) {
  getInstance();
  activateCout(verboseLevel <= Pitch::debug);
  activateClog(verboseLevel <= Pitch::info);
  activateCerr(verboseLevel <= Pitch::error);
}

void Pitch::activateCout(bool isActive) {
  getInstance();
  if (isActive)
    std::cout.rdbuf(stdout->rdbuf());
  else
    std::cout.rdbuf(voidout->rdbuf());
}

bool Pitch::coutIsActive() {
  getInstance();
  return std::cout.rdbuf() != voidout->rdbuf();
}

void Pitch::activateCerr(bool isActive) {
  getInstance();
  if (isActive)
    std::cerr.rdbuf(stdlog->rdbuf());
  else
    std::cerr.rdbuf(voidout->rdbuf());
}

bool Pitch::cerrIsActive() {
  getInstance();
  return std::cerr.rdbuf() != voidout->rdbuf();
}

void Pitch::activateClog(bool isActive) {
  getInstance();
  if (isActive)
    std::clog.rdbuf(stdlog->rdbuf());
  else
    std::clog.rdbuf(voidout->rdbuf());
}

bool Pitch::clogIsActive() {
  getInstance();
  return std::clog.rdbuf() != voidout->rdbuf();
}

void Pitch::initialiseOutputs() {
  // assume they are all uninitialised if one is (we always initialise together)
  voidout = new std::ofstream();  // this points at /dev/null by default
  stdout  = new std::ofstream();
  stdout->rdbuf(std::cout.rdbuf());
  stdlog  = new std::ofstream();
  stdlog->rdbuf(std::clog.rdbuf());
  stderr  = new std::ofstream();
  stderr->rdbuf(std::cerr.rdbuf());
}

std::ostream& Pitch::nullOut() {
  getInstance();
  return *voidout;
}

std::ostream& Pitch::coutOut() {
  getInstance();
  return *stdout;
}

std::ostream& Pitch::clogOut() {
  getInstance();
  return *stdlog;
}

std::ostream& Pitch::cerrOut() {
  getInstance();
  return *stderr;
}

std::string Pitch::levelStr(errorLevel level) {

  switch (level) {
      case debug : return "DEBUG";
      case info : return "INFO";
      case warning : return "WARNING";
      case error : return "ERROR";
      case fatal : return "FATAL";
      case log : return "LOG";
  }

  return "UNKNOWN";
}

void Pitch::print(errorLevel level, std::string message) {

  mout(level) << "[" << levelStr(level) << "] " << message << std::endl;
}

void Pitch::print(errorLevel level, std::string message, std::string func) {

  mout(level) << "[" << levelStr(level) << "][" << func << "] " << message << std::endl;
}
