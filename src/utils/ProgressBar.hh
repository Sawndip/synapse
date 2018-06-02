#ifndef PROGRESSBAR_HH
#define PROGRESSBAR_HH

// C++ includes
#include <iostream>
#include <unistd.h>
#include <chrono>
#include <sys/ioctl.h>

// Additional includes
#include "Pitch.hh"

typedef std::chrono::time_point<std::chrono::system_clock> time_stamp;

/** @brief Provides a progess bar to be included in a for loop
 *
 * 	   It defines a progress bar for that adapts to the terminal width,
 *	   prints the rate of progress and the estimated time of arrival.
 */
class ProgressBar {
 public:

  /** @brief Default constructor, sets the output to Pitch::info */
  ProgressBar(Pitch::errorLevel level=Pitch::info);

  /** @brief Copy constructor */
  ProgressBar(const ProgressBar& pbar);

  /** @brief Equality operator */
  ProgressBar& operator=(const ProgressBar& pbar);

  /** @brief Destructor */
  ~ProgressBar();

  /** @brief Sets the start of execution time */
  void SetStartTime(time_stamp start)			{ _start = start; }

  /** @brief Gets the start of execution time */
  const time_stamp& GetStartTime()			{ return _start; }

  /** @brief Sets the current time */
  void SetCurrentTime(time_stamp current)		{ _current = current; }

  /** @brief Gets the current time */
  const time_stamp& GetCurrentTime()			{ return _current; }

  /** @brief Sets the width of the terminal */
  void SetTerminalWidth(size_t twidth)			{ _twidth = twidth; }

  /** @brief Gets the width of the terminal */
  const size_t& GetTerminalWidth()			{ return _twidth; }

  /** @brief Sets the verbosity level */
  void SetVerbosityLevel(Pitch::errorLevel level)	{ _level = level; }

  /** @brief Gets the verbosity level */
  const Pitch::errorLevel& GetVerbosityLevel()		{ return _level; }

  /** @brief Prints the current progress in the terminal 
   *
   *  @param	n	Current entry being processed (0->(N-1))
   *  @param	N	Total number of entries to process
   */
  void GetProgress(const size_t n, const size_t N);

 private:
  /** @brief Fetch the terminal width from the OS */
  void FetchTerminalWidth();

  /** @brief Prints the bar for a given percentage of progress
   *
   *  @param	pc	Percentage of progress
   *  @param	n	Number of entries processed
   */
  void PrintProgress(const size_t pc, const size_t n);

  time_stamp 		_start;		///< Start time
  time_stamp 		_current;	///< Current time
  size_t 		_twidth;	///< Terminal width
  Pitch::errorLevel	_level;		///< Verbosity level
};

#endif
