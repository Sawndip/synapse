#include "ProgressBar.hh"

ProgressBar::ProgressBar(Pitch::errorLevel level) :
  _start(std::chrono::system_clock::now()),
  _current(std::chrono::system_clock::now()),
  _twidth(0), _level(level) {

  FetchTerminalWidth();
  PrintProgress(0, 0);
}

ProgressBar::ProgressBar(const ProgressBar& pbar) {
  *this = pbar;
}

ProgressBar& ProgressBar::operator=(const ProgressBar& pbar) {
  if ( this == &pbar )
      return *this;

  this->SetStartTime(pbar._start);
  this->SetCurrentTime(pbar._current);
  this->SetTerminalWidth(pbar._twidth);
  this->SetVerbosityLevel(pbar._level);

  return *this;
}

ProgressBar::~ProgressBar () {}

void ProgressBar::GetProgress(size_t n, size_t N) {

  // Progress of the entire task in percent (integer lower bound)
  size_t pc = 100*(n+1)/N;

  // Only print the progress bar if it has incremented since the previous entry
  size_t pc_1 = 100*n/N;
  if ( pc > pc_1 )
      PrintProgress(pc, n);
}

void ProgressBar::FetchTerminalWidth() {
  struct winsize size;
  ioctl(STDERR_FILENO, TIOCGWINSZ, &size);
  _twidth = size.ws_col;
}

void ProgressBar::PrintProgress(const size_t pc, const size_t n) {

  // Get the current time to compare it with the start time (integer lower bound)
  _current = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed = _current-_start;
  int elapsed_s = (int)elapsed.count();
  int elapsed_m = elapsed_s/60;
  elapsed_s %= 60;

  // Calculate the amount of entries process per second (integer lower bound)
  size_t eps = n/elapsed.count();

  // Infer how long it will take until it's finished (ETA)
  int left_s(0);
  if ( pc > 0 )
    left_s = elapsed.count()*(100.-pc)/pc;
  int left_m = left_s/60;
  left_s %= 60;

  // Define and fill a string with the essential information first
  std::string info_string;
  if ( pc < 100 ) {
    info_string = " "+std::to_string(pc)+"%    "+std::to_string(eps)+"e/s    ETA "+
				std::to_string(left_m)+"m "+std::to_string(left_s)+"s\r";
  } else {
    info_string = " "+std::to_string(pc)+"%    "+std::to_string(eps)+"e/s    in "+
		    		std::to_string(elapsed_m)+"m "+std::to_string(elapsed_s)+"s\r";
  }

  // Fill the rest of the available space with a progression bar
  size_t ntaken = info_string.size();
  size_t navail = _twidth-ntaken-2;
  std::string bar_string;
  bar_string.append("[");
  for (size_t k = 0; k < navail; k++) {
    if ( k < pc*navail/100 ) bar_string.append("=");
    else if ( k == pc*navail/100 ) bar_string.append(">");
    else if ( k > pc*navail/100 ) bar_string.append(" ");
  }
  bar_string.append("]");

  // Print the string to the terminal
  Pitch::mout(_level) << bar_string+info_string;
  if ( pc < 100 ) {
    Pitch::mout(_level).flush();
  } else {
    Pitch::mout(_level) << std::endl;
  }
}
