template<typename T> void Trim(std::vector<T>& vx,
			       const double& frac) {

  // If the sample is empty, throw
  Assert::IsNotEmpty("Math::Trim", "Input sample", vx);

  // Handle the special values of n
  if ( !frac ) {
    return;
  } else if ( frac == .5 ) {
    vx = {Median(vx)};
    return;
  } else if ( frac < 0. || frac > .5 ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	               		"Invalid trimming fraction (0 < frac < .5)",
				"Math::Trim"));
  }

  // Sort the points and discard the smallest and largest values
  std::vector<double> vnew;
  std::sort(vx.begin(), vx.end());
  size_t ntrim = frac*vx.size();
  size_t nfinal = vx.size()-2*ntrim;
  size_t i;
  for (i = ntrim; i < ntrim+nfinal; i++)
      vnew.push_back(vx[i]);

  // Randomly shuffle the elements to not return a sorted vector
  std::srand(std::time(NULL));
  std::random_shuffle(vnew.begin(), vnew.end());
  vx = vnew;
}

template<typename T> std::vector<T> Trimmed(std::vector<T> vx,
			       		    const double& frac) {

  Trim(vx, frac);
  return vx;
}

template<typename T> void Winsorize(std::vector<T>& vx,
			       	    const double& frac) {

  // If the sample is empty, throw
  Assert::IsNotEmpty("Math::Winsorize", "Input sample", vx);

  // Handle the special values of n
  if ( !frac ) {
    return;
  } else if ( frac == .5 ) {
    double med = Median(vx);
    for (size_t i = 0; i < vx.size(); i++)
	vx[i] = med;
    return;
  } else if ( frac < 0. || frac > .5 ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	               		"Invalid fraction (0 < frac < .5)",
				"Math::Winsorize"));
  }

  // Sort the points and replace the smallest and largest values
  std::sort(vx.begin(), vx.end());
  size_t id = frac*vx.size();
  double lbound(vx[id]), ubound(vx[vx.size()-1-id]);
  size_t i;
  for (i = 0; i < id; i++) {
    vx[i] = lbound;
    vx[vx.size()-1-i] = ubound;
  }

  // Randomly shuffle the elements to not return a sorted vector
  std::srand(std::time(NULL));
  std::random_shuffle(vx.begin(), vx.end());
}

template<typename T> std::vector<T> Winsorized(std::vector<T> vx,
			       	    	       const double& frac) {

  Winsorize(vx, frac);
  return vx;
}

template<typename I> I RandomElement(I begin, I end) {

    const size_t n = std::distance(begin, end);
    const size_t divisor = RAND_MAX / n;

    unsigned long k;
    do { k = std::rand() / divisor; } while (k >= n);

    std::advance(begin, k);
    return begin;
}

template<typename T> std::vector<T> Resample(const std::vector<T>& vx,
					     size_t m,
					     const bool rep) {

  // If the sample is empty, throw
  Assert::IsNotEmpty("Math::Resample", "Input sample", vx);

  // If m = 0, default it to the sample size
  if ( !m )
      m = vx.size();

  // If no replacement, m must be smaller than n. If not, throw.
  if ( !rep && m >= vx.size() )
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	               		  "Not enough points to resample without replacement",
				  "Math::Resample"));

  // Randomly pick m samples from the input sample. If no replacement, ditch (n-m) points randomly
  std::vector<T> subsample;
  size_t i, pos, divisor;
  std::srand(time(NULL));
  if ( rep ) {
    subsample.resize(m);
    divisor = RAND_MAX/vx.size();
    for (i = 0; i < m; i++) {
      do {
	pos = std::rand()/divisor;
      } while ( pos > vx.size() );
      subsample[i] = vx[pos];
    }
  } else {
    subsample = vx;
    std::random_shuffle(subsample.begin(), subsample.end());
    subsample.resize(m);
  }

  return subsample;
}
