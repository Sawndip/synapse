#include "OptimalBinning.hh"

// Jackknife likelihood Minuit function used by the minimizer
static std::vector<double> lower, range;
static const std::vector<std::vector<double>>* data(NULL);
static void Jackknife(int& npar, double* gin, double& f, double* par, int flag) {

  // Set the alpha parameter, get the dimension of the data
  double alpha = 1.;
  size_t n = data->at(0).size();

  // Calculate the ranges, stepping and the total number of bins
  double vol = 1;
  std::vector<size_t> steps(n), number(n);
  std::vector<double> upper(n);
  double width;
  size_t ntotal = 1;
  size_t k;
  for (k = 0; k < (size_t)n; k++) {
    steps[k] = 1;
    number[k] = (size_t)(par[k]);
    width = range[k]/par[k];
    upper[k] = lower[k]+number[k]*width;
    vol *= width;
    ntotal *= number[k];
    if ( k > 0 )
	steps[k] = number[k-1]*steps[k-1];
  }

  // Get the sum of the current (N+N_{0}*...*N_{n-1}*alpha-1)
  double sumn = data->size() + ntotal*alpha - 1;

  // Bin the data
  Histogram hist(number, lower, upper);
  hist.FillN(*data);

  // Calculate and store the log likelihood
  double Nj;			// Number of entries in the current bin
  f = 0.;			// Log likelihood to increment
  for (const Bin& bin : hist.GetBinArray()) {
    Nj = bin.GetContent();
    if ( Nj )
        f += Nj*log((Nj+alpha-1.)/(vol*sumn));
  }

  f = -f;	// Return -ln(L)
}

static void CrossValidation(int& npar, double* gin, double& f, double* par, int flag) {

  // Get the amount of data and its dimension
  size_t N = data->size();
  size_t n = data->at(0).size();

  // Calculate the ranges, stepping and the total number of bins
  double vol = 1;
  std::vector<size_t> steps(n), number(n);
  std::vector<double> upper(n);
  double width;
  size_t ntotal = 1;
  size_t k;
  for (k = 0; k < (size_t)n; k++) {
    steps[k] = 1;
    number[k] = (size_t)(par[k]);
    width = range[k]/par[k];
    upper[k] = lower[k]+number[k]*width;
    vol *= width;
    ntotal *= number[k];
    if ( k > 0 )
	steps[k] = number[k-1]*steps[k-1];
  }

  // Get the first element of the function
  double first = 2./((N-1)*vol);

  // Bin the data
  Histogram hist(number, lower, upper);
  hist.FillN(*data);

  // Loop over the bins and increment the second element
  double second(0.);
  for (const Bin& bin : hist.GetBinArray())
      second += pow(bin.GetContent(), 2);
  second *= (N+1)/(N*N*(N-1)*vol);

  f = first-second;
}

OptimalBinning::OptimalBinning() :
  _dim(0), _hist(), _int(), _interp(false), _ex(false) {
}

OptimalBinning::OptimalBinning(const std::vector<std::vector<double>>& points,
			       const std::string algo,
			       const bool interp,
			       const bool ex) :
  _dim(0), _hist(), _int(), _interp(interp), _ex(ex) {

  try {
    // Save the point cloud at the global level for the optimizer to access it
    if ( !data )
        delete data;
    data = &points;

    // Initialize the optimal histogram
    Initialize(points, algo);
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "OptimalBinning::OptimalBinning"));
  }
}

OptimalBinning::OptimalBinning(const OptimalBinning& ob) {
  *this = ob;
}

OptimalBinning& OptimalBinning::operator=(const OptimalBinning& ob) {
  if ( this == &ob )
      return *this;

  _dim = ob._dim;
  _hist = ob._hist;
  _int = ob._int;
  _interp = ob._interp;
  _ex = ob._ex;

  return *this;
}

OptimalBinning::~OptimalBinning () {}

void OptimalBinning::Initialize(const std::vector<std::vector<double>>& points,
				const std::string& algo) {

  // If there are no points provided, no way to know what size to build the histo as
  Assert::IsNotEmpty("OptimalBinning::Initialize", "Input vector of point", points);
  Assert::IsNonZero("OptimalBinning::Initialize", "Point dimension", points[0].size());
  _dim = points[0].size();

  // Base the boundaries on the data:
  // Find the mean and spread of the sample to determine its boundaries
  // Using [mu-5*sigma, mu+5*sigma] contains most things
  std::vector<double> samples;
  std::vector<double> mean(_dim), std(_dim), upper(_dim);
  lower.resize(_dim);
  range.resize(_dim);
  size_t i, j;
  for (i = 0; i < _dim; i++) {
    samples.resize(0);
    for(j = 0; j < points.size(); j++)
	samples.push_back(points[j][i]);

    mean[i] = Math::Mean(samples);
    std[i] = Math::RMS(samples);
    lower[i] = mean[i] - 5*std[i];
    upper[i] = mean[i] + 5*std[i];
    range[i] = upper[i]-lower[i];
  }

  // Determine the optimal binning according to the chosen algorithm
  std::vector<size_t> number(_dim);
  if ( algo == "scott" ) {
    // In the normal distribution approximation, the width is fixed by Scott's normal ref. rule
    double delta = pow(24*sqrt(M_PI)/points.size(), 1./(_dim+2));
    for (i = 0; i < _dim; i++)
	number[i] = (size_t)(range[i]/delta/std[i]);

  } else if ( algo == "jackknife") {
    // Create a likelihood function to minimize to find the optimal binning
    // Algorithm developed by D.W. Hogg using a simple "jackknife" likelihood
    // Minimize the jackknife likelihood with TMinuit (faster)
    TMinuit* minimizer = new TMinuit(_dim);	// Specify the number of parameters
    minimizer->SetFCN(&Jackknife);		// Sets the function to minimize (jackknife)
    minimizer->SetPrintLevel(-1);		// Set verbosity to 0

    // Set the limits of the parameters to fit. It is preferable to fit the number of
    // bin as a double as it makes the optimization function continuous
    double start(10.5), step(2), pmin(2), pmax(1000);
    int ierflg;
    for (i = 0; i < _dim; i++)
	minimizer->mnparm(i, TString::Format("N_%d", (int)i), start, step, pmin, pmax, ierflg);

    // Call MIGRAD to minimize
    std::vector<double> arglist = {100};
    minimizer->mnexcm("SIM", &(arglist[0]), 1, ierflg);

    // Extract the optimized binning
    std::vector<double> values(_dim), error(_dim);
    for (i = 0; i < _dim; i++) {
      minimizer->GetParameter(i, values[i], error[i]);
      number[i] = (size_t)(values[i]);
    }
    
    range.clear();
    delete minimizer;

  } else if ( algo == "cross") {
    // Create a cross validation function J to minimize to find the optimal binning
    // Algorithm developed by M. Solis by minimizing the MISE of the leave-one-out estimator
    // Minimize the function J with TMinuit (faster)
    TMinuit* minimizer = new TMinuit(_dim);	// Specify the number of parameters
    minimizer->SetFCN(&CrossValidation);	// Sets the function to minimize (J)
    minimizer->SetPrintLevel(-1);		// Set verbosity to 0

    // Set the limits of the parameters to fit. It is preferable to fit the number of
    // bin as a double as it makes the optimization function continuous
    double start(10.5), step(2), pmin(2), pmax(1000);
    int ierflg;
    for (i = 0; i < _dim; i++)
	minimizer->mnparm(i, TString::Format("N_%d", (int)i), start, step, pmin, pmax, ierflg);

    // Call MIGRAD to minimize
    std::vector<double> arglist = {100};
    minimizer->mnexcm("SIM", &(arglist[0]), 1, ierflg);

    // Extract the optimized binning
    std::vector<double> values(_dim), error(_dim);
    for (i = 0; i < _dim; i++) {
      minimizer->GetParameter(i, values[i], error[i]);
      number[i] = (size_t)(values[i]);
    }
    
    range.clear();
    delete minimizer;

  } else if ( algo == "optbins") {
    // Create a likelihood function to minimize to find the optimal binning
    // Alogrithm developed by K.H. Knuth from bayesian statistics
    std::vector<double> likelihood;
    std::vector<size_t> number(_dim);
    size_t nmax = 100; // TODO
    int k, temp;
    for (i = 0; i < pow(nmax, _dim); i++) {

      // Extract the number of bins in each dimension
      temp = i;
      for (k = _dim-1; k > -1; k++) {
        number[k] = temp/pow(100, k);
	temp -= number[k]*pow(100, k);
      }

      // Bin out the data into a histogram
//      Histogram hist(data, number, lower, upper);

      // Calculate and store the log likelihood
//      size_t N = data->size();
//      double logM = 
    }
    
  } else {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Optimization algorithm not recognized: "+algo,
	  "OptimalBinning::Initialize"));
  }

  // If debug mode, print out the optimal amount of bins
//  if ( number.size() > 1 )
//      number[1] = number[0];
  std::string bins_str;
  for (i = 0; i < number.size(); i++)
      bins_str += std::to_string(number[i])+"  ";
  Pitch::print(Pitch::debug, "Optimal binning: "+bins_str);

  // Initialize the histogram, scale it to density
  _hist = Histogram(number, lower, upper);
  _hist.FillN(points);
  _hist.ScaleToDensity();

  // Initialize the interpolator if requested
  if ( _interp )
      _int = Interpolator(_hist.GetInterpolationGrid(), _ex);
}

double OptimalBinning::Evaluate(const std::vector<double>& v) const {

  // If the interpolation is requested, leave it to the interpolator
  if ( _interp )
      return _int(v);

  // If not, find the best bin and return its density
  return _hist.GetBinContent(_hist.GetBinID(v));
}

double OptimalBinning::Evaluate(const double& v) const {

  return Evaluate(std::vector<double>({v}));
}

double OptimalBinning::Evaluate(const double* v) const {

  std::vector<double> vv(v, v+_dim);
  return Evaluate(vv);
}

std::vector<TPolyLine*> OptimalBinning::Meshing() const {

  // Throw if there is no interpolation on this estimate
  try {
    std::vector<TPolyLine*> meshing;
    if ( _interp ) {
      throw(Exceptions::Exception(Exceptions::recoverable,
	    "No interpolation requested",
	    "OptimalBinning::Meshing"));
      return meshing;
    }

    return _int.Meshing();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::recoverable,
	  "Could not return the interpolation mesh",
	  "OptimalBinning::Meshing"));
  }
}
