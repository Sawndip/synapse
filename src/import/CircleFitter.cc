#include "CircleFitter.hh"

bool FitCircle(const std::vector<double>& x, const std::vector<double>& y,
               double& x0, double& y0, double& rad) {

  auto Chi2Function = [&x, &y](const double *par) {
    // Minimisation function computing the sum of squares of residuals
    // looping over the points
    double f = 0.0;
    for (size_t i = 0; i < x.size(); i++) {
        double u = x[i] - par[0];
        double v = y[i] - par[1];
        double dr = par[2] - std::sqrt(u*u+v*v);
        f += dr*dr;
    }
    return f;
  };

  // Wrap chi2 funciton in a function object for the fit
  // 3 is the number of fit parameters (size of array par)
  ROOT::Math::Functor fcn(Chi2Function, 3);
  ROOT::Fit::Fitter fitter;
  double pStart[3] = {0, 0, 1};
  fitter.SetFCN(fcn, pStart);
  fitter.Config().ParSettings(0).SetName("x0");
  fitter.Config().ParSettings(1).SetName("y0");
  fitter.Config().ParSettings(2).SetName("R");
  
  // do the fit 
  bool ok = fitter.FitFCN();
  if ( !ok ) {
    return ok;
  }   
  const ROOT::Fit::FitResult & result = fitter.Result();
  result.Print(std::cout);
  
  // Return the circle parameters
  x0 = result.Parameter(0);
  y0 = result.Parameter(1);
  rad = result.Parameter(2);
  return true;
}
