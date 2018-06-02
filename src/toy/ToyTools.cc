#include "ToyTools.hh"

std::map<std::string, std::vector<double>> GaussianBunch(const double& m,
							 const double& p,
							 const double& n,
							 const double& eps,
							 const double& beta,
							 const double alpha) {

  std::map<std::string, std::vector<double>> samples;
  if ( !alpha ) {
    // Spread in position and momentum of the input beam
    double sig = sqrt(m*eps*beta/p);
    double sigp = (m*eps)/sig;

    // Produce the sample as an ideal gaussian beam with requested parameters
    int time = std::chrono::duration_cast<std::chrono::nanoseconds>			
			(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
    TRandom3 rdmzer(time);
    size_t i(0);
    double x, y, px, py;
    while ( i < n ) {
      x = rdmzer.Gaus(0, sig);
      y = rdmzer.Gaus(0, sig);
      px = rdmzer.Gaus(0, sigp);
      py = rdmzer.Gaus(0, sigp);
      if ( pow(px, 2) + pow(py, 2) > pow(p, 2) )
	  continue;
      samples["x"].push_back(x);
      samples["y"].push_back(y);
      samples["px"].push_back(px);
      samples["py"].push_back(py);
      samples["pz"].push_back(p);
      i++;
    }
  } else {
    // If alpha is non zero, simply sample from a 2D correlated gaussian twice
    double gamma = (1.+alpha*alpha)/beta;
    Matrix<double> cov({{beta*m*eps/p, -alpha*m*eps}, {-alpha*m*eps, gamma*m*eps*p}});
    std::vector<double> means = {0, 0};
    DGaus fgaus(means, cov);

    size_t i(0);
    std::vector<double> temp;
    double x, y, px, py;
    while ( i < n ) {
      temp = fgaus.RandomVector();
      x = temp[0];
      px = temp[1];
      temp = fgaus.RandomVector();
      y = temp[0];
      py = temp[1];
      if ( pow(px, 2) + pow(py, 2) > pow(p, 2) )
	  continue;

      samples["x"].push_back(x);
      samples["y"].push_back(y);
      samples["px"].push_back(px);
      samples["py"].push_back(py);
      samples["pz"].push_back(p);
      i++;
    }
  } 

  return samples;
}
