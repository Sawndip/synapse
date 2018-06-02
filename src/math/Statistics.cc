#include "Statistics.hh"

namespace Math {

double Mean(const std::vector<double>& vx) {

  Assert::IsNotEmpty("Math::Mean", "Sample", vx);

  double mean(0.);
  for (const double &x : vx)
      mean += x;
  mean /= vx.size();
  return mean;
}

double MeanError(const std::vector<double>& vx,
		 const std::vector<double>& vex) {

  Assert::IsNotEmpty("Math::MeanError", "Sample", vx);
  if ( !vex.size() )
      return MeanSError(vx);

  Assert::SameSize("Math::MeanError", "Sample and errors", vx, vex);
  return sqrt(pow(MeanMError(vex), 2) + pow(MeanSError(vx), 2));
}

double MeanMError(const std::vector<double>& vex) {

  Assert::IsNotEmpty("Math::MeanMError", "Error vector", vex);

  double mean_error(0.);
  for (const double &ex : vex)
      mean_error += ex*ex;
  mean_error = sqrt(mean_error);
  return mean_error/vex.size();
}

double MeanSError(const std::vector<double>& vx) {

  Assert::IsNotEmpty("Math::MeanSError", "Sample", vx);
  return sqrt(Variance(vx)/vx.size());
}

double RobustMean(const std::vector<double>& vx) {

  TRobustEstimator robust;
  std::vector<double> vec(vx);
  double mean, sigma;
  robust.EvaluateUni(vec.size(), &vec[0], mean, sigma);
  return mean;
}

double TMean(const std::vector<double>& vx,
	     const double& frac) {

  Assert::IsNotEmpty("Math::TMean", "Sample", vx);
  return Mean(Trimmed(vx, frac));
}

double TMeanSError(const std::vector<double>& vx,
	     	   const double& frac) {

  Assert::IsNotEmpty("Math::TMeanSError", "Sample", vx);
  return WRMS(vx, frac)/(1.-2*frac)/sqrt(vx.size());
}

double WMean(const std::vector<double>& vx,
	     const double& frac) {

  Assert::IsNotEmpty("Math::WMean", "Sample", vx);
  return Mean(Winsorized(vx, frac));
}

void DecrementMean(const size_t& n,
		   double& mean,
		   const double& x) {

  Assert::IsNonZero("Math::DecrementMean", "Current number of elements", n);
  if ( n < 2 ) {
    mean = 0.;
  } else {
    mean = (n*mean-x)/(n-1);
  }
}

void IncrementMean(const size_t& n,
		   double& mean,
		   const double& x) {

  mean = (n*mean+x)/(n+1);
}

double Variance(const std::vector<double>& vx) {

  Assert::HasAtLeast("Math::Variance", "Sample", vx, 2);

  double xmean(Mean(vx));
  double var(0.);
  size_t i;
  for (i = 0; i < vx.size(); i++)
      var += (vx[i]-xmean)*(vx[i]-xmean);
  var /= (vx.size()-1);
  return var;
}

double VarianceError(const std::vector<double>& vx,
		     const std::vector<double>& vex) {

  Assert::HasAtLeast("Math::VarianceError", "Sample", vx, 2);
  if ( !vex.size() )
      return VarianceSError(vx);

  Assert::SameSize("Math::VarianceError", "Sample and errors", vx, vex);
  return sqrt(pow(VarianceMError(vx, vex), 2) + pow(VarianceSError(vx), 2));
}

double VarianceMError(const std::vector<double>& vx,
		      const std::vector<double>& vex) {

  Assert::HasAtLeast("Math::VarianceMError", "Sample", vx, 2);
  Assert::SameSize("Math::VarianceMError", "Sample and errors", vx, vex);

  double xmean(Mean(vx));
  double var_error(0.);
  size_t i;
  for (i = 0; i < vx.size(); i++)
      var_error += pow((vx[i]-xmean)*vex[i], 2);
  var_error = sqrt(var_error);
  var_error *= 2./(vx.size()-1);
  return var_error;
}

double VarianceSError(const std::vector<double>& vx) {

  Assert::HasAtLeast("Math::VarianceSError", "Sample", vx, 2);
  return Variance(vx)*sqrt(2./(vx.size()-1));
}

double RMS(const std::vector<double>& vx) {

  Assert::HasAtLeast("Math::RMS", "Sample", vx, 2);
  return sqrt(Variance(vx));
}

double RMSError(const std::vector<double>& vx,
		const std::vector<double>& vex) {

  Assert::HasAtLeast("Math::RMSError", "Sample", vx, 2);
  if ( !vex.size() )
      return RMSSError(vx);

  Assert::SameSize("Math::RMSError", "Sample and errors", vx, vex);
  return sqrt(pow(RMSMError(vx, vex), 2) + pow(RMSSError(vx), 2));
}

double RMSMError(const std::vector<double>& vx,
		 const std::vector<double>& vex) {

  Assert::HasAtLeast("Math::RMSMError", "Sample", vx, 2);
  return VarianceMError(vx, vex)/(2*RMS(vx));
}

double RMSSError(const std::vector<double>& vx) {

  Assert::HasAtLeast("Math::RMSSError", "Sample", vx, 2);
  return RMS(vx)*sqrt(1./(2*(vx.size()-1)));
}

double WRMS(const std::vector<double>& vx,
	    const double& frac) {

  Assert::HasAtLeast("Math::WRMS", "Sample", vx, 2);
  return RMS(Winsorized(vx, frac));
}

double Covariance(const std::vector<double>& vx,
		  const std::vector<double>& vy) {

  Assert::HasAtLeast("Math::Covariance", "Sample", vx, 2);
  Assert::SameSize("Math::Covariance", "Samples", vx, vy);

  double xmean(Mean(vx)), ymean(Mean(vy));
  double cov(0.);
  size_t i;
  for (i = 0; i < vx.size(); i++)
      cov += (vx[i]-xmean)*(vy[i]-ymean);
  cov /= (vx.size()-1);
  return cov;
}

double CovarianceError(const std::vector<double>& vx,
		       const std::vector<double>& vy,
		       const std::vector<double>& vex,
		       const std::vector<double>& vey) {

  Assert::HasAtLeast("Math::CovarianceError", "Sample", vx, 2);
  Assert::SameSize("Math::CovarianceError", "Samples", vx, vy);
  if ( !vex.size() || !vey.size() )
      return CovarianceSError(vx, vy);

  Assert::SameSize("Math::CovarianceError", "First sample and errors", vx, vex);
  Assert::SameSize("Math::CovarianceError", "Second sample and errors", vy, vey);
  return sqrt(pow(CovarianceMError(vx, vy, vex, vey), 2) + pow(CovarianceSError(vx, vy), 2));
}

double CovarianceMError(const std::vector<double>& vx,
			const std::vector<double>& vy,
		        const std::vector<double>& vex,
			const std::vector<double>& vey) {

  Assert::HasAtLeast("Math::CovarianceMError", "Sample", vx, 2);
  Assert::SameSize("Math::CovarianceMError", "Samples", vx, vy);
  Assert::SameSize("Math::CovarianceMError", "First sample and errors", vx, vex);
  Assert::SameSize("Math::CovarianceMError", "Second sample and errors", vy, vey);

  double xmean(Mean(vx)), ymean(Mean(vy));
  double cov_error(0.);
  size_t i;
  for (i = 0; i < vx.size(); i++) {
    cov_error += pow((vy[i]-ymean)*vex[i], 2);
    cov_error += pow((vx[i]-xmean)*vey[i], 2);
  }
  cov_error = sqrt(cov_error);
  cov_error /= (vx.size()-1);
  return cov_error;
}

double CovarianceSError(const std::vector<double>& vx,
			const std::vector<double>& vy) {

  Assert::HasAtLeast("Math::CovarianceSError", "Sample", vx, 2);
  Assert::SameSize("Math::CovarianceSError", "Samples", vx, vy);
  return Covariance(vx, vy)*sqrt(2./(vx.size()-1));
}

void DecrementCovariance(const size_t& n,
		     	 double& cov,
		     	 const double& xmean,
		     	 const double& ymean,
		     	 const double& x,
		     	 const double& y) {

  Assert::IsNonZero("Math::DecrementCovariance", "Current number of elements", n);
  if ( n < 3 ) {
    cov = 0.;
  } else {
    cov = (n-1)*cov/(n-2)-n*(x-xmean)*(y-ymean)/(n-1)/(n-2);
  }
}

void IncrementCovariance(const size_t& n,
		     	 double& cov,
		     	 const double& xmean,
		     	 const double& ymean,
		     	 const double& x,
		     	 const double& y) {

  if ( !n ) {
    cov = 0.;
  } else {
    cov = (n-1)*cov/n+(x-xmean)*(y-ymean)/(n+1);
  }
}

Matrix<double> CovarianceMatrix(const Matrix<double>& mat) {

  // If the matrix is empty or has a single element per sample, throw
  Assert::IsNotEmpty("Math::CovarianceMatrix", "Matrix", mat.std());
  Assert::HasAtLeast("Math::CovarianceMatrix", "Sample", mat[0], 2);

  // Fill the covariance matrix
  Matrix<double> S(mat.Nrows(), mat.Nrows());
  size_t i, j;
  for (i = 0; i < mat.Nrows(); i++)
    for (j = 0; j < mat.Nrows(); j++)
      S[i][j] = (j < i) ? S[j][i] : Covariance(mat[i], mat[j]);
      
  return S;
}

Matrix<double> RobustCovarianceMatrix(const Matrix<double>& mat,
				      const double hh) {

  // If the matrix is empty or has a single element per sample, throw
  // Copy the matrix as ROOT is not const safe...
  Assert::IsNotEmpty("Math::CovarianceMatrix", "Matrix", mat.std());
  Assert::HasAtLeast("Math::CovarianceMatrix", "Sample", mat[0], 2);
  Matrix<double> sample(mat);

  // Feed the data to the ROOT implementation of MCD
  TRobustEstimator robust(sample.Ncols(), sample.Nrows(), hh*sample.Ncols());
  size_t i, j;
  for (i = 0; i < sample.Nrows(); i++)
      robust.AddColumn(&sample[i][0]);
  robust.Evaluate();

  // Extract the covariance matrix
  TMatrixDSym covmat = *robust.GetCovariance();
  Matrix<double> S(sample.Nrows(), sample.Nrows());
  for (i = 0; i < sample.Nrows(); i++)
    for (j = 0; j < sample.Nrows(); j++)
      S[i][j] = (j < i) ? S[j][i] : covmat[i][j];

  return S;
}

void DecrementCovarianceMatrix(const size_t& n,
		     	       Matrix<double>& covmat,
		     	       const std::vector<double>& means,
		     	       const std::vector<double>& v) {

  Assert::IsNonZero("Math::DecrementCovarianceMatrix", "Current number of elements", n);

  // Update each of the covariance elements on the upper triangle, copy to the lower one
  size_t i, j;
  for (i = 0; i < covmat.Nrows(); i++)
    for (j = 0; j < covmat.Nrows(); j++)
      if ( j < i ) {
	covmat[i][j] = covmat[j][i];
      } else {
	DecrementCovariance(n, covmat[i][j], means[i], means[j], v[i], v[j]);
      }
}

void IncrementCovarianceMatrix(const size_t& n,
		     	       Matrix<double>& covmat,
		     	       const std::vector<double>& means,
		     	       const std::vector<double>& v) {

  // Update each of the covariance elements on the upper triangle, copy to the lower one
  size_t i, j;
  for (i = 0; i < covmat.Nrows(); i++)
    for (j = 0; j < covmat.Nrows(); j++)
      if ( j < i ) {
	covmat[i][j] = covmat[j][i];
      } else {
	IncrementCovariance(n, covmat[i][j], means[i], means[j], v[i], v[j]);
      }
}

double DeterminantError(const Matrix<double>& mat,
			const Matrix<double>& mat_error) {

  // If the matrix is empty or has a single element per sample, throw
  Assert::IsNotEmpty("Math::DeterminantError", "Matrix", mat.std());
  Assert::HasAtLeast("Math::DeterminantError", "Sample", mat[0], 2);
  if ( !mat_error.Nrows() )
      return DeterminantSError(mat);

  Assert::SameSize("Math::DeterminantError", "Matrix and error matrix", mat.std(), mat_error.std());
  return sqrt(pow(DeterminantMError(mat, mat_error), 2)+pow(DeterminantSError(mat), 2));
}

double DeterminantMError(const Matrix<double>& mat,
			 const Matrix<double>& mat_error) {

  // If the matrix is empty or has a single element per sample, throw
  Assert::IsNotEmpty("Math::DeterminantMError", "Matrix", mat.std());
  Assert::HasAtLeast("Math::DeterminantMError", "Sample", mat[0], 2);
  Assert::SameSize("Math::DeterminantMError", "Matrix and error matrix", mat.std(), mat_error.std());

  // Fetch the means
  std::vector<double> means;
  size_t i, j, k;
  for (i = 0; i < mat.Nrows(); i++)
      means.push_back(Mean(mat[i]));

  // Get the cofactor matrix
  Matrix<double> C = CovarianceMatrix(mat).CofactorMatrix();

  // Compute the variance
  double var(0.);
  Matrix<double> E(mat.Nrows(), mat.Nrows(), 0.);
  Matrix<double> CEC(mat.Nrows(), mat.Nrows());
  for (i = 0; i < mat.Ncols(); i++) {

    // Produce an nxn diagonal matrix with the errors on the variables
    for (j = 0; j < mat.Nrows(); j++)
        E[j][j] = pow(mat_error[j][i], 2);

    // Produce C^TEC that you need the elements of (C=C^T, symmetric matrix)
    CEC = C*E*C;

    // Loop over all the combination of keys and add
    for (j = 0; j < mat.Nrows(); j++)
      for (k = 0; k < mat.Nrows(); k++)
	var += CEC[j][k]*(mat[j][i]-means[j])*(mat[k][i]-means[k]);
  }

  // Add the factor
  return 2*sqrt(var)/(mat.Ncols()-1);
}

double DeterminantSError(const Matrix<double>& mat) {

  // If the matrix is empty or has a single element per sample, throw
  Assert::IsNotEmpty("Math::DeterminantSError", "Matrix", mat.std());
  Assert::HasAtLeast("Math::DeterminantSError", "Sample", mat[0], 2);

  // Return the determinant statistical error
  return CovarianceMatrix(mat).Determinant()*sqrt(2.*mat.Nrows()/(mat.Ncols()-1));
}

double Median(const std::vector<double>& vx) {

  Assert::IsNotEmpty("Math::Median", "Sample", vx);

  size_t n = vx.size();
  std::vector<double> values = vx; // Preserve the original vector
  std::sort(values.begin(), values.end());
  if ( values.size() % 2 )
      return values[n/2];

  return .5*(values[n/2-1]+values[n/2]);
}

double MedianSError(const size_t n, const double& density) {

  Assert::IsNonZero("Math::MedianSError", "Number of points", n);
  return sqrt(1./(4.*(n+2)))/density;
}

double LWidth(const std::vector<double>& vx, const double& p) {

  Assert::HasAtLeast("Math::LWidth", "Sample", vx, 2);

  double xmed(Median(vx));
  double width(0.);
  size_t i;
  for (i = 0; i < vx.size(); i++)
      width += pow(fabs(vx[i]-xmed), p);
  width /= (vx.size()-1);

  return pow(width, 1./p);
}

double Quantile(const std::vector<double>& vx, const double& alpha) {

  Assert::IsNotEmpty("Math::Quantile", "Sample", vx);
  Assert::IsProbability("Math::Quantile", "alpha", alpha, false, false);

  std::vector<double> values = vx; // Preserve the original vector
  std::sort(values.begin(), values.end());
  double pos = alpha*(values.size()-1.);
  size_t rank = (int)pos;
  if ( pos == rank )
      return values[rank];

  return values[rank]+(pos-rank)*(values[rank+1]-values[rank]);
}

double QuantileSError(const size_t n, const double& alpha, const double& density) {

  Assert::IsNonZero("Math::QuantileSError", "Number of points", n);
  Assert::IsProbability("Math::QuantileSError", "alpha", alpha, false, false);

  return sqrt(alpha*(1.-alpha)/(n+2))/density;
}

double Mode(const std::vector<double>& vx) {

  // Check that the input sample is not empty
  Assert::IsNotEmpty("Math::Mode", "Sample", vx);

  // Define the L0 (Lebesgue space) penalty function necessary to find the mode
  auto penalty_L0 =
  [] (double* x, double* par) {
    size_t n = par[0];
    size_t i;
    double sum(0);
    for (i = 0; i < n; i++)
        sum += fabs(x[0]-par[i+1])/(1+fabs(x[0]-par[i+1]));

    return sum;
  };

  // Define a function TF1 as the penalty function
  size_t N = vx.size();
  double min(Min(vx)), max(Max(vx));
  TF1* fdist0 = new TF1("fdist0", penalty_L0, min, max, N+1);
  double params[N+1];
  params[0] = N;
  size_t i;
  for (i = 0; i < N; i++)
      params[i+1] = vx[i];
  fdist0->SetParameters(params);

  // Minimize the penalty function
  return fdist0->GetMinimumX(min, max);
}

double Min(const std::vector<double>& vx) {

  // Check that the input sample is not empty
  Assert::IsNotEmpty("Math::Min", "Sample", vx);

  // std::min_element returns the iterator, the distance returns the rank
  std::vector<double>::const_iterator it = std::min_element(vx.begin(), vx.end());
  size_t rank = std::distance(vx.begin(), it);
  return vx[rank];
}

double Max(const std::vector<double>& vx) {

  // Check that the input sample is not empty
  Assert::IsNotEmpty("Math::Max", "Sample", vx);

  // std::max_element returns the iterator, the distance returns the rank
  std::vector<double>::const_iterator it = std::max_element(vx.begin(), vx.end());
  size_t rank = std::distance(vx.begin(), it);
  return vx[rank];
}

double Mahalanobis(const Matrix<double>& invcov, const Matrix<double>& vec) {

  return sqrt(MahalanobisSquared(invcov, vec));
}

double MahalanobisSquared(const Matrix<double>& invcov, const Matrix<double>& vec) {

  // Assert that the input is sensible
  Assert::IsNotEmpty("Math::MahalanobisSquared", "Inverse covariance matrix", invcov.std());
  Assert::IsSquare("Math::MahalanobisSquared", "Inverse covariance matrix", invcov.std());
  Assert::IsEqual("Math::MahalanobisSquared",
	"Inverse covariance rank and vector size", invcov.Nrows(), vec.Ncols());
  if ( vec.Ncols() != 1 )
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	               		  "Second input is not a vector",
				  "Math::MahalanobisSquared"));

  // Compute the squared Mahalanobis distance
  Matrix<double> vecT = vec.Transpose();
  return (vecT*invcov*vec)[0][0];
}

double GeneralisedDistance(const double& rhoR,
			   const double& rho0) {

  return sqrt(GeneralisedDistanceSquared(rhoR, rho0));
}

double GeneralisedDistanceSquared(const double& rhoR,
				  const double& rho0) {

  Assert::IsNonZero("Math::GeneralisedDistanceSquared", "Maximum of density", rho0);
  return -2*log(rhoR/rho0);
}

size_t Factorial(const size_t n) {

  if ( n < 2 )
      return 1;

  size_t i, fact(1);
  for (i = 2; i < n+1; i++)
      fact *= i;
  return fact;
}

} // namespace Math
