#include "Geometry.hh"

namespace Math {

std::vector<double> PlaneParameters(const std::vector<Vector<double>>& points) {

  // Define two coplanar vectors from the three points
  Vector<double> a({points[1][0]-points[0][0], points[1][1]-points[0][1], points[1][2]-points[0][2]}),
		 b({points[2][0]-points[0][0], points[2][1]-points[0][1], points[2][2]-points[0][2]});

  // Their cross-product yields a vector normal to the plane, scale it
  Vector<double> normal = a.cross(b);
  normal /= normal.mag();

  // The last parameter can be optained by evaluating in one of the points
  std::vector<double> pars = normal.std();
  pars.push_back(-pars[0]*points[0][0]-pars[1]*points[0][1]-pars[2]*points[0][2]);

  return pars;
}

std::vector<double> PlaneParameters(const Vector<double>& v0,
				    const Vector<double>& v1,
				    const Vector<double>& v2) {

  return PlaneParameters({v0, v1, v2});
}

double PlaneDistance(const std::vector<double>& pars, const std::vector<double>& point) {

  // Increment the numerator, the denominator is scaled to 1 by default (unit-normal)
  Assert::IsEqual("Math::PlaneDistance", "Plane and point dimensions", pars.size()-1, point.size());
  double dist(0.);
  size_t i;
  for (i = 0; i < point.size(); i++)
      dist += pars[i]*point[i];
  dist += pars.back();

  return dist;
}

bool LeftHanded(const Vector<double>& a,
		const Vector<double>& b) {

  // Return true if a is to the left of b in the 2D plane
  if (a[0] >= 0 && b[0] < 0)
      return true;
  if (a[0] < 0 && b[0] >= 0)
      return false;
  if (a[0] == 0 && b[0] == 0) {
    if (a[1] >= 0 || b[1] >= 0)
        return a[1] > b[1];

    return  b[1] > a[1];
  }

  return a[0]*b[1]-a[1]*b[0] < 0;
}

bool RightHanded(const Vector<double>& a,
		 const Vector<double>& b) {

  return !LeftHanded(a, b);
}

void SortCW(std::vector<Vector<double>>& points) {

  // Clear the easy case (all the points are 2D)
  Assert::IsNotEmpty("Math::SortCW", "Vector of input points", points);
  if ( points.size() < 3 ) {
    // No ordering to do
    return;
  }

  // For larger samples, one needs to place oneself in the CMF, first compute the CM
  Vector<double> CM(points[0].size(), 0.);
  for (const Vector<double>& point : points)
      CM += point;
  CM /= points.size();

  // Place all of the points in the CMF
  std::vector<std::pair<size_t, Vector<double>>> cmfpoints(points.size());
  size_t i;
  for (i = 0; i < cmfpoints.size(); i++)
      cmfpoints[i] = std::pair<size_t, Vector<double>>(i, (points[i]-CM));

  // Sort the points clockwise
  std::sort(cmfpoints.begin(), cmfpoints.end(),
    [] (const std::pair<size_t, Vector<double>>& a, const std::pair<size_t, Vector<double>>& b) {
	return LeftHanded(a.second, b.second);
    });

  // Set the sorted points
  std::vector<Vector<double>> temp(points.size());
  for (i = 0; i < points.size(); i++)
      temp[i] = points[cmfpoints[i].first];
  points = temp;
}

void SortACW(std::vector<Vector<double>>& points) {

  // Sort them clockwise then inverse the order
  Assert::IsNotEmpty("Math::SortACW", "Vector of input points", points);
  SortCW(points);
  std::reverse(points.begin(), points.end());
}

double BoundingBox(const std::vector<std::vector<double>>& points,
		   Vector<double>& lower,
		   Vector<double>& upper,
		   double buff) {

  // First find the minimum and the maximum of the range in each of the dimension
  Assert::IsNotEmpty("Math::BoundingBox", "Vector of input points", points);
  Assert::IsNonZero("Math::BoundingBox", "Dimension", points[0].size());
  size_t dim = points[0].size();
  lower.resize(dim);
  upper.resize(dim);
  size_t i, j;
  for (i = 0; i < dim; i++) {
    lower[i] = points[0][i];
    upper[i] = points[0][i];
    for (j = 1; j < points.size(); j++) {
      if ( points[j][i] < lower[i] )
	  lower[i] = points[j][i];
      if ( points[j][i] > upper[i] )
	  upper[i] = points[j][i];
    }
  }

  // If the buffer is non-zero, add a buffer zone between the extrema and the box
  // In any case, calculate the volume of the bounding box
  double vol(1.), temprange;
  for (i = 0; i < dim; i++) {
    if ( buff ) {
      temprange = upper[i]-lower[i];
      lower[i] -= buff*temprange;
      upper[i] += buff*temprange;
    }
    vol *= upper[i]-lower[i];
  }

  return vol;
}

std::vector<std::vector<double>> SimplexVertices(const size_t n) {

  // The cross product between two vertices is always -1/n
  Assert::IsNonZero("Math::SimplexVertices", "Dimension", n);
  double cross = -1./n;

  // Initialize the points
  std::vector<std::vector<double>> vertices(n+1);
  size_t i, j;
  for (i = 0; i < n+1; i++)
      vertices[i].resize(n);

  // Loop over the dimension and set the coordinates
  double tmp(0), ncoord, rcoord;
  for (i = 0; i < n; i++) {
    ncoord = sqrt(1-tmp);
    vertices[i][i] = ncoord;
    rcoord = (cross - tmp)/ncoord;

    for(j = i+1; j < n+1; j++)
        vertices[j][i] = rcoord;

    tmp += rcoord*rcoord;
  }

  return vertices;
}

bool IsInsideSimplex(const std::vector<double>& point,
		     const std::vector<std::vector<double>>& simplex) {

  // Check that point and simplex dimensionally match
  // TODO
  size_t d = point.size();

  // Compute the oriented hypervolumes of parallelotopes associated with the simplex facets.
  // If any of them has an incompatible sign, the point is outside the simplex
  // A_j,k = (x_1-x_k, ..., x_j-1,k, x_j+1,k, ..., x_d+1,k)^T
  // sign(det(A_j,j)) = -sign(det(A_j+1,j+1))
  Matrix<double> A11(d, d);
  size_t m, n, j;
  for (m = 0; m < d; m++)
    for (n = 0; n < d; n++)
	A11[m][n] = simplex[m+1][n]-simplex[0][n];

  double det = A11.Determinant();
  Matrix<double> Aj0(d, d);
  for (j = 0; j < d+1; j++) {
    for (m = 0; m < d; m++)
      for (n = 0; n < d; n++)
	  Aj0[m][n] = m < j ? simplex[m][n] - point[n] : simplex[m+1][n] - point[n];

    if ( (Aj0.Determinant() > 0) != (pow(-1, j)*det > 0) )
	return false;
  }
  
  return true;
}

double LpNormp(const std::vector<double>& v, const double& p) {

  Assert::IsNotEmpty("Math::LpNormp", "Vector", v);
  double norm(0.);
  size_t i;
  for (i = 0; i < v.size(); i++)
      norm += pow(v[i], p);
  return norm;
}

double LpNorm(const std::vector<double>& v, const double& p) {

  Assert::IsNotEmpty("Math::LpNorm", "Vector", v);
  return pow(LpNormp(v, p), 1./p);
}

double UnitBallVolume(const size_t& d, const double p) {

  Assert::IsNonZero("Math::UnitBallVolume", "Dimension", d);
  if ( d == 1 ) 
      return 2.;

  if ( p == 2 )
      return pow(M_PI, (double)d/2)/tgamma((double)d/2+1.);

  return pow(2, d)*pow(tgamma(1./p+1), d)/tgamma((double)d/p+1);
}

double HullVolumeUniformFactor(const double& d, const double& p) {

  if ( p != 2 ) {
    Pitch::print(Pitch::warning, "Non-Euclidean metrics currently not supported");
    return 1.;
  }

  return .5*((d+1.)/(d+3.))
	*(std::tgamma((d+3.)/(d+1.)+d)/Math::Factorial(d-1))
	*pow(2*sqrt(M_PI)*std::tgamma(1.5+d/2)/std::tgamma(1.+d/2), 2./(d+1.));
}

double HullVolumeUniformRelativeBias(const double& d, const double& p, const size_t N) {

  // TODO TODO TODO, needs to dampen high bias
  return 1.-HullVolumeUniformFactor(d, p)*pow(N, -2/(d+1));
}

double HullVolumeTGausRelativeBias(const double& d, const double& alpha, const size_t N) {

  // TODO TODO TODO, function of alpha
  return HullVolumeUniformRelativeBias(d, 2, alpha*N);
}

double HullVolumeUniformRelativeRMS(const double& d, const double& p, const size_t N) {

  return sqrt(2)*pow(N, -(d+3.)/(2*(d+1.)));
}

double HullVolumeTGausRelativeRMS(const double& d, const double& alpha, const size_t N) {

  // TODO TODO TODO, breaks down for high alpha
  return sqrt(exp(TMath::ChisquareQuantile(alpha, d)/2)*(1.-alpha)/(alpha*N));
}
} // namespace Geometry
