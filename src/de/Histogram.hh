#ifndef HISTOGRAM_HH
#define HISTOGRAM_HH

// C++ includes
#include <iostream>
#include <ctime>
#include <cmath>
#include <vector>
#include <map>

// ROOT includes
#include "TROOT.h"
#include "TCanvas.h"
#include "TObject.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

// Other includes
#include "Bin.hh"
#include "Grid.hh"

/** @brief Data stucture for histograms of arbitray large dimensions.
 *
 * 	   Produces a bucketing, provided with a set of data points in n dimensions. Returns 
 *	   the binning, their contents and densities.
 *
 *	   The underlying structure is a Grid object. A histogram of \f$N_1\times...\times N_n\f$
 *  	   bins is simply a grid of \f$(N_1+1)\times...\times(N_n+1)\f$ points that separates
 *	   the right amount of cells.
 */
class Histogram {
 public:
  /** @brief Default constructor, initilizes everything to 0 */
  Histogram();

  /** @brief Empty constructor, builds a histogram from the requested binning
   *
   *  @param	number		Number of bins in each projection
   *  @param	lower		Lower boundaries in each projection
   *  @param	upper		Upper boundaries in each projection
   */
  Histogram(std::vector<size_t> number,
	    const std::vector<double>& lower,
	    const std::vector<double>& upper);

  /** @brief Copy constructor */
  Histogram(const Histogram& hist);

  /** @brief Equality operator */
  Histogram& operator=(const Histogram& hist);

  /** @brief Destructor */
  ~Histogram();

  /** @brief Initializes the histogram from the produced grid */
  void Initialize();

  /** @brief Fills the histogram with a given point and weight w */
  void Fill(const std::vector<double>& point, const double w=1.);

  /** @brief Fills the histogram with an array of points with weights w */
  void FillN(const std::vector<std::vector<double>>& points,
	     const std::vector<double> w=std::vector<double>());

  /** @brief Returns the bin as a Bin object (knows its indices, bounds and content) */
  Bin GetBin(const std::vector<size_t>& indices) const;

  /** @brief Returns all the bins as an array of Bin objects */
  std::vector<Bin> GetBinArray() const;

  /** @brief Returns the bin barycentre in the dimension of the space */
  std::vector<double> GetBinBarycentre(const std::vector<size_t>& indices) const;

  /** @brief Returns the bin centre in a particular dimension */
  double GetBinCentre(const size_t index, const size_t i=0) const;

  /** @brief Sets the content of bin indices */
  void SetBinContent(const std::vector<size_t>& indices, const double value);

  /** @brief Returns the content of bin indices */
  const double& GetBinContent(const std::vector<size_t>& indices) const;

  /** @brief Returns an estimation of the probability density function in bin indices */
  double GetBinDensity(const std::vector<size_t>& indices) const;

  /** @brief Returns the ids of the bin that contains the requested point */
  std::vector<size_t> GetBinID(const std::vector<double>& point) const;

  /** @brief Returns the volume of a bin in this bucketting */
  double GetBinVolume() const;

  /** @brief Returns the bin width in a particular dimension */
  const double& GetBinWidth(const size_t i=0) const	{ return _grid.GetPeriod(i); }

  /** @brief Returns the amount of entries */
  const size_t GetEntries() const			{ return _points.size(); }

  /** @brief Returns the dimension of the space */
  const size_t& GetDimension()	const			{ return _grid.GetDimension(); }

  /** @brief Returns the array of lower bounds in each dimension */
  const std::vector<double>& GetLowerBoundArray() const	{ return _grid.GetLowerBoundArray(); }

  /** @brief Returns the lower bound in dimension i */
  const double& GetLowerBound(const size_t i=0) const	{ return _grid.GetLowerBound(i); }

  /** @brief Returns the array of upper bounds in each dimension */
  const std::vector<double>& GetUpperBoundArray() const	{ return _grid.GetUpperBoundArray(); }

  /** @brief Returns the upper bound in dimension i */
  const double& GetUpperBound(const size_t i=0) const	{ return _grid.GetUpperBound(i); }

  /** @brief Returns the grid that encompasses the bins */
  const Grid& GetGrid() const				{ return _grid; }

  /** @brief Returns a grid of the bin centres on which to interpolate */
  Grid GetInterpolationGrid() const;

  /** @brief Returns the amount of total amount of bins (n-dimensional buckets) */
  size_t GetNbins() const;

  /** @brief Returns the amount of bins in one projection */
  size_t GetNbins(const size_t i=0) const		{ return _grid.GetNpoints(i)-1; }

  /** @brief Scales the histogram by a certain factor */
  void Scale(const double factor);

  /** @brief Scales the histogram contents to the estimated density */
  void ScaleToDensity(const double norm=1.);

  /** @brief Returns a ROOT histogram for dimensions below 3 */
  TObject* ToROOT(const std::string name="") const;

 private:

  Grid 					_grid;		///< Grid that contains the space division
  std::vector<std::vector<double>>	_points;	///< Points that filled the histogram 
  std::vector<double>			_contents;	///< Binned data
  std::vector<double>			_outflows;	///< Under/overflow bins (3^_dim-1 bins)
};

#endif
