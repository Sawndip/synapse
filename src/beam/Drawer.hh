#ifndef DRAWER_HH
#define DRAWER_HH

// C++ includes
#include <iostream>
#include <string>
#include <map>

// ROOT includes
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TList.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TF1.h"
#include "TGaxis.h"
#include "TLine.h"
#include "TBox.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TCanvas.h"

// Additional includes
#include "Pitch.hh"
#include "Globals.hh"
#include "ProgressBar.hh"
#include "InfoBox.hh"
#include "ScatterGraph.hh"

/** @brief Environment of beam related objects and methods */
namespace Beam {

/** @brief Structure that stores the drawing options of each data type */
struct DrawOpt {
  const char* title = "";	///< Long name of the data type
  int line_style = 1;		///< Line style (plain, dashed, dotted, etc.)
  int line_color = 1;		///< Line color (TColor code)
  int line_width = 2;		///< Line width
  int fill_style = 0;		///< Fill pattern (TAttFill code)
  int fill_color = 0;		///< Fill color (TColor code)
  int fill_alpha = 1;		///< Fill opacity (From 0=transparent to 1=opaque)
  int marker_style = 0;		///< Marker shape (TAttMarker code)
  int marker_color = 1;		///< Marker color (TColor code)
  int marker_size  = 1;		///< Marker size
  const char* hist_opt = "";	///< ROOT histogram drawing option
  const char* graph_opt = "";	///< ROOT graph drawing option
  const char* hist_leg = "";	///< ROOT histogram legend drawing option
  const char* graph_leg = "";	///< ROOT graph legend drawing option
};

/** @brief Structure that stores the drawing options for each variable */
struct VarSet {
  std::string label;		///< Label of the variable to be drawn
  std::string unit;		///< Unit in which the variable is expressed
  double llim;			///< Lower limit of the corresponding axis
  double ulim;			///< Upper limit of the corresponding axis
  size_t nbins = 100;		///< Number of bins
};

/** @brief Stores the default drawing options and methods to produce ROOT plots
 *
 *	   Defines a list of default colors, fill styles, etc. for the different data
 * 	   types. Contains a list of methods to produce preformatted plots of the
 *   	   of the phase space data.
 */
class Drawer {
 public:
  /** @brief Default constructor, intializes the default drawing options */
  Drawer();

  /** @brief Copy constructor */
  Drawer(const Drawer& dra);

  /** @brief Equality operator */
  Drawer& operator=(const Drawer& dra);

  /** @brief Destructor */
  ~Drawer();

  /** @brief Sets custom settings for the specified variable
   *
   *  @param	var	Name of the variable
   *  @param	varset	Variable settings
   */
  void SetVariableSettings(const std::string& var,
		      	   const VarSet& varset)	{ _varsets[var] = varset; }

  /** @brief Sets the limits of the specified variable
   *
   *  @param	var	Name of the variable
   *  @param	llim	Lower limit
   *  @param	ulim	Upper limit
   */
  void SetVariableLimits(const std::string& var,
		 	 double llim,
	         	 double ulim);

  /** @brief Returns an appropriate TH1F to draw one of the phase space component
   *
   *  @param	name	Name of the histogram
   *  @param	var	Variable to be drawn
   *
   *  @return		Pointer to the initialized histogram
   */
  TH1F* New1DHistogram(const std::string& name,
		       const std::string& var) const;

  /** @brief Defines an appropriate TH1F to draw two of the phase space component
   *  
   *  @param	name	Name of the histogram
   *  @param	varx	First variable to be drawn
   *  @param	vary	Second variable to be drawn
   *
   *  @return		Pointer to the initialized histogram
   */
  TH2F* New2DHistogram(const std::string& name,
		       const std::string& varx,
		       const std::string& vary) const;

  /** @brief Returns a scatter plot with an amplitude color scale
   *  
   *  @param	name	Name of the scatter plot
   *  @param	varx	First variable to be drawn
   *  @param	vary	Second variable to be drawn
   *
   *  @return		Pointer to the initialized scatter plot
   */
  ScatterGraph* NewAmplitudeScatter(const std::string& name,
		       	            const std::string& varx,
		       	            const std::string& vary) const;

  /** @brief Initializes a stack of histograms of different data types
   *
   *  @param	hists		Map of histograms that form the stack
   *
   *  @return			Pointer to the initialized stack
   */
  THStack* NewStack(std::map<std::string, TH1F*> hists) const;

  /** @brief Intializes stat boxes for a stack of histogram of different data types
   *
   *  @param	hists		Map of histograms that form the stack
   *  @param	draw_opt	Drawing options ("E" for entries, "M" for mean, "R" for RMS)
   *  @param	hoffset		Height offset due to the presence of a statbox
   *
   *  @return			Map of pointers to each initialized stack stat box
   */
  std::map<std::string, TPaveStats*> NewStackStats(std::map<std::string, TH1F*> hists,
		       			 	   const std::string draw_opt="EMR",
			 		 	   double* hoffset=NULL) const;

  /** @brief Intializes a legend for a stack of histogram of different data types
   *
   *  @param	hists		Map of histograms that form the stack
   *  @param	hoffset		Height offset due to the presence of a statbox
   *
   *  @return			Pointer to the initialized stack legend
   */
  TLegend* NewStackLegend(std::map<std::string, TH1F*> hists,
			  const double hoffset = 0) const;

  /** @brief Draws superimposed data types histograms
   *
   *  @param	hists		Map of histograms that form the stack
   *  @param	name		File name
   *
   *  @return			Pointer to the drawn stack
   */
  THStack* DrawStack(std::map<std::string, TH1F*> hists,
		     const std::string& name) const;

  /** @brief Saves superimposed data types histograms
   *
   *  @param	hists		Map of histograms that form the stack
   *  @param	name		File name
   *  @param	ratio		Draws the ratio data/simulation
   */
  void SaveStack(std::map<std::string, TH1F*> hists,
		 const std::string& name,
		 bool ratio=false) const;

  /** @brief Initializes a multigraph of graphs of different data types
   *
   *  @param	graphs		Map of graphs that form the multigraph
   *
   *  @return			Pointer to the initialized multigraph
   */
  TMultiGraph* NewMultiGraph(std::map<std::string, TGraphErrors*> graphs) const;

  /** @brief Intializes a legend for a multigraph of different data types
   *
   *  @param	graphs		Map of graphs that form the multigraph
   *  @param	corner		Corner in which to draw the legend
   *
   *  @return			Pointer to the initialized multigraph legend
   */
  TLegend* NewMultiGraphLegend(std::map<std::string, TGraphErrors*> graphs,
			       const std::string corner="tr") const;

  /** @brief Intializes an emittance axis for fractional quantities
   *
   *  @param	func		Function of the emittance
   *  @param	miny	Minimum value on the Y axis
   *  @param	maxy	Maximum value on the Y axis
   *  @param	maxx	Maximum value on the X axis
   *
   *  @return			Pointer to the initialized multigraph legend
   */
  TGaxis* NewMultiGraphEmittanceAxis(TF1* func,
				     const double& miny,
				     const double& maxy,
				     const double& maxx) const;

  /** @brief Draws superimposed data types histograms
   *
   *  @param	graphs		Map of graphs that form the multigraph
   *  @param	name		File name
   *  @param	func		Proportionality function for emittance estimate
   *
   *  @return			Pointer to the drawn multigraph
   */
  TMultiGraph* DrawMultiGraph(std::map<std::string, TGraphErrors*> graphs,
		     	      const std::string& name,
			      TF1* func = NULL);

  /** @brief Saves superimposed data types histograms
   *
   *  @param	graphs		Map of graphs that form the multigraph
   *  @param	name		File name
   *  @param	func		Proportionality function for emittance estimate
   */
  void SaveMultiGraph(std::map<std::string, TGraphErrors*> graphs,
		      const std::string& name,
		      TF1* func = NULL);

  /** @brief Produces a TH1F of the mean Y as a function of X from a TH2F
   *
   *  @param	hist		Input 2D histogram
   *  @param	fit		Flag for Gaussian peak fitting
   *
   *  @return			Pointer to the initialized projected histogram
   */
  TH1F* ProjectionMeanX(TH2F* hist,
		        bool fit = false) const;

  /** @brief Sets the info box
   *
   *  @param	obj	MAUS info box object
   */
  void SetInfoBox(InfoBox* info)			{ _info = info; }

  /** @brief Returns a pointer to the info box */
  InfoBox* GetInfoBox() const				{ return _info; }

 private:

  /** @brief Initializes the default drawing options **/
  void Initialize();

  /** @brief Sets the histogram style according to the predefined options
   *
   *  @param	obj	ROOT histogram (may be TH1, TH2 or TH3)
   *  @param	type	Data type
   */
  void SetStyle(TH1* obj,
	   	const std::string& type) const;

  /** @brief Sets the graph style according to the predefined options
   *
   *  @param	obj	ROOT graph (may be TGraph or TGraphErrors)
   *  @param	type	Data type
   */
  void SetStyle(TGraph* obj,
	   	const std::string& type) const;

  /** @brief Draws the MICE beam line elements onto an evolution graph
   *
   *  @param	miny	Minimum value on the Y axis
   *  @param	maxy	Maximum value on the Y axis
   */
  void DrawMICEModules(const double& miny,
		       const double& maxy) const;

  /** @brief Finds the corner that leasts overlaps the graph
   *
   *  @param	mg	Map of graphs to not overlap
   *
   *  @return		Returns "tl", "tr", "bl" or "br" to specify the corner
   */
  std::string FindBestCorner(std::map<std::string, TGraphErrors*> graphs) const;

  /** @brief Saves the canvases, in order, to a single GIF
   *
   *  @param	name	Name of the object to print
   *  @param	canv	Array of canvases to save
   */
  void SaveGIF(const std::string& name,
	       const std::vector<TCanvas*>& canv) const;

  std::vector<std::string>		_types;		///< Data types	
  std::vector<std::string>		_vars;		///< Phase space variables
  std::map<std::string, DrawOpt>	_options;	///< Drawing options for each data type
  std::map<std::string, VarSet>		_varsets;	///< Drawing options for each data type
  InfoBox*				_info;		///< MAUS info box to be drawn
};
} // namespace Beam

#endif
