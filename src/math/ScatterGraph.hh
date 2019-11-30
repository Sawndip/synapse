#ifndef SCATTERGRAPH_HH
#define SCATTERGRAPH_HH

// C++ includes
#include <iostream>
#include <vector>
#include <algorithm>

// ROOT includes
#include "TStyle.h"
#include "TGraph2D.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TCanvas.h"
#include "TH2F.h"

/** @brief Enumerates the possible ordering for the points 
  *
  *        List of the possible order in which the points are drawn:
  *	     - ascending: ascending order of mapping
  *	     - descending: descending order of mapping
  */
enum Order {ascending, descending};

/** @brief Defines functions that allow to plot a ROOT scatter graph with a color scale.
 *
 *         It is a class derives from the TGraph2D class drawn with the PCOLZ drawing option,
 *	   viewed from top, with updated TGaxis that fit the canvas in an elegant way.
 */

class ScatterGraph {
 public:

  /** @brief Default constructor - initialises to 0/NULL */
  ScatterGraph();

  /** @brief Copy constructor - any pointers are deep copied */
  ScatterGraph(const ScatterGraph& scat);

  /** @brief Equality operator - any pointers are deep copied */
  ScatterGraph& operator=(const ScatterGraph& scat);

  /** @brief Constructor with number of points to be hosted */
  ScatterGraph(size_t n);

  /** @brief Constructor with data */
  ScatterGraph(size_t n, double* x, double* y, double* z);

  /** @brief Full constructor with name, title and data */
  ScatterGraph(const char *name, const char *title, size_t n, double* x, double* y, double* z);

  /** @brief Initializes the underlying TGraph2D from another TGraph2D */
  ScatterGraph(const TGraph2D& graph);

  /** @brief Destructor - any member pointers are deleted */
  virtual ~ScatterGraph();

  /** @brief Sets one of the points to a certain value */
  void AddPoint(double x, double y, double w=1);

  /** @brief Sets one of the points to a certain value */
  void SetPoint(size_t i, double x, double y, double w=1)	{ _graph->SetPoint(i, x, y, w); }

  /** @brief Sets one of the points to a certain value */
  void SetPoints(size_t n, double* x, double* y, double* w);

  /** @brief Sets the order in which the points are drawn */
  void SetOrder(const Order& order);

  /** @brief Sets the maximum amount of points to be drawn */
  void SetMaxSize(const size_t& size);

  /** @brief Draws the scatter plot */
  void Draw();

  /** @brief Saves the scatter plot */
  void Save();

  /** @brief Sets the name of the underlying TGraph2D */
  void SetName(const char* name) 				{ _graph->SetName(name); }

  /** @brief Returns the name of the underlying TGraph2D */
  virtual const char* GetName() const 				{ return _graph->GetName(); }

  /** @brief Sets the name of the underlying TGraph2D */
  void SetTitle(const char* title) 				{ _graph->SetTitle(title); }

  /** @brief Returns the title of the underlying TGraph2D */
  virtual const char* GetTitle() const 				{ return _graph->GetTitle(); }

  /** @brief Returns the number of entries in the underlying TGraph2D */
  size_t GetN() const	    					{ return _graph->GetN(); }

  /** @brief Returns a pointer to the underlying array of x values */
  double* GetX() const	    					{ return _graph->GetX(); }

  /** @brief Returns a pointer to the underlying array of y values */
  double* GetY() const	    					{ return _graph->GetY(); }

  /** @brief Returns a pointer to the underlying array of z values */
  double* GetZ() const	    					{ return _graph->GetZ(); }

  /** @brief Returns the the underlying TGraph2D */
  TGraph2D* GetGraph2D()	    				{ return _graph; }

  /** @brief Returns a pointer to the x axis, in order to set its attributes */
  TGaxis* GetXaxis()						{ return _xaxis; }

  /** @brief Returns a pointer to the y axis, in order to set its attributes */
  TGaxis* GetYaxis()						{ return _yaxis; }

  /** @brief Returns a pointer to the y axis, in order to set its attributes */
  TAxis* GetZaxis()						{ return _graph->GetZaxis(); }

  /** @brief Sets the minimum of the underlying TGraph2D */
  void SetMinimum(const double min)				{ _graph->SetMinimum(min); }

  /** @brief Sets the minimum of the underlying TGraph2D */
  void SetMaximum(const double max)				{ _graph->SetMaximum(max); }

  /** @brief Sets the limits of the x and y axes */
  void SetRange(double minx, double miny, double maxx, double maxy);

  /** @brief Sets the marker style of the underlying TGraph2D */
  void SetMarkerStyle(int style)				{ _graph->SetMarkerStyle(style); }

  /** @brief Sets the marker size of the underlying TGraph2D */
  void SetMarkerSize(double size)				{ _graph->SetMarkerSize(size); }

  /** @brief Sets the title of the x axis */
  void SetXaxisTitle(const char *title="")			{ _xaxis->SetTitle(title); }

  /** @brief Sets the title of the y axis */
  void SetYaxisTitle(const char *title="")			{ _yaxis->SetTitle(title); }

  /** @brief Save the histogram to the root file */
  void Write()							{ _graph->Write(); }

 private:

  TGraph2D*	_graph;		///< Underlying TGraph2D object
  TGaxis*  	_xaxis; 	///< Underlying TGaxis object defined for the x axis
  TGaxis*  	_yaxis;		///< Underlying TGaxis object defined for the y axis
};

#endif
