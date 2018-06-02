#include "ScatterGraph.hh"

ScatterGraph::ScatterGraph()
  : _graph(), _xaxis(), _yaxis() {
  _graph = new TGraph2D();
  _xaxis = new TGaxis();
  _yaxis = new TGaxis();
}

ScatterGraph::ScatterGraph(size_t n)
  : _graph(), _xaxis(), _yaxis() {
  _graph = new TGraph2D(n);
  _xaxis = new TGaxis();
  _yaxis = new TGaxis();
}

ScatterGraph::ScatterGraph(size_t n, double* x, double* y, double* z)
  : _graph(), _xaxis(), _yaxis() {
  _graph = new TGraph2D(n, x, y, z);
  _xaxis = new TGaxis();
  _yaxis = new TGaxis();
}

ScatterGraph::ScatterGraph(const char *name, const char *title, size_t n, double* x, double* y, double* z)
  : _graph(), _xaxis(), _yaxis() {
  _graph = new TGraph2D(name, title, n, x, y, z);
  _xaxis = new TGaxis();
  _yaxis = new TGaxis();
}

ScatterGraph::ScatterGraph(const TGraph2D& graph)
  : _graph(), _xaxis(), _yaxis() {

  _graph = new TGraph2D(graph);
  _xaxis = new TGaxis();
  _yaxis = new TGaxis();
}

ScatterGraph::ScatterGraph(const ScatterGraph& scat) {

  *this = scat;
}

ScatterGraph& ScatterGraph::operator=(const ScatterGraph& scat) {
  if (this == &scat)
    return *this;

  if ( _graph )
      delete _graph;
  _graph = new TGraph2D(*scat._graph);
  if ( _xaxis )
      delete _xaxis;
  _xaxis = new TGaxis();
  if ( _yaxis )
      delete _yaxis;
  _yaxis = new TGaxis();

  return *this;
}

ScatterGraph::~ScatterGraph() {

  if ( _graph )
      delete _graph;
  if ( _xaxis )
      delete _xaxis;
  if ( _yaxis )
      delete _yaxis;
}

void ScatterGraph::AddPoint(double x, double y, double w) {

  size_t N = GetN();
  _graph->SetPoint(N, x, y, w);
}

void ScatterGraph::SetPoints(size_t n, double* x, double* y, double* w) {

  size_t N = GetN();
  size_t i;
  for (i = 0; i < n; i++)
      _graph->SetPoint(N+i, x[i], y[i], w[i]);
}

void ScatterGraph::SetRange(double minx, double miny, double maxx, double maxy) {

  TH2F* hdummy = new TH2F("dummy", "", 1, minx, maxx, 1, miny, maxy);
  _graph->SetHistogram(hdummy);
}

void ScatterGraph::SetOrder(const Order& order) {

  // Order the points according to the requested order
  std::vector<std::vector<double>> sorted(GetN());
  for (size_t i = 0; i < sorted.size(); i++)
        sorted[i] = {_graph->GetX()[i], _graph->GetY()[i], _graph->GetZ()[i]};

  if ( order == descending ) {
    std::sort(sorted.begin(), sorted.end(),
	      [] (const std::vector<double>& a,  const std::vector<double>& b)
	      { return a[2] > b[2]; });
  } else {
    std::sort(sorted.begin(), sorted.end(),
	      [] (const std::vector<double>& a,  const std::vector<double>& b)
	      { return a[2] < b[2]; });
  }

  std::vector<double> x(GetN()), y(GetN()), z(GetN());
  for (size_t i = 0; i < x.size(); i++) {
    x[i] = sorted[i][0];
    y[i] = sorted[i][1];
    z[i] = sorted[i][2];
  }

  SetPoints(x.size(), &x[0], &y[0], &z[0]);
}

void ScatterGraph::SetMaxSize(const size_t& size) {

  // If the cloud is larger than size, truncate the sample
  if ( (size_t)_graph->GetN() > size )
      _graph->Set(size);
}

void ScatterGraph::Draw() {

  // Set the color palette
  gStyle->SetPalette(1);

  // Draw the graph, set top view
  _graph->Draw("PCOLZ");
  _graph->SetMarkerStyle(20);
  _graph->SetMarkerSize(1);
  _graph->GetZaxis()->SetLabelFont(42);
  _graph->GetZaxis()->SetTitleFont(42);
  gPad->SetTheta(90);
  gPad->SetPhi(0);

  // Make the basic default axes disappear as they are misplaced by default
  _graph->GetXaxis()->SetTitleOffset(999);
  _graph->GetXaxis()->SetLabelOffset(999);  
  _graph->GetXaxis()->SetTickSize(0);  

  _graph->GetYaxis()->SetTitleOffset(999);
  _graph->GetYaxis()->SetLabelOffset(999);  
  _graph->GetYaxis()->SetTickSize(0);

  // Define new axes and draw them
  double wmin = _graph->GetXaxis()->GetXmin();
  double wmax = _graph->GetXaxis()->GetXmax();
//  size_t ndiv = _graph->GetXaxis()->GetNdivisions();
  _xaxis->SetLabelFont(42);
  _xaxis->SetTextFont(42);
  _xaxis->DrawAxis(-.5775, -.5775, .5775, -.5775, wmin, wmax, 505);

  wmin = _graph->GetYaxis()->GetXmin();
  wmax = _graph->GetYaxis()->GetXmax();
//  ndiv = _graph->GetYaxis()->GetNdivisions();
  _yaxis->SetLabelFont(42);
  _yaxis->SetTextFont(42);
  _yaxis->DrawAxis(-.5775, -.5775, -.5775, .5775, wmin, wmax, 505);
}

void ScatterGraph::Save() {

  // Initialize a canvas, draw the graph on it and save it
  TCanvas *c = new TCanvas("c", "c", 1200, 800);
  gPad->SetLogz();
  Draw();
  c->SaveAs(_graph->GetName()+TString(".pdf"));
  delete c;
}
