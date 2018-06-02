#include "DFunction.hh"

DFunction::DFunction() :
  _dim(0), _lower(0), _upper(0), _name(""), _title(""), _rdmzer() {
}

DFunction::DFunction(const DFunction& df) {
  *this = df;
}

DFunction& DFunction::operator=(const DFunction& df) {
  if ( this == &df )
      return *this;

  _dim = df._dim;
  _lower = df._lower;
  _upper = df._upper;
  _name = df._name;
  _title = df._title;
  _rdmzer = df._rdmzer;

  return *this;
}

DFunction::~DFunction () {}

double DFunction::Evaluate(const double& v) const {

  return Evaluate(std::vector<double>({v}));
}

double DFunction::Evaluate(const double* v) const {

  std::vector<double> vv(v, v+_dim);
  return Evaluate(vv);
}

void DFunction::SetLowerBound(const double lower, const size_t i) {

  Assert::IsWithinBounds("DFunction::SetLowerBound", _lower.size(), i);
  _lower[i] = lower;
}

void DFunction::SetUpperBound(const double upper, const size_t i) {

  Assert::IsWithinBounds("DFunction::SetUpperBound", _upper.size(), i);
  _upper[i] = upper;
}

void DFunction::Range(std::vector<double>& lower,
	     	      std::vector<double>& upper) const {

  lower = _lower;
  upper = _upper;
}

void DFunction::SetRange(const double lower,
	      		 const double upper,
			 const size_t i) {

  Assert::IsWithinBounds("DFunction::SetRange", _lower.size(), i);
  if ( lower >= upper )
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	    "The boundaries define a zero or negative interval",
	    "DFunction::SetRange"));

  _lower[i] = lower;
  _upper[i] = upper;
}

void DFunction::SetRange(const std::vector<double>& lower,
	      		 const std::vector<double>& upper) {

  Assert::SameSize("DFunction::SetRange", "Lower and upper boundaries", lower, upper);
  Assert::SameSize("DFunction::SetRange", "Function space and boundaries", _lower, lower);
  for (size_t i = 0; i < lower.size(); i++)
    if ( lower[i] >= upper[i] )
        throw(Exceptions::Exception(Exceptions::nonRecoverable,
	      "The boundaries of axis "+std::to_string(i)+" define a zero or negative interval",
	      "DFunction::SetRange"));

  _lower = lower;
  _upper = upper;
}

TGraph* DFunction::Graph(size_t idx, std::vector<double> x) const {

  // Make sure the the requested id is smaller than the dimension
  Assert::IsGreater("DFunction::Graph", "Dimension", _dim, idx, true);

  // Make sure the vector is of the right dimension. If not provided, fill with 0
  if ( !x.size() )
      x = std::vector<double>(_dim, 0);
  Assert::IsEqual("DFunction::Graph", "Dimension and vector size", _dim, x.size());
  
  // For a certain amount of points, initialize a TGraph
  const size_t npoints = 100;
  TGraph *graph = new TGraph(npoints);
  graph->SetTitle(_title.c_str());
  graph->GetXaxis()->SetTitle(TString::Format("x_{%d}", (int)idx));

  // Fill a TGraph with the function values
  size_t i;
  for (i = 0; i < npoints; i++) {
    x[idx] = _lower[idx]+i*(_upper[idx]-_lower[idx])/(npoints-1);
    graph->SetPoint(i, x[idx], Evaluate(x));
  }
    
  // Set the style and return
  graph->SetLineWidth(2);
  graph->SetLineColor(2);
  return graph;
}

TGraph* DFunction::GraphRadial() const {

  try {
  // For a certain amount of points, initialize a TGraph
    const size_t npoints = 100;
    TGraph *graph = new TGraph(npoints);
    graph->SetTitle(_title.c_str());
    graph->GetXaxis()->SetTitle("R");

  // Fill a TGraph with the radial function values
    double r;
    size_t i;
    for (i = 0; i < npoints; i++) {
      r = i*_upper[0]/(npoints-1);
      graph->SetPoint(i, r, Radial(r));
    }
    
    // Set the style and return
    graph->SetLineWidth(2);
    graph->SetLineColor(2);
    return graph;
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not compute the radial function"+std::string(e.what()),
	  "DFunction::GraphRadial"));
  }
}

TH2F* DFunction::Graph2D(size_t idx, size_t idy, std::vector<double> x) const {

  // Make sure the the requested ids are smaller than the dimension
  Assert::IsGreater("DFunction::Graph2D", "Dimension", _dim, idx, true);
  Assert::IsGreater("DFunction::Graph2D", "Dimension", _dim, idy, true);

  // Make sure the vector is of the right dimension. If not provided, fill with 0
  if ( !x.size() )
      x = std::vector<double>(_dim, 0);
  Assert::IsEqual("DFunction::Graph2D", "Dimension and vector size", _dim, x.size());
  
  // For a certain amount of points, compute the function and fill a TH2F
  const size_t npoints = 100;
  TH2F *graph = new TH2F(_name.c_str(), _title.c_str(), npoints, _lower[idx],
			 _upper[idx], npoints, _lower[idy], _upper[idy]);
  graph->GetXaxis()->SetTitle(TString::Format("x_{%d}", (int)idx));
  graph->GetYaxis()->SetTitle(TString::Format("x_{%d}", (int)idy));
  size_t i, j;
  for (i = 0; i < (size_t)graph->GetNbinsX(); i++) {
    for (j = 0; j < (size_t)graph->GetNbinsY(); j++) {
      x[idx] = graph->GetXaxis()->GetBinCenter(i+1);
      x[idy] = graph->GetYaxis()->GetBinCenter(j+1);
      graph->SetBinContent(i+1, j+1, Evaluate(x));
    }
  }

  // Return
  return graph;
}

TGraph* DFunction::GraphCDFRadial() const {

  try {
  // For a certain amount of points, initialize a TGraph
    const size_t npoints = 100;
    TGraph *graph = new TGraph(npoints);
    graph->SetTitle(_title.c_str());
    graph->GetXaxis()->SetTitle("R");

  // Fill a TGraph with the radial function values
    double r;
    size_t i;
    for (i = 0; i < npoints; i++) {
      r = i*_upper[0]/(npoints-1);
      graph->SetPoint(i, r, CDFRadial(r));
    }
    
    // Set the style and return
    graph->SetLineWidth(2);
    graph->SetLineColor(2);
    return graph;
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not compute the radial CDF"+std::string(e.what()),
	  "DFunction::GraphRadial"));
  }
}

void DFunction::Draw(const std::string opt, int idx, int idy, std::vector<double> x) const {

  // Sort dimension wise and produce the relevant functions
  if ( _dim == 1 || idy < 0 ) {
    TGraph* graph = Graph(idx, x);
    graph->Draw(opt.c_str());
    return;
  }
  
  TH2F* graph = Graph2D(idx, idy, x);
  graph->Draw(opt.c_str());
  graph->SetLineColor(2);
  graph->SetLineWidth(2);
  return;
}

void DFunction::DrawRadial(const std::string opt) const {

  try {
    TGraph* graph = GraphRadial();
    graph->Draw(opt.c_str());
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not draw the radial function"+std::string(e.what()),
	  "DFunction::DrawRadial"));
  }
}

void DFunction::Print(const std::string opt, int idx, int idy, std::vector<double> x) const {

  try {
    gStyle->SetOptStat(0);
    TCanvas *c = new TCanvas("c", "c", 1200, 800);
    Draw(opt, idx, idy, x);
    c->SaveAs(TString::Format("%s.pdf", _name.c_str()));
    delete c;
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not print the function"+std::string(e.what()),
	  "DFunction::Print"));
  }
}

void DFunction::PrintRadial(const std::string opt) const {

  try {
    gStyle->SetOptStat(0);
    TCanvas *c = new TCanvas("c", "c", 1200, 800);
    DrawRadial(opt);
    c->SaveAs(TString::Format("%s_radial.pdf", _name.c_str()));
    delete c;
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not print the radial function"+std::string(e.what()),
	  "DFunction::PrintRadial"));
  }
}
