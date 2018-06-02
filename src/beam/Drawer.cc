#include "Drawer.hh"

namespace Beam {

Drawer::Drawer() :
  _types(), _vars(), _options(), _varsets(), _info(NULL) {

  try {
    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "Drawer::Drawer"));
  }
}

Drawer::Drawer(const Drawer& dra) {
  *this = dra;
}

Drawer& Drawer::operator=(const Drawer& dra) {
  if ( this == &dra )
      return *this;

  _types = dra._types;
  _vars = dra._vars;
  _options = dra._options;
  _varsets = dra._varsets;
  if ( _info )
      delete _info;
  _info = new InfoBox(*dra._info);

  return *this;
}

Drawer::~Drawer () {

  if ( _info )
      delete _info;
}

void Drawer::SetVariableLimits(const std::string& var,
		 	       double llim,
	         	       double ulim) {
  try {
    _varsets.at(var).llim = llim;
    _varsets.at(var).ulim = ulim;

  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not set the limits of "+var+std::string(e.what()),
	  "Drawer::SetVariableLimits"));
  }
}

TH1F* Drawer::New1DHistogram(const std::string& name,
		     	     const std::string& var) const {

  VarSet set = _varsets.at(var);
  std::string x_title = set.label+" ["+set.unit+"]";

  TH1F* h = new TH1F((name+"_"+var).c_str(), (";"+x_title).c_str(), set.nbins, set.llim, set.ulim);
  return h;
}

TH2F* Drawer::New2DHistogram(const std::string& name,
		       	     const std::string& varx,
		     	     const std::string& vary) const {

  VarSet setx = _varsets.at(varx);
  VarSet sety = _varsets.at(varx);
  std::string x_title = setx.label+" ["+setx.unit+"]";
  std::string y_title = sety.label+" ["+sety.unit+"]";

  TH2F* h = new TH2F((name+"_"+varx+vary).c_str(), (";"+x_title+";"+y_title).c_str(),
		     setx.nbins, setx.llim, setx.ulim, sety.nbins, sety.llim, sety.ulim);
  return h;
}

ScatterGraph* Drawer::NewAmplitudeScatter(const std::string& name,
		       	      	      	  const std::string& varx,
		       	      	      	  const std::string& vary) const {

  VarSet setx = _varsets.at(varx);
  VarSet sety = _varsets.at(varx);
  std::string x_title = setx.label+" ["+setx.unit+"]";
  std::string y_title = sety.label+" ["+sety.unit+"]";

  ScatterGraph* s = new ScatterGraph();
  s->SetName((name+"_"+varx+vary).c_str());
  s->SetXaxisTitle(x_title.c_str());
  s->SetYaxisTitle(y_title.c_str());
  s->GetZaxis()->SetTitle("A_{#perp}  [mm]");
  return s;
}

THStack* Drawer::NewStack(std::map<std::string, TH1F*> hists) const {

  // Set the style of each histogram according to its data type,
  // add them to a stack with the appropriate drawing options
  THStack* hs = new THStack();
  for (const std::string& type : {"recmc", "truth", "data"})
    if ( hists.find(type) != hists.end() ) {
      SetStyle(hists[type], type);
      hs->Add(hists[type], _options.at(type).hist_opt);
    }

  // Set the axes of the stack
  TH1F* hist = (*hists.begin()).second;
  hs->SetTitle(TString::Format(";%s;%s;Density",
	hist->GetXaxis()->GetTitle(), hist->GetYaxis()->GetTitle()));

  return hs;
}

TLegend* Drawer::NewStackLegend(std::map<std::string, TH1F*> hists,
				const double hoffset) const {

  // Build the legend
  double height = hists.size()*.08;
  TLegend *leg = new TLegend(.64, .88-height-hoffset, .9, .88-hoffset);
  leg->SetFillStyle(0);
  leg->SetLineColorAlpha(0, 0);
  for (const std::pair<std::string, TH1F*> hist : hists)
      leg->AddEntry(hist.second, _options.at(hist.first).title, _options.at(hist.first).hist_leg);

  return leg;
}

std::map<std::string, TPaveStats*> Drawer::NewStackStats(std::map<std::string, TH1F*> hists,
		       			       		 const std::string draw_opt,
			 		       		 double* hoffset) const {

  // If the statistics are requested, initialize them
  std::map<std::string, TPaveStats*> stats;
  double h = 0.;
  if ( draw_opt.size() ) {
    h = .05*(draw_opt.size()+1);
    size_t id(0);
    for (const std::pair<std::string, TH1F*> hist : hists) {
      stats[hist.first] = new TPaveStats(.64, .9-(id+1)*h, .9, .9-id*h, "NDC");
      stats[hist.first]->SetTextFont(42);
      stats[hist.first]->SetFillStyle(0);
       
      TText *title = stats[hist.first]->AddText(TString::Format("%s = ", _options.at(hist.first).title));
      title->SetTextFont(62);
      if ( draw_opt.find("E") != std::string::npos )
          stats[hist.first]->AddText(TString::Format("Entries = %0.0lf", hist.second->GetEntries()));
      if ( draw_opt.find("M") != std::string::npos )
          stats[hist.first]->AddText(TString::Format("Mean = %0.2lf #pm %0.2lf",
					hist.second->GetMean(), hist.second->GetMeanError()));
      if ( draw_opt.find("R") != std::string::npos )
          stats[hist.first]->AddText(TString::Format("RMS = %0.2lf #pm %0.2lf",
					hist.second->GetRMS(), hist.second->GetRMSError()));
      id++;
    }
  }

  if ( hoffset != NULL ) 
      *hoffset = h*hists.size();

  return stats;
}

THStack* Drawer::DrawStack(std::map<std::string, TH1F*> hists,
		      	   const std::string& name) const {

  // Draw the THStack
  THStack* hs = NewStack(hists);
  hs->SetName(name.c_str());
  hs->Draw("NOSTACK");

  // Draw the stat boxes
  double hoffset = 0;
  for ( std::pair<std::string, TPaveStats*> stat : NewStackStats(hists, "EM", &hoffset) )
      stat.second->Draw("SAME");

  // Draw the legend
  NewStackLegend(hists, hoffset)->Draw("SAME");

  // Draw the info box if defined
  if ( _info )
      _info->Draw();

  return hs;
}

void Drawer::SaveStack(std::map<std::string, TH1F*> hists,
		       const std::string& name,
		       bool ratio) const {

  // Initialize the TCanvas and partition it if requested
  TCanvas* c = new TCanvas("c", "c", 1200, 800);
  TPad *tpad, *bpad;
  if ( ratio && hists.find("data") != hists.end() && hists.find("recmc") != hists.end() ) {
    tpad = new TPad("tpad", "", 0., .25, 1., 1.);
    tpad->SetBottomMargin(0.);
    tpad->Draw();
    tpad->cd();
  }

  // Draw stack and surrounding info boxes
  THStack* hs = DrawStack(hists, name);

  // If requested, initialize the ratio and draw it on the second pad
  if ( ratio && hists.find("data") != hists.end() && hists.find("recmc") != hists.end() ) {
    c->cd();
    bpad = new TPad("bpad", "", 0., .0, 1., .25);
    bpad->SetTopMargin(0.);
    bpad->SetBottomMargin(0.3);
    bpad->Draw();
    bpad->cd();

    hs->SetMinimum(1e-4);

    TH1F* hist_ratio = (TH1F*)hists["data"]->Clone();
    hist_ratio->Divide(hists["recmc"]);
    hist_ratio->SetMinimum(0);
    hist_ratio->SetMaximum(2);
    hist_ratio->Draw("EP");

    double xmin = hist_ratio->GetXaxis()->GetXmin();
    double xmax = hist_ratio->GetXaxis()->GetXmax();
    double ymin = hist_ratio->GetMinimum();
    double ymax = hist_ratio->GetMaximum();

    hist_ratio->GetXaxis()->SetTitleSize(.105);
    hist_ratio->GetXaxis()->SetTickSize(0);
    hist_ratio->GetXaxis()->SetLabelOffset(999);
    TGaxis *xaxis = new TGaxis(xmin, ymin, xmax, ymin, xmin+1e-3, xmax-1e-3, 510, "");
    xaxis->SetLabelSize(0.105);
    xaxis->SetLabelFont(42);
    xaxis->Draw("SAME");

    hist_ratio->GetYaxis()->SetTitle("N_{D} /N_{S}");
    hist_ratio->GetYaxis()->SetTitleSize(.105);
    hist_ratio->GetYaxis()->SetTitleOffset(.3);
    hist_ratio->GetYaxis()->SetTickSize(0);
    hist_ratio->GetYaxis()->SetLabelOffset(999);
    TGaxis *yaxis = new TGaxis(xmin, ymin, xmin, ymax, ymin+1e-3, ymax-1e-3, 505, "");
    yaxis->SetLabelSize(0.105);
    yaxis->SetLabelFont(42);
    yaxis->Draw("SAME");

    TLine* line = new TLine(xmin, 1., xmax, 1.);
    line->Draw("SAME");
  }

  // Save the canvas
  c->SaveAs(TString::Format("%s.pdf", name.c_str()));
  delete c;
}

TMultiGraph* Drawer::NewMultiGraph(std::map<std::string, TGraphErrors*> graphs) const {

  // Set the style of each graph according to its data type,
  // add them to a multigraph with the appropriate drawing options
  TMultiGraph* mg = new TMultiGraph();
  for (const std::string& type : {"truth", "recmc", "data"})
    if ( graphs.find(type) != graphs.end() ) {
      SetStyle(graphs[type], type);
      mg->Add(graphs[type], _options.at(type).graph_opt);
    }

  // Set the axes of the stack
  TGraphErrors* graph = (*graphs.begin()).second;
  mg->SetTitle(TString::Format(";%s;%s",
	graph->GetXaxis()->GetTitle(), graph->GetYaxis()->GetTitle()));

  return mg;
}

TLegend* Drawer::NewMultiGraphLegend(std::map<std::string, TGraphErrors*> graphs,
				     const std::string corner) const {

  // Build the legend
  std::map<std::string, double> xmax = {{"tl",.475}, {"tr",.725}, {"bl",.475}, {"br",.725}};
  std::map<std::string, double> ymax = {{"tl",.88}, {"tr",.88}, {"bl",.28}, {"br",.28}};
  double height = graphs.size()*.16/3.;
  TLegend *leg = new TLegend(xmax[corner]-.2, ymax[corner]-height, xmax[corner], ymax[corner]);
  leg->SetFillStyle(1000);
  leg->SetLineColorAlpha(0, 0);
  for (const std::pair<std::string, TGraphErrors*> graph : graphs)
      leg->AddEntry(graph.second, _options.at(graph.first).title, _options.at(graph.first).graph_leg);

  return leg;
}

TGaxis* Drawer::NewMultiGraphEmittanceAxis(TF1* func,
				    	   const double& miny,
				      	   const double& maxy,
					   const double& maxx) const {

  func->SetRange(func->GetX(miny), func->GetX(maxy));
  TGaxis *axis = new TGaxis(maxx, miny, maxx, maxy, func->GetName(), 510, "+L");
  axis->SetTitle("#hat{#epsilon}_{#perp}  [mm]");
  axis->SetLabelFont(42);
  axis->SetTextFont(42);

  return axis;
}

TMultiGraph* Drawer::DrawMultiGraph(std::map<std::string, TGraphErrors*> graphs,
		      	   	    const std::string& name,
			      	    TF1* func) {

  // Draw the TMultiGraph
  TMultiGraph* mg = NewMultiGraph(graphs);
  mg->SetName(name.c_str());
  mg->Draw("A");

  // Draw the MICE beam line elements if it is the MICE beam line
  double miny(mg->GetYaxis()->GetXmin()), maxy(mg->GetYaxis()->GetXmax());
  if ( Globals::GetInstance()["mice"] )
      DrawMICEModules(miny, maxy);

  // Draw the emittance axis if the proportionality function is provided
  if ( func )
      NewMultiGraphEmittanceAxis(func, miny, maxy, gPad->GetUxmax())->Draw("SAME");

  // Find the best corner to drawn the stats
  std::string corner = FindBestCorner(graphs);

  // Draw the legend
  if ( corner.size() )
      NewMultiGraphLegend(graphs, corner)->Draw("SAME");

  // Draw the info box if defined
  if ( _info && corner.size() ) {
    _info->SetPosition(corner);
    _info->Draw();
  }

  return mg;
}

void Drawer::SaveMultiGraph(std::map<std::string, TGraphErrors*> graphs,
		    	    const std::string& name,
			    TF1* func) {

  // Initialize the TCanvas, set the Y grid
  TCanvas* c = new TCanvas("c", "c", 1200, 800);
  gPad->SetGridy();

  // Draw multigraph and surrounding info boxes
  DrawMultiGraph(graphs, name, func);

  // Save the canvas
  c->SaveAs(TString::Format("%s.pdf", name.c_str()));
  delete c;
}

TH1F* Drawer::ProjectionMeanX(TH2F* hist,
		      	      bool fit) const {

  TH1F* h1d = new TH1F(TString::Format("%s_1d", hist->GetName()),
		       TString::Format("%s;%s;%s", hist->GetTitle(), 
		       hist->GetXaxis()->GetTitle(), hist->GetYaxis()->GetTitle()),
		       hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());

  TF1* fgaus;
  if ( fit )
      fgaus = new TF1("fgaus", "[0]*TMath::Gaus(x, [1], [2])", -1e3, 1e3);

  size_t k;
  TH1D* hproj;
  double mean, error;
  for (k = 0; k < (size_t)hist->GetNbinsX(); k++) {
    hproj = hist->ProjectionY("bin", k+1, k+1);

    if ( !fit ) {
      mean = hproj->GetMean();
      error = hproj->GetMeanError();
    } else {
      fgaus->SetRange(hproj->GetMean()-hproj->GetRMS(), hproj->GetMean()+hproj->GetRMS());
      fgaus->SetParameters(hproj->GetMaximum(), hproj->GetMean(), hproj->GetRMS());
      hproj->Fit(fgaus, "RQ");
      mean = fgaus->GetParameter(1);
      error = fgaus->GetParError(1);
    }

    if ( error > mean )
	continue;

    h1d->SetBinContent(k+1, mean);
    h1d->SetBinError(k+1, error);
    delete hproj;
  }

  return h1d;
}

void Drawer::Initialize() {

  // Set the default drawing options
  gStyle->SetOptStat(0);

  // Set the known data types
  _types = {"truth", "recmc", "data"};

  // Set the known variables and their characteristics
  _vars = {"x", "y", "px", "py", "pz"};
  std::vector<std::string> labels = {"x", "y", "p_{x}", "p_{y}", "p_{z}"};
  std::vector<std::string> units = {"mm", "mm", "MeV/c", "MeV/c", "MeV/c"};
  std::vector<double> llims = {-250, -250, -100, -100, 0};
  std::vector<double> ulims = {250, 250, 100, 100, 250};
  for (size_t i = 0; i < _vars.size(); i++) {
    std::string var = _vars[i];
    _varsets[var].label = labels[i];
    _varsets[var].unit = units[i];
    _varsets[var].llim = llims[i];
    _varsets[var].ulim = ulims[i];
  }

  // Initialize the default style for each data type
  _options["truth"].title = "Truth";
  _options["truth"].line_color = kBlue+2;
  _options["truth"].fill_style = 1000;
  _options["truth"].fill_color = kBlue+2;
  _options["truth"].fill_alpha = .5;
  _options["truth"].marker_style = 27;
  _options["truth"].hist_opt = "hist";
  _options["truth"].graph_opt = "lpe3";
  _options["truth"].hist_leg = "f";
  _options["truth"].graph_leg = "lpf";

  _options["recmc"].title = "Simulation";
  _options["recmc"].fill_style = 1000;
  _options["recmc"].fill_color = kOrange-2;
  _options["recmc"].marker_style = 21;
  _options["recmc"].marker_color = kRed+2;
  _options["recmc"].hist_opt = "hist";
  _options["recmc"].graph_opt = "pe";
  _options["recmc"].hist_leg = "f";
  _options["recmc"].graph_leg = "lp";

  _options["data"].title = "Data";
  _options["data"].marker_style = 20;
  _options["data"].hist_opt = "pe";
  _options["data"].graph_opt = "pe";
  _options["data"].hist_leg = "lp";
  _options["data"].graph_leg = "lp";
}

void Drawer::SetStyle(TH1* hist,
	      	      const std::string& type) const {

  // Set the line style
  hist->SetLineStyle(_options.at(type).line_style);
  hist->SetLineColor(_options.at(type).line_color);
  hist->SetLineWidth(_options.at(type).line_width);

  // Set the fill style
  hist->SetFillStyle(_options.at(type).fill_style);
  hist->SetFillColorAlpha(_options.at(type).fill_color, _options.at(type).fill_alpha);

  // Set the marker style
  hist->SetMarkerStyle(_options.at(type).marker_style);
  hist->SetMarkerColor(_options.at(type).marker_color);
  hist->SetMarkerSize(_options.at(type).marker_size);
}

void Drawer::SetStyle(TGraph* graph,
	      	      const std::string& type) const {

  // Set the line style
  graph->SetLineStyle(_options.at(type).line_style);
  graph->SetLineColor(_options.at(type).line_color);
  graph->SetLineWidth(_options.at(type).line_width);

  // Set the fill style
  graph->SetFillStyle(_options.at(type).fill_style);
  graph->SetFillColorAlpha(_options.at(type).fill_color, _options.at(type).fill_alpha);

  // Set the marker style
  graph->SetMarkerStyle(_options.at(type).marker_style);
  graph->SetMarkerColor(_options.at(type).marker_color);
  graph->SetMarkerSize(2*_options.at(type).marker_size);
}

void Drawer::DrawMICEModules(const double& miny,
		       	     const double& maxy) const {

  // Create and draw TLines at the location of the tracker stations
  GeometryHandler geoh = Globals::GetInstance().GetGeometryHandler();
  for (const std::string det : {"tku", "tkd"}) { 
    for (size_t st = 0; st < 5; st++) {
	TLine *line = new TLine(geoh[det][st].z(), miny, geoh[det][st].z(), maxy);
	line->SetLineColorAlpha(4, .5);
	line->Draw("SAME");
    }
  }

  // Create and draw a TBox for the absorber
  double abs_dz = 65;
  TBox *abs_box = new TBox(geoh["abs"].z()-abs_dz/2, miny, geoh["abs"].z()+abs_dz/2, maxy);
  abs_box->SetLineColor(0);
  abs_box->SetFillStyle(1000);
  abs_box->SetFillColorAlpha(kGreen+2, .5);
  abs_box->Draw("SAME");
}

std::string Drawer::FindBestCorner(std::map<std::string, TGraphErrors*> graphs) const {

  // Loop over the corners, find one that does not have any points drawn. If it cannot be
  // found, abort and do not draw the stat and info boxes.
  std::vector<double> xllim = {.72, .12, .12, .72};
  std::vector<double> xulim = {.88, .28, .28, .88};
  std::vector<double> yllim = {.72, .72, .12, .12};
  std::vector<double> yulim = {.88, .88, .72, .72};
  std::vector<std::string> corners = {"tr", "tl", "bl", "br"};
  double xndc, yndc;
  size_t nin;
  for (size_t i = 0; i < corners.size(); i++) {
    nin = 0;
    for (const std::pair<std::string, TGraphErrors*>& element : graphs) {
      TGraphErrors* graph = element.second;
      for (size_t j = 0; j < (size_t)graph->GetN(); j++) {
	xndc = (graph->GetX()[j] - gPad->GetX1())/(gPad->GetX2()-gPad->GetX1());
	yndc = (graph->GetY()[j] - gPad->GetY1())/(gPad->GetY2()-gPad->GetY1());
        if ( xndc > xllim[i] && xndc < xulim[i] && yndc > yllim[i] && yndc < yulim[i] ) {
	  nin++;
	  break;
	}
      }
      if ( nin )
	  break;
    }
    if ( !nin )
	return corners[i];
  }

  return "";
}

void Drawer::SaveGIF(const std::string& name,
		     const std::vector<TCanvas*>& canv) const {

  // Silence ROOT to not flood the output
  gErrorIgnoreLevel = kWarning;

  // Print all the canvases
  Pitch::print(Pitch::debug, "Printing the individual canvases", "Drawer::SaveGIF");
  size_t i;
  ProgressBar pbar(Pitch::debug);
  for (i = 0; i < canv.size(); i++) {
    // Convert i in a 3 char string
    std::string id = std::to_string(i);
    while ( id.size() < 3 )
        id.insert(0, "0");

    // Save the current canvas
    canv[i]->SaveAs(TString::Format("%s_%s.gif", name.c_str(), id.c_str()));

    // Display the progress in %
    pbar.GetProgress(i, canv.size());
  }

  // Merge with ImageMagick, this is much faster than ROOT
  std::string sys_cmd = "convert -delay 10 -loop 0 "+name+"_*.gif "+name+".gif";
  if ( system(sys_cmd.c_str()) )
      throw(Exceptions::Exception(Exceptions::recoverable,
	    "Failed to merge the files into a single animated gif",
	    "Drawer::SaveGIF"));

  // Clean up
  sys_cmd = "rm "+name+"_*.gif";
  if ( system(sys_cmd.c_str()) )
      throw(Exceptions::Exception(Exceptions::recoverable,
	    "Failed to delete the temporary files used to produce an animated gif",
	    "Drawer::SaveGIF"));

  // Restore ROOT verbosity
  gErrorIgnoreLevel = kInfo;
}
} // namespace Beam
