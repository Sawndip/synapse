#include "InfoBox.hh"

InfoBox::InfoBox(const std::string type,
	         const std::string maus,
	         const std::string run,
	         const std::string cycle) :
  _text(), _type(type), _maus(maus), _run(run), _cycle(cycle) {

  // Set it in the top right corner by default
  _text = new TPaveText(.7, .7, .875, .875, "NDC");
  _text->SetLineColorAlpha(0, 0);
  _text->SetFillStyle(1000);
  _text->SetFillColor(0);
  _text->SetTextAlign(12);
  _text->SetTextFont(42);
}

InfoBox::InfoBox(const InfoBox& ib) {
  *this = ib;
}

InfoBox& InfoBox::operator=(const InfoBox& ib) {
  if ( this == &ib )
      return *this;

  if ( _text )
      delete _text;
  _text = ib._text;
  _type = ib._type;
  _maus = ib._maus;
  _run = ib._run;
  _cycle = ib._cycle;

  return *this;
}

InfoBox::~InfoBox () {

//  delete _text;
}

void InfoBox::SetPosition(const std::string& pos,
			  const double alpha,
			  const double beta) {

  // Get the limits of the text box from the pos variable
  double xmin = pos.find("r") < pos.size() ? .9-alpha-beta : .1+alpha;
  double xmax = pos.find("r") < pos.size() ? .9-alpha : .1+alpha+beta;
  double ymin = pos.find("t") < pos.size() ? .9-alpha-beta : .1+alpha;
  double ymax = pos.find("t") < pos.size() ? .9-alpha : .1+alpha+beta;

  if ( _text )
      delete _text;
  _text = new TPaveText(xmin, ymin, xmax, ymax, "NDC");
  _text->SetLineColorAlpha(0, 0);
  _text->SetFillStyle(1000);
  _text->SetFillColor(0);
  _text->SetTextAlign(12);
  _text->SetTextFont(42);
}

void InfoBox::SetPosition(const double& xmin,
		   	  const double& ymin,
		   	  const double& xmax,
		   	  const double& ymax) {

  if ( _text )
      delete _text;
  _text = new TPaveText(xmin, ymin, xmax, ymax, "NDC");
  _text->SetLineColorAlpha(0, 0);
  _text->SetFillStyle(1000);
  _text->SetFillColor(0);
  _text->SetTextAlign(12);
  _text->SetTextFont(42);
}

void InfoBox::SetOpacity(const double& opacity) {

  _text->SetFillColorAlpha(0, opacity);
}

void InfoBox::Draw() {

  _text->Clear();

  TText* t = _text->AddText(TString::Format("#scale[1.5]{MICE} %s", _type.c_str()));
  t->SetTextFont(62);
  if ( _cycle.size() )
      _text->AddText(TString::Format("ISIS Cycle %s", _cycle.c_str()));
  if ( _run.size() )
      _text->AddText(TString::Format("Run setting %s", _run.c_str()));
  if ( _maus.size() )
      _text->AddText(TString::Format("MAUS v%s", _maus.c_str()));

  _text->Draw("SAME");
}
