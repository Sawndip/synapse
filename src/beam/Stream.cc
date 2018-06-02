#include "Stream.hh"

namespace Beam {

Stream::Stream() :
  _bunches(), _plane_ids(), _reaches(), _fraction(1.) {
}

Stream::Stream(const std::map<size_t, Bunch>& bunches,
	       const std::vector<size_t> reaches) :
  _bunches(bunches), _plane_ids(), _reaches(reaches), _fraction(1.) {

  try {
    Initialize();
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "Stream::Stream"));
  }
}

Stream::Stream(const Stream& stream) {
  *this = stream;
}

Stream& Stream::operator=(const Stream& stream) {
  if ( this == &stream )
      return *this;

  _bunches = stream._bunches;
  _plane_ids = stream._plane_ids;
  _reaches = stream._reaches;
  _fraction = stream._fraction;

  return *this;
}

Stream::~Stream () {}

const Bunch& Stream::Front() const {

  double minpos(DBL_MAX);
  size_t id;
  for (const size_t& i : _plane_ids)
    if ( _bunches.at(i).Position() < minpos ) {
      minpos = _bunches.at(i).Position();
      id = i;
    }

  return _bunches.at(id);
}

const Bunch& Stream::Back() const {

  double maxpos(-DBL_MAX);
  size_t id;
  for (const size_t& i : _plane_ids)
    if ( _bunches.at(i).Position() > maxpos ) {
      maxpos = _bunches.at(i).Position();
      id = i;
    }

  return _bunches.at(id);
}

double Stream::Transmission(const size_t plane_id) const {

  // Find the amount of particles in the first bunch
  size_t size = Front().Size();

  // Get the value of the transmission
  return (double)_bunches.at(plane_id).Size()/size;
}

void Stream::SetCoreFraction(const double frac) {

  // Find the amount of particles to be selected in the core
  size_t size = Front().Size()*frac;
  _fraction = frac;

  // Execute
  SetCoreSize(size);
}

TGraphErrors* Stream::EvolutionGraph(const SumStat& stat) const {

  try {
    // If the transmission is requested, call the relevant function
    if ( stat == trans )
	return TransmissionGraph();

    // Initialize the graph
    TGraphErrors* graph = new TGraphErrors();

    // Get the values of the varibale at each of the z position
    size_t id(0);
    double pos;
    Variable val;
    for (const std::pair<size_t, Bunch> bunch : _bunches) {
      pos = bunch.second.Position();
      val = bunch.second.SummaryStatistic(stat);
      graph->SetPoint(id, pos, val.GetValue());
      graph->SetPointError(id, 1, val.GetError());
      id++;
    }

    // Set titles
    graph->SetName(SumStatDict[stat].name.c_str());
    graph->GetXaxis()->SetTitle("z [mm]");
    std::string axis_title = SumStatDict[stat].label;
    if ( SumStatDict[stat].frac )
        axis_title += "_{"+std::to_string((int)(100*_fraction))+"}";
    if ( SumStatDict[stat].unit.size() )
        axis_title += "  ["+SumStatDict[stat].unit+"]";
    graph->GetYaxis()->SetTitle(axis_title.c_str());

    return graph;
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::recoverable,
	  "Could not produce a evolution graph"+std::string(e.what()),
	  "Stream::EvolutionGraph"));
  }
}

TGraphErrors* Stream::TransmissionGraph() const {

  try {
    // Initialize the graph
    TGraphErrors* graph = new TGraphErrors();

    // Find the amount of particles in the first bunch
    size_t size = Front().Size();

    // Get the values of the transmission at each of the z position
    size_t id(0);
    double pos, trans, err;
    for (const std::pair<size_t, Bunch> bunch : _bunches) {
      pos = bunch.second.Position();
      trans = (double)bunch.second.Size()/size;
      err = trans*sqrt((1.-trans)/(trans*size));
      graph->SetPoint(id, pos, 100*trans);
      graph->SetPointError(id, 1, 100*err);
      id++;
    }

    // Set titles
    graph->SetName("trans");
    graph->GetXaxis()->SetTitle("z [mm]");
    graph->GetYaxis()->SetTitle("Transmission [%]");

    return graph;
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::recoverable,
	  "Could not produce a transmission graph"+std::string(e.what()),
	  "Stream::TransmissionGraph"));
  }
}

std::map<SumStat, TGraphErrors*> Stream::FractionalGraphs(const size_t idu,
						    	  const size_t idd) {

  try {
    // Check that the stream contains the requested plane ids
    if ( !Contains(idu) || !Contains(idd) )
        throw(Exceptions::Exception(Exceptions::nonRecoverable,
	      "One of the plane ids involved in the comparison is missing",
	      "Stream::FractionalGraphs"));

    // Checks that idu is strictly upstream of idd
    if ( _bunches.at(idu).Position() >= _bunches.at(idd).Position() )
        throw(Exceptions::Exception(Exceptions::nonRecoverable,
	      "The downstream reference plane is upstream",
	      "Stream::FractionalGraphs"));

    // Initialize the graphs
    std::map<SumStat, TGraphErrors*> graphs;
    for (const SumStat& stat : {amp, subeps, vol})
        graphs[stat] = new TGraphErrors();

    // Loop over increasing core fractions
    double fracmin(.05), fracstep(0.025), fracmax(.5);
    Variable vu, vd;
    double change, error;
    size_t id(0), size;
    for (double frac = fracmin; frac < fracmax+1e-6; frac += fracstep) {
      // If the fraction is larger than the transmission, abort
      if ( frac >= Transmission(idu) || frac >= Transmission(idd) )
	  break;

      // Set the current core fraction, compute the volumes
      size = frac*_bunches.at(idu).Size();
      _bunches.at(idu).SetCoreSize(size);
      _bunches.at(idd).SetCoreSize(size);

      _bunches.at(idu).SetCoreVolume();
      _bunches.at(idd).SetCoreVolume();

      // Compute the relative change, fill the graphs
      for (const SumStat& stat : {amp, subeps, vol}) {
	vu = _bunches.at(idu).SummaryStatistic(stat);
	vd = _bunches.at(idd).SummaryStatistic(stat);

	change = 100.*(vd.GetValue()/vu.GetValue()-1.);
	error = change*sqrt(pow(vu.GetError()/vu.GetValue(), 2)+pow(vd.GetError()/vd.GetValue(), 2));

	graphs[stat]->SetPoint(id, 100.*frac, change);
	graphs[stat]->SetPointError(id, 0., error);
      }
      id++;
    }

    // Set titles
    for (const SumStat& stat : {amp, subeps, vol}) {
      graphs[stat]->SetName((SumStatDict[stat].name+"_frac").c_str());
      graphs[stat]->GetXaxis()->SetTitle("Fraction #alpha [%]");
      std::string axis_title = "#Delta"+SumStatDict[stat].label+"_{#alpha}/"+
				SumStatDict[stat].label+"_{#alpha}^{u}  [%]";
      graphs[stat]->GetYaxis()->SetTitle(axis_title.c_str());
    }

    return graphs;
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::recoverable,
	  "Could not produce a evolution graph"+std::string(e.what()),
	  "Stream::EvolutionGraph"));
  }
}

TF1* Stream::FractionalFunction(const SumStat& stat) const {

  // Throw if it is not a fractional quantity
  if ( !SumStatDict[stat].frac )
      throw(Exceptions::Exception(Exceptions::recoverable,
	    "Not a fractional quantity: "+SumStatDict[stat].name,
	    "Bunch::FractionalFunction"));

  // Get the dimensionality of the beam and mass of the particles
  double mass = Front().Mass();
  double dim = Front().Dimension();

  // Set the parameters
  TF1* func = NULL;
  double R2 = TMath::ChisquareQuantile(_fraction, dim);
  if ( stat == amp ) {
    func = new TF1(("f"+SumStatDict[stat].name).c_str(), "[0]*x", 0, 100);
    func->SetParameter(0, R2);
  } else if ( stat == subeps ) {
    func = new TF1(("f"+SumStatDict[stat].name).c_str(), "[0]*x", 0, 100);
    func->SetParameter(0, TMath::Gamma(3., R2/2.)/_fraction);
  } else if ( stat == vol ) {
    func = new TF1(("f"+SumStatDict[stat].name).c_str(),
	"pow(TMath::Pi()*[0]*[1]*x, [2]/2)/TMath::Gamma([2]/2+1)", 0, 100);
    func->SetParameters(mass, R2, dim);
  }

  return func;
}

void Stream::Initialize() {

  // Set the list of plane ids that compose the stream
  for (const std::pair<size_t, Bunch>& bunch : _bunches)
      _plane_ids.push_back(bunch.first);
}

void Stream::SetCoreSize(const size_t size) {

  try {
    // Set the core size in each of the bunches
    for (const size_t& i : _plane_ids) {
      _bunches[i].SetFraction(_fraction);
      _bunches[i].SetCoreSize(size);
    }

  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not set"+std::string(e.what()),
	  "Stream::SubSample"));
  }
}

bool Stream::Contains(const size_t plane_id) const {

  return _bunches.find(plane_id) != _bunches.end();
}
} // namespace Beam
