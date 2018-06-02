#include "Aperture.hh"

namespace Beam {

Aperture::Aperture() :
  _apertures({-DBL_MAX, DBL_MAX, DBL_MAX}), _starts(0), _ends(0), _radii(0), _names(0) {
}

Aperture::Aperture(std::vector<double> starts,
		   std::vector<double> ends,
		   std::vector<double> radii,
		   std::vector<std::string> names) :
  _apertures(0), _starts(starts), _ends(ends), _radii(radii), _names(names) {

  try {
    Initialize(starts, ends, radii, names);
  } catch ( Exceptions::Exception& e ) {
    throw(Exceptions::Exception(Exceptions::nonRecoverable,
	  "Could not initialize"+std::string(e.what()),
	  "Aperture::Aperture"));
  }
}

Aperture::Aperture(const Aperture& aperture) {
  *this = aperture;
}

Aperture& Aperture::operator=(const Aperture& aperture) {
  if ( this == &aperture )
      return *this;

  _apertures = aperture._apertures;
  _starts = aperture._starts;
  _ends = aperture._ends;
  _radii = aperture._radii;
  _names = aperture._names;

  return *this;
}

Aperture::~Aperture () {}

void Aperture::Initialize(std::vector<double> starts,
		   	  std::vector<double> ends,
		   	  std::vector<double> radii,
		   	  std::vector<std::string> names) {

  // Check that all the vectors are of the same dimensions (accept 0)
  Assert::SameSize("Aperture::Initialize", "Starts and ends vectors", starts, ends);
  Assert::SameSize("Aperture::Initialize", "Starts, ends and radii vectors", starts, radii);
  if ( names.size() && names.size() != starts.size() )
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	    "Wrong amount of names provided",
	    "Aperture::Initialize"));

  // Initially, there are no apertures, i.e. DBL_MAX limit everywhere
  _apertures = {-DBL_MAX, DBL_MAX, DBL_MAX};

  // Loop over the apertures, edit the apertures vector accordingly
  size_t N = starts.size();
  size_t i;
  for (i = 0; i < N; i++)
    if ( names.size() ) {
      Add(starts[i], ends[i], radii[i], names[i]);
    } else {
      Add(starts[i], ends[i], radii[i]);
    }
}

void Aperture::SetMICEDefault(const std::string& filename) {

  // Get the position of the trackers
  GeometryHandler geoh(filename, {"tku", "tkd"});

  // Magnet bore (200\,mm inner radius of the vacuum vessel)
  Add(geoh["tku"].z()-1e3, geoh["tkd"].z()+1e3, 200.);

  // Tracker fiducials
  Add(geoh["tku"][4].z(), geoh["tku"][0].z(), 150.);
  Add(geoh["tkd"][0].z(), geoh["tkd"][4].z(), 150.);

  // TODO Add the MMB apertures and other windows
}

void Aperture::Add(double start, double end, double radius, std::string name) {

  // Add the aperture to the current list
  _starts.push_back(start);
  _ends.push_back(end);
  _radii.push_back(radius);
  if ( name.size() ) {
      _names.push_back(name);
  } else {
      _names.push_back("");
  }

  // Loop over the current limits in the aperture vector
  std::vector<double> temp = {_apertures.front()};
  size_t i, j;
  for (i = 0; i < _apertures.size(); i+=2) {
    if ( start > _apertures[i] && start < _apertures[i+2] ) {

      // If the radius is smaller, the initial aperture ends at the start
      if ( radius < _apertures[i+1] ) {
	temp.push_back(_apertures[i+1]);
	temp.push_back(start);
      }

      // If the start is contained, find the end, add aperture when necessary
      for (j = i+2; j < _apertures.size(); j+=2) {
        if ( end > _apertures[j] ) {
	  // If the end is further than this limit, check before the limit
 	  // If the radius is larger, add current, if it is smaller, keep looking
	  if ( radius > _apertures[j-1] ) {
	    temp.push_back(_apertures[j-1]);
	    temp.push_back(_apertures[j]);
          }
        } else {
	  if ( radius < _apertures[j-1] ) {
	    temp.push_back(radius);
	    temp.push_back(end);
          }
	  i = j - 2;	// Update to not go over the same apertures twice
	  break;
        }
      }
    } else if ( i == _apertures.size()-1 || end < _apertures[i] ) {

      // If after the end point, just add the current aperture 
      temp.push_back(_apertures[i-1]);
      temp.push_back(_apertures[i]);
    } else if ( _apertures[i+2] < start ) {

      // If before the starting point, just add the current aperture 
      temp.push_back(_apertures[i+1]);
      temp.push_back(_apertures[i+2]);
    }
  }

  _apertures = temp;
}

const double& Aperture::Radius(double z) const {

  size_t i;
  for (i = 0; i < _apertures.size()-2; i+=2)
    if ( z < _apertures[i+2] )
        return _apertures[i+1];

  return _apertures.back();
}

bool Aperture::IsIn(double x, double y, double z) const {

  return sqrt(pow(x, 2)+pow(y, 2)) < Radius(z);
}

void Aperture::Draw(std::string name, double dz) const {

  // If there are less than 7 elements, there are no apertures
  if ( _apertures.size() < 7 ) {
    Pitch::print(Pitch::warning, "No apertures added to the space, nothing to draw");
    return;
  }

  // Initialize the TGraph
  TGraph graph;
  graph.SetTitle("Apertures");
  graph.SetMarkerStyle(4);
  
  // Loop from the start of the first to the end of the last aperture
  size_t N = _apertures.size();
  size_t i(0);
  double z, radius;
  double temp = Radius(_apertures[2]);
  for (z = _apertures[2]-dz; z < _apertures[N-3]+dz; z += dz) {
    radius = std::min(Radius(z), 1000.);
    if ( radius != temp ) {
      graph.SetPoint(i, z, temp);
      i++;
      graph.SetPoint(i, z, radius);
      i++;
      temp = radius;
    }
  }

  // If there are to the apertures, add them to the canvas
  std::vector<TText*> labels;
  for (i = 0; i < _names.size(); i++)
    if ( _names[i].size() )
      labels.push_back(new TText((_starts[i]+_ends[i])/2, _radii[i], _names[i].c_str()));

  // Draw on a canvas, clean up
  TCanvas *canv = new TCanvas("c", "c", 1200, 800);
  graph.Draw("APL");
  for (TText* label : labels)
      label->Draw("SAME");
  graph.GetXaxis()->SetTitle("z [mm]");
  graph.GetYaxis()->SetTitle("R [mm]");

  std::string cname = (bool)name.size() ? "apertures_"+name+".pdf" : "apertures.pdf";
  canv->SaveAs(cname.c_str());
  delete canv;
  for (TText* label : labels)
      delete label;
}
} // namespace Beam
