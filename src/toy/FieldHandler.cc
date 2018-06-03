#include "FieldHandler.hh"

FieldHandler::FieldHandler() :
  _mu0(4*M_PI*1e-7) {
}

FieldHandler::FieldHandler(const FieldHandler& fh) {
  *this = fh;
}

FieldHandler& FieldHandler::operator=(const FieldHandler& fh) {
  if ( this == &fh )
      return *this;

  _mu0 = fh._mu0;
  _fields = fh._fields;

  return *this;
}

FieldHandler::~FieldHandler () {}

bool FieldHandler::AddSolenoid(const std::string& name,
		  	       const double& a,
		   	       const double& L,
		   	       const double& I) {

  // Create a field element
  FieldElement fel;
  fel.magtype = "solenoid";

  // Create a map with the characteristics of the solenoid
  std::map<std::string, double> chars;
  chars["a"] = a;
  chars["L"] = L;
  chars["I"] = I;
  fel.chars = chars;
  
  // Append the list of magnets
  if ( _fields.find(name) != _fields.end() )
      Pitch::print(Pitch::warning, "Overwriting field with name "+name, "FieldHandler::AddSolenoid");
  _fields[name] = fel;
  return true;
}

TVector3 FieldHandler::GetMagneticField(const double& x,
			    		const double& y,
			    		const double& z) const {

  // Loop over the available fields, fetch them and return
  TVector3 Bfield(0., 0., 0.);
  for ( auto &fel : _fields ) {
    if ( fel.second.magtype == "solenoid" ) {
      Bfield += GetSolenoidField(fel.second, x, y, z);
    } else {
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	    "Magnet of type \""+fel.second.magtype+"\" not supported",
	    "FieldHandler::GetMagneticField"));
    } 
  }

  return Bfield;
}

TVector3 FieldHandler::GetSolenoidField(const FieldElement& fel,
			    		const double& x,
			    		const double& y,
			    		const double& z) const {

  // Get the cylindrical coordinates from the x, y, z
  double rho = sqrt(x*x+y*y);
  double th(0.);
  if ( y )
      th = (.5+(y < 0))*M_PI - atan(x/y);
  th = (x < 0)*M_PI;

  // Get the cylindrical components of the field
  double Bperp = SolenoidBperp(fel.chars.at("a"), fel.chars.at("L"), fel.chars.at("I"), rho, z);
  double Blong = SolenoidBlong(fel.chars.at("a"), fel.chars.at("L"), fel.chars.at("I"), rho, z);

  // Convert to cartesian and return
  TVector3 Bcyl(Bperp, 0., Blong);
  Bcyl.RotateZ(th);
  return Bcyl;
}

double FieldHandler::SolenoidBperp(const double& a,
		       		   const double& L,
		       		   const double& I,
		       		   const double& rho,
		       		   const double& z) const {

  // The radial component of the field at rho=0 is 0
  if ( rho == 0 )
      return 0.;

  // Get the constants to use in the computation
  double zetap = z+L/2;
  double zetam = z-L/2;
  double kp = sqrt(4*a*rho/(pow(a+rho, 2)+pow(zetap, 2)));
  double km = sqrt(4*a*rho/(pow(a+rho, 2)+pow(zetam, 2)));

  // Compute the Bfield
  double fact = _mu0*I*sqrt(a/rho)/(2*M_PI*L);
  double plus = (kp*kp-2)*K(kp*kp)/kp+2*E(kp*kp)/kp;
  double minus = (km*km-2)*K(km*km)/km+2*E(km*km)/km;
  return fact*(plus-minus);
}

double FieldHandler::SolenoidBlong(const double& a,
		       		   const double& L,
		       		   const double& I,
		       		   const double& rho,
		       		   const double& z) const {

  // Get the constants to use in the computation
  double zetap = z+L/2;
  double zetam = z-L/2;
  double h = sqrt(4*a*rho/pow(a+rho, 2));
  double kp = sqrt(4*a*rho/(pow(a+rho, 2)+pow(zetap, 2)));
  double km = sqrt(4*a*rho/(pow(a+rho, 2)+pow(zetam, 2)));

  // Compute the Bfield
  double fact = _mu0*I/(4*M_PI*L);
  double plus = zetap*(K(kp*kp)+(a-rho)*Pi(h*h,kp*kp)/(a+rho))/sqrt(pow(a+rho, 2)+pow(zetap, 2));
  double minus = zetam*(K(km*km)+(a-rho)*Pi(h*h,km*km)/(a+rho))/sqrt(pow(a+rho, 2)+pow(zetam, 2));
  return fact*(plus-minus);
}

double FieldHandler::K(const double& m) const {

  // Define the function to be integrated
  TF1 func("func", "1./sqrt(1.-[0]*sin(x)*sin(x))", 0, M_PI/2);
  func.SetParameter(0, m);
  
  // Return its integral between 0 and pi/2. Lack of backward compatibility
  // for in ROOT 5 in ROOT 6 prevents using user defined accurary (TODO)
  return func.Integral(0, M_PI/2);
}

double FieldHandler::E(const double& m) const {

  // Define the function to be integrated
  TF1 func("func", "sqrt(1.-[0]*sin(x)*sin(x))", 0, M_PI/2);
  func.SetParameter(0, m);
  
  // Return its integral between 0 and pi/2. Lack of backward compatibility
  // for in ROOT 5 in ROOT 6 prevents using user defined accurary (TODO)
  return func.Integral(0, M_PI/2);
}

double FieldHandler::Pi(const double& n, const double& m) const {

  // Define the function to be integrated
  TF1 func("func", "1./((1-[0]*sin(x)*sin(x))*sqrt(1.-[1]*sin(x)*sin(x)))", 0, M_PI/2);
  func.SetParameter(n, m);
  
  // Return its integral between 0 and pi/2. Lack of backward compatibility
  // for in ROOT 5 in ROOT 6 prevents using user defined accurary (TODO)
  return func.Integral(0, M_PI/2);
}
