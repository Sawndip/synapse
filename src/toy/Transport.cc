#include "Transport.hh"

Transport::Transport() :
  _dim(0), _M(0, 0) {
}

Transport::Transport(const size_t dim) :
  _dim(dim), _M(dim, dim) {

  // Flag if the dimension requested is not supported
  if ( _dim != 2 && _dim != 4 )
      throw(Exceptions::Exception(Exceptions::nonRecoverable,
	    "Only supports 2 and 4 dimension transverse samples",
	    "Transport::Transport"));

  // Initialize the transfer matrix to the unit matrix
  size_t i;
  for (i = 0; i < _dim; i++)
        _M[i][i] = 1.;
}

Transport::Transport(const Transport& transport) {
  *this = transport;
}

Transport& Transport::operator=(const Transport& transport) {
  if ( this == &transport )
      return *this;

  _dim = transport._dim;
  _M = transport._M;

  return *this;
}

Transport::~Transport () {}

void Transport::TransportParticle(Matrix<double>& particle) {

  particle = _M*particle;
}

void Transport::TransportBunch(std::map<std::string, std::vector<double>>& beam) {

  // Fill an array of the variables to consider depending on the dimension
  std::vector<std::string> vars;
  if ( _dim == 2 )
      vars = {"x", "px"};
  if ( _dim == 4 )
      vars = {"x", "px", "y", "py"};

  // Transport the particles in the beam
  Matrix<double> particle(_dim, 1);
  size_t i, j;
  for (i = 0; i < beam["pz"].size(); i++) {
    // Fill the column matrix that corresponds to the particle
    for (j = 0; j < _dim; j++)
	particle[j][0] = beam[vars[j]][i];

    // Transport the particle
    TransportParticle(particle);

    // Fill the beam with the transported particle
    for (j = 0; j < _dim; j++)
	beam[vars[j]][i] = particle[j][0];
  }
}

void Transport::AddDriftSpace(const double& L, const double& p) {

  if ( _dim == 2 )
      _M *= Matrix<double>({{1, 	L/p},
			    {0, 	1}});

  if ( _dim == 4 )
      _M *= Matrix<double>({{1, 	L/p,	0, 	0},
			    {0, 	1,	0,	0},
			    {0,		0,	1,	L/p},
			    {0,		0,	0,	1}});
}

void Transport::AddSolenoid(const double& L, const double& K) {

  double C(cos(K*L)), S(sin(K*L));
  if ( _dim == 2 )
      _M *= Matrix<double>({{C, 	S/K},
			    {-K*S, 	C}});

  if ( _dim == 4 )
/*      _M *= Matrix<double>({{C, 	S/K,	0, 	0},
			    {-K*S, 	C,	0,	0},
			    {0,		0,	C,	S/K},
			    {0,		0,	-K*S,	C}});*/
      _M *= Matrix<double>({{C*C, 	S*C/K,	S*C, 	S*S/K},
			    {-K*S*C, 	C*C,	-K*S*S,	S*C},
			    {-S*C,	-S*S/K,	C*C,	S*C/K},
			    {K*S*S,	-S*C,	-K*S*C,	C*C}});
}
