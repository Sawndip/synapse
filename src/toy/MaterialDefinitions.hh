// C++ includes
#include <string>

  /** @brief Data structure that stores the element properties.
   *
   * 	   Needs to be provided with the PDG characteristics.
   */
  struct Element {
    std::string	name;	///< Name of the element/molecule/crystal
    double	ZoA;	///< Ratio of the atomic number over the atomtic mass Z/A [mol/g]
    double	rho;	///< Density [g/cm^3]
    double	I;	///< Mean excitation potential [eV]
    double 	a;	///< Prefactor of the power low that describes the density effects
    double 	k;	///< Exponent of the power low that describes the density effects
    double 	X0;	///< Value of log10(bg) below which the density effects are negligible
    double 	X1;	///< Value above which the density effects attain their asymptotic value
    double	C;	///< \f$2\ln(I/h\nu_p)+1\f$
    double 	delta0;	///< Relevant to the density effects in conductors
  };

  /** @brief Data structure that stores the material properties */
  struct Material {
    std::string 		name;		///< Name of the material
    size_t	 		number;		///< Number of elements
    std::vector<Element> 	elements;	///< List of materials that makes it up
    std::vector<double>		fractions;	///< Mass fraction represented by each eleement
  };

  // Absorbers
  static Element MiceLiH    =  {"MiceLiH",	// Name
			       	 0.56716,	// Z/A [mol/g]
			   	 0.693,		// rho [g/cm^3]
			   	 36.5,		// I [eV]
			   	 0.90567,	// a 
			   	 2.5849,	// k
			   	 -0.0988,	// X0
			      	 1.4515,	// X1
			   	 2.3580,	// C
			   	 0.};		// delta0

  static Material AbsLiH    =  	{"AbsLiH",	// Name
				 1,		// Number
			       	 {MiceLiH},	// Elements
			   	 {1.}};		// Fractions

  // Time-of-flight hodoscopes material
  static Element PVT = 		{"PVT",		// Name
			       	 0.54141,	// Z/A [mol/g]
			   	 1.032,		// rho [g/cm^3]
			   	 64.7,		// I [eV]
			   	 0.16101,	// a 
			   	 3.2393,	// k
			   	 0.1464,	// X0
			      	 2.4855,	// X1
			   	 3.1997,	// C
			   	 0.};		// delta0

  static Material TofPVT    =  	{"TofPVT",	// Name
				 1,		// Number
			       	 {PVT},		// Elements
			   	 {1.}};		// Fractions

  // Tracker fibres material
  static Element PS = 		{"PS",		// Name
			       	 0.53768,	// Z/A [mol/g]
			   	 1.060,		// rho [g/cm^3]
			   	 68.7,		// I [eV]
			   	 0.16454,	// a 
			   	 3.2224,	// k
			   	 0.1647,	// X0
			      	 2.5031,	// X1
			   	 3.2999,	// C
			   	 0.};		// delta0

  static Material TrackerPS =  	{"TrackerPS",	// Name
				 1,		// Number
			       	 {PS},		// Elements
			   	 {1.}};		// Fractions

  // Tracker He window
  static Element Al = 		{"Aluminium",	// Name
			       	 0.48181,	// Z/A (Z = 13, A = 26.9815385) [mol/g]
			   	 2.699,		// rho [g/cm^3]
			   	 166.0,		// I [eV]
			   	 0.08024,	// a 
			   	 3.6345,	// k
			   	 0.1708,	// X0
			      	 3.0127,	// X1
			   	 4.2395,	// C
			   	 0.12};		// delta0

  static Material WindowAl =  	{"WindowAl",	// Name
				 1,		// Number
			       	 {Al},		// Elements
			   	 {1.}};		// Fractions

  // Diffuser materials
  static Element Ni    =  	{"Nickel",	// Name
			       	 0.47706,	// Z/A (Z = 28, A = 58.6934) [mol/g]
			   	 8.902,		// rho [g/cm^3]
			   	 311.0,		// I [eV]
			   	 0.16496,	// a 
			   	 2.8430,	// k
			   	 -0.0566,	// X0
			      	 3.1851,	// X1
			   	 4.3115,	// C
			   	 0.10};		// delta0

  static Element Cu =  		{"Copper",	// Name
			       	 0.45636,	// Z/A (Z = 29, A = 63.546) [mol/g]
			   	 8.960,		// rho [g/cm^3]
			   	 322.0,		// I [eV]
			   	 0.14339,	// a 
			   	 2.9044,	// k
			   	 -0.0254,	// X0
			      	 3.2792,	// X1
			   	 4.4190,	// C
			   	 0.08};		// delta0

  static Element Zn =  		{"Zinc",	// Name
			       	 0.45886,	// Z/A (Z = 30, A = 65.38) [mol/g]
			   	 7.133,		// rho [g/cm^3]
			   	 330.0,		// I [eV]
			   	 0.14714,	// a 
			   	 2.8652,	// k
			   	 0.0049,	// X0
			      	 3.3668,	// X1
			   	 4.6906,	// C
			   	 0.08};		// delta0

  static Element W =  		{"Tungsten",	// Name
			       	 0.40252,	// Z/A (Z = 74, A = 183.84) [mol/g]
			   	 19.300,	// rho [g/cm^3]
			   	 727.0,		// I [eV]
			   	 0.15509,	// a 
			   	 2.8447,	// k
			   	 0.2167,	// X0
			      	 3.4960,	// X1
			   	 5.4059,	// C
			   	 0.14};		// delta0

  static Material Brass    =  	{"Brass",		// Name
				 2,			// Number
			       	 {Cu, Zn},		// Elements
			   	 {.65, .35}};		// Fractions

  static Material TAM1000    =  {"TAM1000",		// Name
				 3,			// Number
			       	 {W, Ni, Cu},		// Elements
			   	 {.9, .075, 0.025}};	// Fractions
