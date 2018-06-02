// C++ includes
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <chrono>

// ROOT includes
#include "TF2.h"
#include "TRandom3.h"

// Other includes
#include "DGaus.hh"
#include "Matrix.hh"

/** @file ToyTools.hh
    @brief Tools used by the Toy Monte Carlo algorithm
 */

/** @brief Generates a Gaussian bunch from the given beam parameters
 *
 *  @param	m	Particle mass [MeV/c^2]
 *  @param	p	Total beam momentum [MeV/c]
 *  @param	n	Number of particles to generate
 *  @param	eps	Transverse normalised emittance [mm]
 *  @param	beta	Transverse beta function [mm]
 *  @param	alpha	Transverse alpha function (correlation)
 *
 *  @return		Input beam: x, (y), px, (py), pz [mm, MeV/c]
 */
std::map<std::string, std::vector<double>> GaussianBunch(const double& m,
							 const double& p,
							 const double& n,
							 const double& eps,
							 const double& beta,
							 const double alpha=0.);
