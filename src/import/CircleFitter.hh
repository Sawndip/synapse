#ifndef ROOTHEADERS_HPP_INCLUDED
#define ROOTHEADERS_HPP_INCLUDED
#ifdef __GNUC__
#pragma GCC system_header
#endif
#include "Math/Functor.h"
#include "Fit/Fitter.h"
#endif

/** @file Fitter.hh
    @brief Contains the definition of the circle fitter.
 */

/** @brief Fits a circle to a set of points.
 *
 *  @param	x	Set of x coordinates
 *  @param	y	Set of y coordinates
 *  @param	x0	x coordinate of the circle centre
 *  @param	y0	y coordinate of the circle centre
 *  @param	rad	Raidus of the circle
 */
bool FitCircle(const std::vector<double>& x,
	       const std::vector<double>& y,
               double& x0, double& y0, double& rad);
