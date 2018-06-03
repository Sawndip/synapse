![synapse](https://francois-drielsma.github.io/synapse/logo.png)

# SynAPSE

# 1. About

The Synergistic Analysis of Phase Space Evolution (**SynAPSE**)
framework provides support for a broad array of particle beam
phase space density characterization methods, both common and
novel in accelerator physics.

This package was developed with the intention of measuring the
evolution of the transverse phase space of muons across an absorber
in the context of the 
[Muon Ionization Cooling Experiment (MICE)](http://mice.iit.edu).
The experiment measures the phase space of 
individual beam muons before and after passing through the absorber 
and builds ensembles at the analysis level. Beam non linearties and
transmission losses experienced in the magnetic channel
motivated the development of non-standard phase space density 
estimation techniques. The techniques rely on the concept of
particle amplitude and 
[nonparametric density estimation](https://en.wikipedia.org/wiki/Nonparametric_statistics).
Support was added for G4Beamline ASCII output in order
to broaden the scope of the package.

This code is aimed at measuring and outputting a broad array 
of beam parameters (thoroughly described in the following),
provided with a set phase space coordinate measurements,
i.e. (x, y, px, py, pz), of a beam at a single or a set
of z positions along the beam line. The framework measures the evolution
of those parameters along the beam line and produces summary
graphs for each of them.

The package uniquely supports a wide array of non parametric
density estimators that may be used to evaluate the probability
density functions of the beams in the 2D or 4D transverse phase space. This
allows for the unbiased evaluation of beam volumes, local densities and
distributions of filamented beams in phase space.

--------------------

# 2. Installation

## 2.1. Prerequisites

The following are required to build the SynAPSE framework:
  - GCC 4.6.4 and above (C++11);
  - CMake 2.8.0 and above;
  - ROOT 5.34/36 and above.

The following are provided with the SynAPSE framework:
  - [Qhull](http://www.qhull.org/) 2015.2;
  - [NanoFLANN](https://github.com/jlblancoc/nanoflann/releases) 1.2.3.

The following is optional but required to import MAUS data:
  - [MAUS](https://launchpad.net/maus) 2.8.3 and above.

## 2.2. Build

First, if you want to import MAUS data, source the MAUS environment

    source ${MAUS_ROOT_DIR}/env.sh

If it is not included, the MAUS importers will not be built but
the rest of the code will work fine.

To build this code, simply run the all purpose builder:

    ./build.sh -j N --use-maus-gcc flag

with N the number of cores you want to build with. Set the flag to
'**true**' if you require to use the MAUS third party gcc to build
the code.

## 2.3. C++ API reference
Browse the [Doxygen documentation](https://francois-drielsma.github.io/synapse).

--------------------

# 3. Structure

## 3.1. Global variables

The default global variables used by the framework are specified in the
ConfiguratioDefaults.txt file. If the variables are changed in
the defaults, it will affect all the algorithms.

Each variable may be overridden individually at run time by adding
an option to the command line of the form

    ./program --option0 value0 --option1 value1 ... [datafile.root]

or

    ./program --option0=value0 --option1=value1 ... [datafile.root]

If the parameter cannot be found in the defaults, it is not known by the
SynAPSE framework and the exception handler will throw.

A custom set of cards may be used to overwrite the default by running

    ./program --configuration_file custom_config.txt [datafile.root]

The list of available cards and their definition may be printed by running
using the flag --help or -h.

Running the program with the option --verbose or -v will print out the
debugging information when the program is ran by rerouting the Pitch::debug
to std::out.

## 3.2. Input

Currently there are two supported data structure:
  - MAUS ROOT TTree data structure;
  - G4BL ASCII output structure.

### 3.2.1. MAUS

The MAUS output data structure may be imported in three different ways.\n

ImportMAUSData.cc is exectuted as follows:

    ./import_data [options] maus_data.root

It imports the data and the reconstructed simulation in an identical fashion,
including the cuts requested in the datacards:
  - Upstream momentum selection;
  - Track quality, aperture selection;
  - Particle species slsection.

If the Monte Carlo truth is present, it is recorded for each particle that
makes the particle selection cuts.

ImportMAUSSimulation.cc is exectuted as follows:

    ./import_sim [options] maus_sim.root

It imports the Monte Carlo truth only and applies the selection criteria
to the truth information rather than the digitized information.

ImportMAUSMinimal.cc is exectuted as follows:

    ./import_minimal [options] maus_sim.root

It imports the Monte Carlo truth only and does not apply any selection criteria.
This is useful for minimal simulations that e.g. only include virtual planes, fields
and an absorber.

### 3.2.2. G4Beamline

ImportG4BL.cc is exectuted as follows:

    ./import_g4bl [options] g4bl_ascii.txt

It imports the G4Beamline simulation truth and does not apply any selection criteria.

### 3.2.3. Output
The output of the importers is a single ROOT file structured as follows
  - One TNtuple (Data) that contains the real MICE data;
  - One TNtuple (RecMC) that contains the reconstructed Monte Carlo (digitized);
  - One TNtuple (Truth) that contains the Monte Carlo truth information corresponding to RecMC;
  - One TNtuple (UncutTruth) that contains the uncut Monte Carlo truth information;
  - Additional variables depending on the input (e.g. MAUS version).

The Truth and UncutTruth NTuples are structured as follows:

    SpillID EventID VirtualPlaneID x y z px py pz

The other two NTuples are structed as follows:

    SpillID EventID TrackerID StationID x y z px py pz xe ye ze pxe pye pze

with xe, ye, ... the reconstruction uncertainties on the phase space coordinates.

Any of the four trees may be left empty when filled and the
emittance code will proceed accordingly. Specify which data type is
requested in the datacards.

--------------------

# 4. Main scripts

All the following algorithms support an identical data structure
and are strictly ran on the imported data.

## 4.1. Phase space

PhaseSpace.cc is executed as follows

    ./phase_space [options] import.root

It reconstructs the evolution of several beam phase space summary
statistics along the beam line. The parameters include
  - Twiss parameters: &alpha;, &beta; and &gamma;;
  - mecanical angular momentum L;
  - normalised transverse emittance &epsilon;<sub>&perp;</sub>;
  - transmission;
  - mean total momentum;
  - &alpha;-amplitude A<sub>&alpha;</sub>;
  - &alpha;-subemittance e<sub>&alpha;</sub>;
  - &alpha;-fractional emittance &epsilon;<sub>&alpha;</sub>.

All the quantities are provided with a measurement uncertainty (if
the data has been reconstructed) and a statistical uncertainty. The 
quantities are computed and compiled in a single Beam class built
for each z position at which the beam is sampled.

If the '**mice**' flag is set to 1 in the datacards, the default MICE beam line
elements will be drawn on the canvases along with a information box that specifies
the MAUS version, the ISIS user cycle, the run number and the data type.

If the '**de**' flag is set to 1 in the datacards, the fractional emittance
is reconstructed using non parametric density estimation. The evolution
graphs are compiled in a single folder of the name of the input file.

The other datacards that apply to the algorithm are defined as
  - '**frac**' specifies the fraction of the beam considered in the fractional quantities;
  - '**min_vid**' ID of the first virtual plane to draw;
  - '**tku_vid**' ID of the upstream reference virtual plane;
  - '**tkd_vid**' ID of the downstream reference virtual plane;
  - '**max_vid**' ID of the last virtual plane to draw

For a definition of &alpha;-amplitude, subemittance and fractional
emittance, see https://pos.sissa.it/295/099/.

## 4.2. Amplitudes

Amplitudes.cc is executed as follows

    ./amplitudes [options] import.root

It reconstructs the amplitude distribution upstream and downstream of
the absorber. If applied to a simulation, it compares the distributions
measured at the virtual planes specified by '**tku_vid**' and '**tkd_vid**'
in the datacards.

A series of flags are available in the datcards:
  - '**corrected**' computes the amplitudes by recursively removing high amplitudes;
  - '**mcd**' computes the amplitudes based on the MCD covariance matrix;
  - '**generalised**' uses non parametric density estimation to produce amplitudes;
  - '**significance**' computes the significance of the difference between each amplitude bin;
  - '**rebin**' bins the data in buckets of equal phase space volume;
  - '**poincare**' produces Poincare sections of the phase space.

If the flag '**voronoi**' is true, the distribution of Voronoi tesselation
cell volumes is reconstructed upstream and downstream of the absorber. It
is time intensive in four dimensions.

## 4.3. Phase space profiles

Profiles.cc is executed as follows

    ./profiles [options] import.root

It reconstructs the phase space profiles upstream and downstream of
the absorber. If applied to a simulation, it reconstructs the distributions
measured at the virtual planes specified by '**tku_vid**' and '**tkd_vid**'
in the datacards.

It creates a 1D histogram for each of the phase space variable and a 2D
histogram for each possible combination of two of them and stores them to
a ROOT file.

It also produces comparison between the different data types of 1D histograms.
If the data and simulation are both provided, it produces a ratio between
simulaiton and data.

## 4.4. Density profiles

DensityProfiles.cc is executed as follows

    ./density [options] import.root

It reconstructs the nonparametric density profile upstream and downstream of
the absorber. If applied to a simulation, it reconstructs the profiles
measured at the virtual planes specified by '**tku_vid**' and '**tkd_vid**'
in the datacards.

The density profile shows the probability contour level of the density
estimator as a function of the fraction of the beam that the contour
encompasses, i.e.

![equation](https://latex.codecogs.com/gif.latex?%5Crho_%5Calpha%20%3D%20%5Crho%28F%5E%7B-1%7D%28%5Calpha%29%29%2C)

with F the cumulative density function (CDF) and &alpha; the fraction of
the beam included inside the contour. For a Gaussian beam, the profile
follows a function of the form

![equation](https://latex.codecogs.com/gif.latex?%5Crho_%5Calpha%20%3D%20%282%5Ed%5Cpi%5Ed%7C%5Cmathbf%7B%5CSigma%7D%7C%29%5E%5Cfrac12%20e%5E%7B-%5Cchi_d%5E2%28%5Calpha%29/2%7D%2C)

which scales with the size of the RMS ellipse. If the beam density has
increased, the profile is scaled up everywhere. If the beam has partially
scraped, the core increase is preserved but the tail density is reduced.

This is a good choice of variable because it is independant of the scale of
the maximum density which is prone to statistical fluctuation.

If the '**poincare**' flag is set to true, the algorithm creates poincaré
sections of the density estimate upstream and downstream of the absorber.

## 4.5. Beam reweighting routine

ReweightBeam.cc is executed as follows

    ./density [options] import_data.root import_sim.root

This algorithm computes the amplitude distributions of the source
(import_sim.root) and the target (import_data.root) and reweights the 
source to fit the target distribution. This allows to match simulations
to the observed distributions.

Two reweighting algorithm may be used:
  - Most distant bin (0): Finds the bin that is most under the target, reweight
    everything with respect to it;
  - Bin by bin (1): Scale the bins to the target from the highest to the lowest.

Both algorithms iterate until the Kolmogorov-Smirnov and Chi squared tests
show that the source and target histograms agree.

--------------------

# 5. Test scripts

The test scripts are designed to test and demonstrate the performance of
all the internal components of the SynAPSE framework. It is mostly focused
on the nonparametric density esitmators developed for the package.

## 5.1 Estimator quality

TestEstimatorQuality.cc is executed as follows

    ./test_quality [options]

This algorithm qualitatively tests the performance of a nonparametric
density estimator at reproducing the true underlying density profile of 
a known distribution. The input parameters are specified in the 
configuration defaults:
  - '**de_algo**' specifies the class of density estiator to test.

## 5.2 Estimator convergence

TestEstimatorConvergence.cc is executed as follows

    ./test_convergence [options]

This algorithm quantitatively tests the performance of a nonparametric
density estimator at converging towards the true underlying density
profile of a known distribution in terms of Mean Integrated Squared Error
(MISE). The input parameters are specified in the configuration defaults:
  - '**de_algo**' specifies the class of density estiator to test.

## 5.3 Estimator parameter optimization

TestOptimalParameters.cc is executed as follows

    ./test_optimal [options]

This algorithm uses a golden search to minimize the Mean Integrated Squared
Error (MISE) of a nonparametric density estimator. It optimizes the
free parameters of the estimators, i.e.
  - k in the k Nearest Neighbour (kNN) estimator;
  - J in the Penalised Bootstrap Aggregate Tesselation Density Estimator (PBATDE).

The input parameters are specified in the configuration defaults:
  - '**de_algo**' specifies the class of density estiator to optimize.

## 5.4 Convex hull volume convergence

TestConvexHull.cc is executed as follows

    ./test_hull [options]

This algorithm tests the bias and uncertainty on the volume reconstructed
from the convex hull of points sampled inside it. It may be used to test
the rate of convergence of the method for points uniformly sampled inside
a d-ball or inside an &alpha;-contour of a d-Gaussian distribution.

## 5.5 Alpha-complex volume convergence

TestAlphaComplex.cc is executed as follows

    ./test_alphacomplex [options]

This algorithm tests the bias and uncertainty on the volume reconstructed
from the alpha complex of points sampled inside it. It may be used to test
the rate of convergence of the method for points uniformly sampled inside
a d-ball or inside an &alpha;-contour of a d-Gaussian distribution.

## 5.6 Amplitude reconstruction

TestAmplitudeRecon.cc is executed as follows

    ./test_amplitude [options]

The algorithm tests the different modes of amplitude reconstruction and their
performance in non-linear scenarios. It also compares the different 
amplitude reconstruction scheme with each other.

A series of flags are available in the datcards:
  - '**corrected**' computes the amplitudes by recursively removing high amplitudes;
  - '**mcd**' computes the amplitudes based on the MCD covariance matrix;
  - '**generalised**' uses non parametric density estimation to produce amplitudes.


--------------------

# 6. Toy Monte Carlo

The toy Monte Carlo code allows for the fast production and transport of ideal beams.
ToySimulation.cc is executed as follows

    ./toy_sim [options]

## 6.1 Generation

The generator produces a Gaussian beam based on the cards specified in the 
ConfigurationsDefaults.txt or in the command line arguments:
  - '**toy_mass**' specifies the mass of the beam particles in MeV/c^2;
  - '**toy_n**' specifies the amount of particles to be generated;
  - '**toy_seed**' sets the seem of the random number generator;
  - '**toy_eps**' sets the normalised emittance of the beam in mm;
  - '**toy_mom**' sets the longitudinal momentum of the particles in the beam in MeV/c;
  - '**toy_alpha**' sets the beam &alpha; function at production;
  - '**toy_beta**' sets the beam &beta; function at production in mm.

The beam is generated just upstream of the toy absorber.

## 6.1 Abosrber

The absorber characteristics are specified in the data cards. The energy loss is
computed in a purely deterministic fashion by integrating the Bethe-Bloch formula
over the full length of the absorber (EnergyLoss):

![equation](https://latex.codecogs.com/gif.latex?K%5Crho%5Cfrac%7BZ%7D%7BA%7D%5Cfrac1%7B%5Cbeta%5E2%7D%5Cleft%28%5Cln%5Cleft%28%5Cfrac%7B2m_ec%5E2%5Cbeta%5E2%5Cgamma%5E2%7D%7BI%7D%5Cright%29-%5Cbeta%5E2%20%5Cright%20%29.)

The scattering is a stochastic process which is approximatively described
by a Gaussian distribution of RMS angle:

![equation](https://latex.codecogs.com/gif.latex?%5Ctheta_0%20%3D%20%5Cfrac%7B13.6%5Ctext%7BMeV%7D%7D%7Bp%5Cbeta%20c%7D%5Csqrt%7B%5Cfrac%7BL%7D%7BX_0%7D%7D%5Cleft%5B1&plus;0.038%5Cln%5Cleft%28L/X_0%20%5Cright%20%29%20%5Cright%20%5D%2C)

with X0 the radiation length and L and the length of the absorber.

The physics processes are controlled through:
  - '**toy_scat**' turns the scattering on or off (0/1);
  - '**toy_eloss**' turns the energy loss on or off (0/1).

## 6.2. Transport

Beam line elements may be added to the toy simulation right downstream of the
toy absorber. One may add:
 - '**toy_drift**' adds a drift space of required length in mm;
 - '**toy_sol**' adds a solenoid of required length in mm.

The beam line elements are reprensted by unitary transfer matrices that rotate 
the phase space vectors of the beam particles in a linear fashion through

![equation](https://latex.codecogs.com/gif.latex?%5Cmathbf%7Bx%7D%27%20%3D%20%5Cmathbf%7BM%7D%5Cmathbf%7Bx%7D%2C) 

which in turn transforms the covariance matrix throuh

![equation](https://latex.codecogs.com/gif.latex?%5Cmathbf%7B%5CSigma%7D%27%20%3D%20%5Cmathbf%7BM%7D%5Cmathbf%7B%5CSigma%7D%5Cmathbf%7BM%7D%5ET.)

## 6.3 Output
The toy simulation algorithm outputs the beam at the initial reference plane (0)
and after all the beam line elements (1). The output file may be used in the main
scripts by specifying the **tku_vid** and **tkd_vid**.

--------------------

# 7. Backend packages

## 7.1 Exception handler

Propriatary exception handling and output stream
  - Produces an accurate description of any exceptions (Exceptions);
  - Redirects the output stream of different classes of messages (Pitch).

## 7.2 Mathematics package

The mathematics package includes
  - Linear algebra (Matrix, Vector, DecompLU, DecompQR);
  - Computational geometry package (Geometry.hh);
  - Statistics package (Statistics.hh);
  - Broad set of n-dimensional distribution types (Gaus, Exponential, Uniform, etc.)

## 7.3. Density estimation package

Density estimation package
  - Broad array of non parametric density estimators:
    + Delaunary Tesselation Field Estimator (DTFE);
    + Kernel Density Esimator (KDE);
    + k-Neareast Neighbour density estimator (KNearestNeighbours);
    + Local Reachability density estimator (LocalReachability);
    + Optimal binning density estimator (OptimalBinning);
    + Tesselation Density Estimator (TDE);
    + Penalised Bootstrap Aggregate Tesselation Denisty Estimator (PBATDE).
  - Support for computational geometric methods to support the estimators
    + Delaunay triangulation (Delaunay);
    + Voronoi tesselation (Voronoi);
    + Alpha-complicies (AlphaComplex).
  - Interpolators to fill the gaps between discrete estimates
    + n-linear interpolator (LinearInterpolator);
    + simplical interpolator (SimplexInterpolator).

 
--------------------

# 8. Latest version updates

Version **0.6.0**
  - Optimized the main algorithms (Amplitudes.cc, PhaseSpace.cc, Profiles.cc, DenistyProfiles.cc)
  - Created a new namespace (Beam) that encompasses every object that deal with beams
  - Created the Stream object, an array of beam bunches at different z positions (Bunch)
  - Created a drawing assistant (Drawer) to replace the old plotting tools (removed)
  - Created a data extractor (Extractor)
  - Does not compute the measurement uncertainties by default anymore (Bunch)
  - Fixed the spiral function, now fully supported (DSpiral)
  - Added verbosity option for the progress bar (ProgressBar)
  - Many bug fixes

Version **0.5.2**
  - Optimized reweighting algorithm (ReweightBeam.cc)
  - Added tool to draw data-sim comparisons (PlotTools.hh)
  - Bug fixes

Version **0.5.1**:

  - Added the default data cards (may be overridden)
  - Added support for custom data cards
  - Added support for command line arguments
  - Bug fixes

Version **0.5.0**:

  - First Doxygen documentation
  - Beam data importers
  - Beam data structure and methods
  - Toy monte carlo package
  - Propriatary exception handling and output stream
  - Mathematics package
  - Density estimation package
  - Global variables handler
  - Geometry handler
  - Beam line aperture handler
  - Beam plotting tools

