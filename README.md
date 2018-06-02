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
density functions of the beams in 2D, 4D or 5D phase space. This
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

# 4. Main algorithms

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

## 4.4. Density profiles

## 4.5. Beam reweighting routine

--------------------

# 5. Beam

--------------------

# 6. Toy Monte Carlo

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

