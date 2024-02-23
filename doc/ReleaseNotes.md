# v02-05-02

* 2024-02-23 tmadlener ([PR#9](https://github.com/iLCSoft/KalTest/pull/9))
  - Update CI to run on latest clicdp nightlies and use the key4hep-build action

* 2022-12-02 Thomas Madlener ([PR#8](https://github.com/iLCSoft/KalTest/pull/8))
  - Migrate CI to github actions and remove obsolete travis CI configuration

# v02-05-01

* 2022-06-27 Daniel Jeans ([PR#5](https://github.com/iLCSoft/KalTest/pull/5))
  - add missing factor 0.5 to density term of Bethe-Bloch parameterisation

* 2020-04-12 Frank Gaede ([PR#4](https://github.com/iLCSoft/KalTest/pull/4))
  - fix issue w/ c++17 
        - this caused KalTest to return a “filtered” state as its “smoothed” or “inverse-filtered” state
        - patch provided by K.Fujii

# v02-05

* 2018-08-21 Bo Li ([PR#3](https://github.com/ilcsoft/KalTest/pull/3))
  - Merged Runge-Kutta (RK) track propagation functions; 
  - Added another trsansport option by using RK; 
  - Cleaned the ct_nonuniform example

# v02-04

* 2018-01-31 Frank Gaede ([PR#2](https://github.com/iLCSoft/KalTest/pull/2))
  - fix all compiler warnings (gcc54-ub1604)

* 2017-11-30 Andre Sailer ([PR#1](https://github.com/iLCSoft/KalTest/pull/1))
  - Running profiling, replace dynamic_cast with static_cast as dynamic_cast return value is never checked, see also AIDAsoft/aidaTT#19

# v02-03

# v02-03

Frank Gaede 2017-04-12 
  - Updating version to v02-03
  - Release Notes for v02-03

Marko Petric 2017-04-08 
  - Coverity integration

Marko Petric 2017-03-23 
  - Add install statement to CI
  - Add CONTRIBUTING.md and PULL_REQUEST_TEMPLATE and fix test script

Frank Gaede 2017-03-21 
  - add AUTHORS file

Marko Petric 2017-03-21 
  - Update release notes
  - Enable CI and add LICENCE

# v02-02
F.Gaede
-  made compatible with c++11 and ROOT6
- removed -ansi -pedantic -Wno-long-long
- added current source dir to include path for dictionary built

# v02-01
F.Gaede
- default measurement dimension (kSdim) set to 5
- added flag BUILD_WITH_T0_FIT (default OFF) setting the measurement dimension to 6 
- TKalDetCradle.cxx
- bug fix for intersection calculation using the newtonian solver such as cones for example (reseting the deflection angle to 0 in TKalDetCradle::Transport()) 

# v02-00
- implemented the optional use of a non-uniform B field ( Li Bo) as described in  Bo Li, Keisuke Fujii, Yuanning Gao, ("Kalman-filter-based track fitting in non-uniform magnetic field with segment-wise helical track model")[http://arxiv.org/abs/1305.7300]


# v01-05-04
- Fixed a field direction bug  (Bo Li)
- Slightly modified the track model class to make the straight track work; Straight track generation is the example ct was implemented.  (Bo Li)

# v01-05-03
- apply bug fix provided by K.Fujii that makes adjusting of the z0 for curlers in ThelicalTrack::MoveTo unnecessary

# v01-05-01
- Fixed orientation for in/out transportation in TKalDetCradle, previously only inward transportation was correct.
- Ensure that material effects are treated correctly by only including material effects from present surface to destination surface.
- THelicalTrack:MoveTo modified to make sure the helix really moves to the new reference point and is not out by more than 2PI.

# v01-05
- Added TCutCone surface with examples of use in hybrid example.
- ScatterBy function corrected in TVTrack class.

# v01-04
- geomlib 
- TCylinder IsOnSurface now checks if point is on the surface (r) not just between +/-z

# v01-03
- added name member and corresponding getName memeber function to TVMeasLayer

# v01-02
- modified to use bounded planes
- examples provided in hybrid example to demonstrate use of bounded planes in VXD and FTD 
- added extra Transport method to transport directly to a measurement layer and modified the  previous Transport to site method to use this method. No change in previous functionality, only additional functionality added
            	
# v01-01-01
- fixed problem in KalTest VERSION


# v01-01
- fixed problem from MacroRootDict.cmake: LD_LIBRARY_PATH unset with rootcint
- improved cmake files
- requires new package ILCUTIL


# v01-00
- first release as part of iLCSoft
- introduced new default units: mm, Tesla (was cm, kGauss)
- introduced cmake build files
- updated examples/kaltest/hybrid/main/EXKalTest.cxx to write track parameters d0 and tnl and errors^2 to ntuple
- Note: default units are not yet changed in example, thus the detector is scaled to 1/10 of its size with the BField 10 times larger ...

