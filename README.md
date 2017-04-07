# KalTest
[![Build Status](https://travis-ci.org/iLCSoft/KalTest.svg?branch=master)](https://travis-ci.org/iLCSoft/KalTest)
[![Coverity Scan Build Status](https://scan.coverity.com/projects/12354/badge.svg)](https://scan.coverity.com/projects/ilcsoft-kaltest)

KalTest is distributed under the [GPLv3 License](http://www.gnu.org/licenses/gpl-3.0.en.html)

[![License](https://www.gnu.org/graphics/gplv3-127x51.png)](https://www.gnu.org/licenses/gpl-3.0.en.html)


## Organization of the KalTest package

- bin: contains shell scripts to configure your system
- conf:contains an Imake template
- doc:documments
- examples/kaltest: sample programs to illustrate usage of the libraries
- src/geomlib:geometry liabrary containing classes to describe track models and measurement layers
- kallib: general purpose abstract base classes to implement Kalman filter algorithm (TVKalSystem, TVKalSite, TVKalState)
- kaltracklib: derived classes that implement pure virtuals of kallib for track fitting purpose and abstract base classes to describe
a. individual measurement layers (TVMeasLayer)
b. detector system as a collection of measurement layers carrying information on the material budget (TVKalDetector)
c. hit point on each measurement layer (TVTrackHit)
- utils: a set of utility classes extracted from LEDA
- Makefile: makefile to build the libraries
- README: this file
- setup: script to setup environmental variables to run the sample programs

## How to Build the Kalman Filter Libraries
- Install ROOT if you haven't.

### How to Run Test Programs
- cd to the examples/kaltest subdirectory.
- Choose one of the subsubdirectories in it:
cdc: track fitting example for a jet-chamber-like central tracker
ct:  track fitting example for a simple cylindrical tracker
hybrid: hybrid track fitting example for VTX+IT(Barell/Fwd/Bwd)+TPC
simple: a simple line fit example using the base kallib libraries only
- cd into the subsubdirectory you chose, build the test program, and run
as follows:

For hybrid:
```
make
cd main/prod
./EXKalTest [-b] [# events] [pt in GeV] [t0 in ns]
//Notice that without -b, you will get into event display mode.
```
For others:
```
xmkmf -a
make
cd prod
./EXKalTest [# events] [pt in GeV]
```


## License and Copyright
Copyright (C), KalTest Authors

KalTest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License long with this program.  If not, see <http://www.gnu.org/licenses/>.
