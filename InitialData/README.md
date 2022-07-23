===============================================================================
Scalar Self-Force Initial Data Code, SSF_init
===============================================================================
 Version 0.1
 Copyright 2017 Niels Warburton
 
 This python code computes initial data for time-domain scalar-field self-force
 calculations. The calculation is made in the frequency domain. The homogeneous
 solutions are computed numerically by integrating the radial field equation.
 The source integration is performed using spectral source integration methods.
 The time-domain solution is reconstructed using the method of homogeneous
 solutions.
 
 Further details on these methods can be found in these papers:
 
 https://arxiv.org/abs/1103.0287
 https://arxiv.org/abs/1506.04742
 https://arxiv.org/abs/0808.2315
 
 The code works for eccentric orbits, both geodesic and accelerated
 
Requirements
------------

Currently this code requires Python version 2.x to run. The python packages it
uses are pretty standard.

Usage
----- 

SSF_init.py p e kappa n lmin lmax, where:
	p is the semi-latus rectum
	e is the orbital eccentricity
	kappa the acceleration coefficient where dt/dchi = kappa*dt/dchi_geo, i.e., set to 1 for geodesic motion
	n is the grid resolution
	lmin is the minimum l-mode to compute
	lmax is the maximum l-mode to compute
	
The code will attempt to load the coordinate data from the data/input/
subdirectory. The output will be placed in data/output.

Output format
-------------

The columns of the output data will be formatted as follows:

1. rho
2. Re[phi]
3. Im[phi]
4. Re[D[phi,r*]]
5. Im[D[phi,r*]]
6. Re[D[phi,t]]
7. Im[D[phi,t]]

License
-------
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.