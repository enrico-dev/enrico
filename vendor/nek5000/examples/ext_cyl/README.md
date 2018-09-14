# Flow past a cylinder in 2D.

A unit-diameter cylinder is centered at (0,0) in a box on [-15,35] x [-15,15].

With unit inflow velocity at x=-15, the Reynolds number is given by

  Re = 1/viscosity

For the current case, we have Re=100, which can be specified in the
.rea file either by setting parameter 2, p2=.01 (the viscosity) or 
p2=-100.  If Nek5000 sees a negative value for parameter 2 it sets
viscosity to the reciprocal, viscosity = 1/|p2|.

## MESH:

The mesh was built by combining a file box.rea generated with genbox, 
using the input file fpcyl.box, with a second file, import.rea, which
contains the cylinder definition.  This file was generated with prenek.

If you wish to change the farfield resolution, you can edit fpcyl.box,
rerun genbox, and then merge box.rea with import.rea.

This case has a relatively fine spectral element mesh and should give
reasonable results for lx1=6 (polynomial order N=5).

The `mkmesh` script shows how to build the box and combine it with
the existing circle mesh.   Just type 

`mkmesh`

to execute the script.

## OUTPUTS:

 * A simple ad-hoc Strouhal number estimator has been added to the .usr file 
 * The drag and lift coefficients are also computed by calls in userchk.
 * Torque can be added by modifying the iftout flag in the torq_calc call.
 * The history of the field values at the probes located at points defined in ext_cyl.his is calculated by call to hpts in userchk. These values are appended to ext_cyl.his file.



