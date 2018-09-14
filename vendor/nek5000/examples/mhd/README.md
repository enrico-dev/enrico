# Helical Forcing Dynamo.

This example gives a cicular Polarized Flow of Galloway & Proctor [1]. 
Magnetic field growth rate (or 1/2 of magnetic energy growth rate)
is equal to ~ 0.3  for K_z=0.57 at Rm=100 (magnetic Reynolds #)

Note: magnetic field growth rate (or 1/2 of magnetic energy growth
rate) is equal to ~ 0.26  for K_z=0.57 at Rm=10 (CKPT95)

Files:
        gpf.*		Base case (Re=Rm=10)
        gpf_re.rea	Case  w/ Re=Rm=100
        gpf_b.*		Base case for GTP solver w/ .box file
        gpf_m.*		Base case for GTP solver w/ .map file

## MHD:

For an MHD case, B-field is stored after temperature array
(and passive scalars if any) so ldimt>=2.  Note that in all routines
w/ NEKUSE common block i.e. uservp, userf, userq, userbc & useric,
B-field can be accessed through ux,uy,uz when
they are called with ifield=ifldmhd (ifield=1 is for velocity,
ifield=2 for temperature or ifield=ifldmhd=2 when IFHEAT=.FALSE.,
etc.) in addition to the direct reference bx(ix,iy,iz,ieg),
by(ix,iy,iz,ieg), etc.
Also  note that a call to outpost in userchk enforces B-field dump
in odd numbers of .f/.fld files in place of velocity. Then for
restart, one has to provide two files in .rea file -- first,
velocity .f/.fld (even number) file name and on the next line, a
file name of the B-field dump (odd .f/.fld number) both preceeded
with a line

      2 PRESOLVE/RESTART OPTIONS

A major limitation of the current version of MHD implementation,
is the requirement that the boundary conditions for magnetic field 
have to be of the same type as for velocity  field, e.g.,
both periodic ('PER'), Direchlet ('v  '), etc.
Also ioinfo control and avg_all routine do not have a support 
for magnetic field yet.

## GTP (Global Tensor Product Solver) :

param(116)=p116, p117 & p118 can be used to activate global
tensor product solver based on global fast-diagonalization method
(for box-like geometries only) which alternatively, can also be
switched on by ndim<0 in .rea file on the line w/ "NEL,NDIM,NELV".
In the latter (recommended) mode of ndim<0, the mesh w/ bcs are read
from .box file, nel should be positive but no mesh w/ bcs should be
in .rea file, no .map file is needed and p116 -- p118 are ignored.
In alternative mode, one may still use .map file (and .rea/.re2 mesh w/ bcs)
with gtp/gfdm solver by setting p116, p117 & p118 to nelx, nely & nelz,
correspondingly, with an option of being multiplied by a minus sign
that indicates a primary coordinate direction for element partitioning.
If p116 -- p118 are positive or ndim<0 then the default partitioning
for gtp/gfdm solver is in nely*nelz pencils of length nelx so nely*nelz
should be greater than a number of processors np or for greater parallel
effeciency, to be a multiple of np.
Note that to use gtp/gfdm solver one also needs lelx, lely, lelz
in SIZE file to be equal or greater than nelx, nely, nelz,
correspondingly, in addition to uncommenting the second pair of
lines w/ lelg_sm & ltfdm2 in nek5_svn/trunk/nek/ZPER (make sure
you use 'make clean' first)
The solver enforces divergence-free condition upto machine accuracy
and is especilly good for meshes w/ elements of high aspect ratio
but the solver is not necessary for MHD problems with Nek5000
so one can use standard Nek5000 solver w/ .map file by keeping
zero p116 -- p118 and positive ndim.

[1] Galloway & Proctor, Nature, 356, pp 692--693, 1992;
Cattaneo, Kim, Proctor & Tao, Phys. Rev. Lett., 75(8),
pp. 1523--1525, 1995 (CKPT95)