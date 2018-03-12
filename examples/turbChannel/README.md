# LES of a turbulent channel flow.

HiOCFD5 Testcase: WS2 / Resolution Coarse (ID:MS3)
LES of a turbulent channel flow with Re_tau=550

Reference: 
[1] https://how5.cenaero.be/content/ws2-les-plane-channel-ret550 

[2] Myoungkyu Lee and Robert D. Moser,  
Direct numerical simulation of turbulent channel flow up to 
Re_tau = 5200, 2015, Journal of Fluid Mechanics, vol. 774, 
pp. 395-415

Note: 1-D statistics are dumped to 
mean_prof.dat and vel_fluc_prof.dat 

Two SGS models are available:
- high pass filtering (set filtering=hpfrt in .par file)
  Typically it's a good idea to use a filterCutOffRation = 0.7
  (2 modes for e.g. N=7) using a filter weight of 40.

- dynamic Smagorinsky (set filtering=none in .par file) 
  This version computes C_s D^2 (mixing length) dynamically and thus 
  does not need an explicit definition of the filter width D. 
  It does, however, need a definition of the grid- to test-filter
  width ratio.  See comp_mij().)

  If preferred, one can also explicitly define a filter width
  and compute C_s using contractions involving this ad-hoc
  filter width, then compute C_s D^2 by multiplying by D^2.
  In this code, simply comment out the early "return" in 
  the routine "set_grid_spacing()".  Preliminary tests indicate 
  no discernable difference between the two approaches and the
  former (i.e., current configuration) should thus be preferred.

  C_s is computed using the Lilly contraction, plus planar and
  temporal averaging of the numerator and denominator.

  Here, the only clipping of the eddy viscosity is simply to
  ensure that C_s > 0.  That is, C_s = max(C_s,0).
