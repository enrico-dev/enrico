import argparse

import numpy as np
import openmc

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--short', action='store_true')
parser.add_argument('-r', '--rings', type=int, default=5)
parser.add_argument('-a', '--axial', type=int, default=0)
args = parser.parse_args()

# Basic model parameters
fuel_or = 0.406
clad_ir = 0.414
clad_or = 0.475
pitch = 1.26
if args.short:
    fuel_length = 10.0
    boundary = 'reflective'
else:
    fuel_length = 200.0
    boundary = 'vacuum'

if args.axial == 0:
    n_axial = int(fuel_length)
else:
    n_axial = args.axial

# Create materials
uo2_density = 10.97
percent_td = 0.96

uo2 = openmc.Material(name='UO2')
uo2.add_element('U', 1.0, enrichment=4.95)
uo2.add_element('O', 2.0)
uo2.set_density('g/cm3', uo2_density*percent_td)

zirc = openmc.Material(name='Zircaloy')
zirc.add_element('Zr', 1.0)
zirc.set_density('g/cm3', 4.5)

water = openmc.Material(name='Water')
water.add_nuclide('H1', 2.0)
water.add_element('O', 1.0)
water.add_s_alpha_beta('c_H_in_H2O')

# Create cylinders
fuel_rings = []
for r in np.linspace(0., fuel_or, args.rings + 1)[1:]:
    fuel_rings.append(openmc.ZCylinder(R=r))
clad_inner = openmc.ZCylinder(R=clad_ir)
clad_outer = openmc.ZCylinder(R=clad_or)

# Division for wedges
xplane = openmc.XPlane(x0=0.0)
yplane = openmc.YPlane(y0=0.0)
xy_bounds = [
    +xplane & +yplane,
    +xplane & -yplane,
    -xplane & -yplane,
    -xplane & +yplane
]

cells = []
for wedge in xy_bounds:
    for annulus in openmc.model.subdivide(fuel_rings)[:-1]:
        c = openmc.Cell(region=annulus & wedge)
        c.fill = [uo2.clone() for i in range(n_axial)]
        cells.append(c)

    gap_annulus = +fuel_rings[-1] & -clad_inner
    gap = openmc.Cell(region=gap_annulus & wedge)
    cells.append(gap)

    clad_annulus = +clad_inner & -clad_outer
    clad = openmc.Cell(fill=zirc, region=clad_annulus & wedge)
    cells.append(clad)

    moderator = openmc.Cell(region=+clad_outer & wedge)
    moderator.fill = [water.clone() for i in range(n_axial)]
    cells.append(moderator)
slice_univ = openmc.Universe(cells=cells)

lattice = openmc.RectLattice()
lattice.lower_left = (-pitch/2, -pitch/2, 0.0)
lattice.pitch = (pitch, pitch, 1.0)
lattice.universes = np.full((n_axial, 1, 1), slice_univ)

# model boundaries
xmin = openmc.XPlane(x0=-pitch/2, boundary_type='periodic')
xmax = openmc.XPlane(x0=pitch/2, boundary_type='periodic')
ymin = openmc.YPlane(y0=-pitch/2, boundary_type='periodic')
ymax = openmc.YPlane(y0=pitch/2, boundary_type='periodic')
zmin = openmc.ZPlane(z0=0.0, boundary_type=boundary)
zmax = openmc.ZPlane(z0=fuel_length, boundary_type=boundary)
box = +xmin & -xmax & +ymin & -ymax & +zmin & -zmax

main_cell = openmc.Cell(fill=lattice, region=box)

model = openmc.model.Model()
model.geometry = openmc.Geometry([main_cell])

model.settings.particles = 1000
model.settings.inactive = 10
model.settings.batches = 50
model.settings.source = openmc.Source(
    space=openmc.stats.Box(
        (-fuel_or, -fuel_or, 0.0),
        (fuel_or, fuel_or, fuel_length),
        True
    )
)

model.export_to_xml()
