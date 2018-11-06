import argparse
from math import pi

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

m5_niobium = 0.01    # http://publications.jrc.ec.europa.eu/repository/bitstream/JRC100644/lcna28366enn.pdf
m5_oxygen = 0.00135  # http://publications.jrc.ec.europa.eu/repository/bitstream/JRC100644/lcna28366enn.pdf
m5_density = 6.494   # 10.1039/C5DT03403E
m5 = openmc.Material(name='M5')
m5.add_element('Zr', 1.0 - m5_niobium - m5_oxygen)
m5.add_element('Nb', m5_niobium)
m5.add_element('O', m5_oxygen)
m5.set_density('g/cm3', m5_density)

# NuScale DCA, Ch. 4, Table 4.1-1
psia = 0.0068947572931683625  # MPa
system_pressure = 1850*psia
core_avg_temperature = (543 - 32)*5/9 + 273.15  # K
water_density = openmc.data.water_density(core_avg_temperature, system_pressure)
water = openmc.Material(name='Water')
water.add_nuclide('H1', 2.0)
water.add_element('O', 1.0)
water.add_s_alpha_beta('c_H_in_H2O')
water.set_density('g/cm3', water_density)

# Create cylinders
radii = np.linspace(0., fuel_or, args.rings + 1)
fuel_rings = [openmc.ZCylinder(R=r) for r in radii[1:]]
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
    for i, annulus in enumerate(openmc.model.subdivide(fuel_rings)[:-1]):
        c = openmc.Cell(region=annulus & wedge)
        c.fill = [uo2.clone() for i in range(n_axial)]
        # Set volume for each material
        volume = pi/4*(radii[i+1]**2 - radii[i]**2)
        for m in c.fill:
            m.volume = volume
        cells.append(c)

    gap_annulus = +fuel_rings[-1] & -clad_inner
    gap = openmc.Cell(region=gap_annulus & wedge)
    cells.append(gap)

    clad_annulus = +clad_inner & -clad_outer
    clad = openmc.Cell(fill=m5, region=clad_annulus & wedge)
    clad.fill.volume = pi/4*(clad_or**2 - clad_ir**2)
    cells.append(clad)

    moderator = openmc.Cell(region=+clad_outer & wedge)
    moderator.fill = [water.clone() for i in range(n_axial)]
    for m in moderator.fill:
        m.volume = (pitch**2 - pi*clad_or**2)/4
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
model.settings.temperature = {
    'default': 523.15,
    'method': 'interpolation',
    'range': (300.0, 1500.0),
    'multipole': True
}

model.export_to_xml()
