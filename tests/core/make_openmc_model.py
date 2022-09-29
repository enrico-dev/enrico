import argparse
from math import pi

import numpy as np
import openmc

parser = argparse.ArgumentParser()
parser.add_argument('-ax', '--assemx', type=int, default=3)
parser.add_argument('-ay', '--assemy', type=int, default=3)
parser.add_argument('-px', '--pinx', type=int, default=2)
parser.add_argument('-py', '--piny', type=int, default=2)
parser.add_argument('-s', '--skip', action='store_true') # skip corners
args = parser.parse_args()

# parse input
na_x = args.assemx
na_y = args.assemy
np_x = args.pinx
np_y = args.piny

# Basic model parameters
fuel_or = 0.406
clad_ir = 0.414
clad_or = 0.475
pitch = 1.26
fuel_length = 10.0

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

# outer universes
water_univ = openmc.Universe(cells=[openmc.Cell(fill=water)])
void_univ = openmc.Universe(cells=[openmc.Cell(name='void')],
                                    name='null')

# Create cylinders
fuel_outer = openmc.ZCylinder(r=fuel_or)
clad_inner = openmc.ZCylinder(r=clad_ir)
clad_outer = openmc.ZCylinder(r=clad_or)

# make pin
pin = openmc.model.pin([fuel_outer, clad_inner, clad_outer],
                       [uo2, void_univ, m5, water])

# set material volumes
for i, c in pin.cells.items():
    m = c.fill
    r = []
    for j, s in c.region.get_surfaces().items():
        r.append(s.r)
    if min(r) == clad_or:
        # moderator
        vol = pitch**2 - pi*min(r)**2
    elif max(r) == fuel_or:
        vol = pi*max(r)**2
    else:
        vol = pi*(max(r)**2 - min(r)**2)
    m.volume = vol

# make assembly lattice
assem_pitch_x = pitch * np_x
assem_pitch_y = pitch * np_y
assem_lat = openmc.RectLattice()
assem_lat.lower_left = (-assem_pitch_x/2., -assem_pitch_y/2.,)
assem_lat.pitch = [pitch, pitch]
assem_lat.universes = np.full((np_y, np_x), pin)
assem_lat.outer = water_univ

# make assembly
assembly_region = openmc.rectangular_prism(assem_pitch_x, assem_pitch_y,
    origin=(0, 0, fuel_length/2.0))
assembly_cell = openmc.Cell(fill=assem_lat, region=assembly_region)
assembly = openmc.Universe(cells=[assembly_cell])

# make core lattice
core_pitch_x = assem_pitch_x * na_x
core_pitch_y = assem_pitch_y * na_y
core_lat = openmc.RectLattice()
core_lat.lower_left = (-core_pitch_x/2., -core_pitch_y/2., 0)
core_lat.pitch = [assem_pitch_x, assem_pitch_y, fuel_length]
core_lat.universes = np.full((1, na_y, na_x), assembly)
core_lat.outer = void_univ

# replace corners with assem outer
if args.skip:
    core_lat.universes[0][0][0] = water_univ
    core_lat.universes[0][na_y - 1][0] = water_univ
    core_lat.universes[0][na_y - 1][na_x - 1] = water_univ
    core_lat.universes[0][0][na_x - 1] = water_univ

# make core
core_region = openmc.rectangular_prism(core_pitch_x, core_pitch_y,
    origin=(0, 0, fuel_length/2.0), boundary_type='vacuum')
core = openmc.Cell(fill=core_lat, region=core_region)

model = openmc.model.Model()
model.geometry = openmc.Geometry([core])

# simulation settings
model.settings.particles = 1000
model.settings.inactive = 10
model.settings.batches = 50
model.settings.source = openmc.Source(
    space=openmc.stats.Box(
        (-core_pitch_x/2., -core_pitch_y/2., 0.0),
        (core_pitch_x/2., core_pitch_y/2., fuel_length),
        True)
)
model.settings.temperature = {
    'default': 523.15,
    'method': 'interpolation',
    'range': (300.0, 1500.0),
    'multipole': True
}

model.export_to_xml()
