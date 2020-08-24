#!/usr/bin/env python

import argparse
import xml.etree.ElementTree as ET

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import openmc
import openmc.lib
from openmc.data import JOULE_PER_EV
from tqdm import tqdm

mpl.rcParams['axes.labelsize'] = 16
mpl.rcParams['axes.xmargin'] = 0
mpl.rcParams['axes.ymargin'] = 0
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = ['Times']
mpl.rcParams['font.size'] = 14
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['pdf.use14corefonts'] = True
mpl.rcParams['savefig.bbox'] = 'tight'
mpl.rcParams['text.usetex'] = True

# Determine reactor power
tree = ET.parse('enrico.xml')
root = tree.getroot()
power = root.find('coupling').find('power')
Q = float(power.text)

pitch = 1.26

parser = argparse.ArgumentParser()
parser.add_argument('-r', '--resolution', type=int, default=500)
parser.add_argument('filename')
parser.add_argument('axial', type=float)
args = parser.parse_args()


with openmc.StatePoint(args.filename, autolink=False) as sp:
    # Get list of material bins and mean values
    t = sp.tallies[1]
    cell_inst = t.filters[0].bins
    heat = JOULE_PER_EV * t.mean.flatten()
    total_heat = np.sum(heat)

    openmc.lib.init()

    index_cell_inst = {tuple(b): i for i, b in enumerate(cell_inst)}

    N = args.resolution
    z = args.axial
    d = pitch/N
    xs = np.linspace(-pitch/2 + d/2, pitch/2 - d/2, N)
    ys = np.linspace(pitch/2 - d/2, -pitch/2 + d/2, N)
    img = np.full((N, N), np.nan)
    for row, y in enumerate(tqdm(ys)):
        for col, x in enumerate(xs):
            cell, inst = openmc.lib.find_cell((x, y, z))

            b = (cell.id, inst)

            idx = index_cell_inst.get(b)
            if idx is not None:
                if heat[idx] > 0:
                    img[row, col] = heat[idx]

    plt.imshow(img, cmap='RdYlBu_r', extent=[-pitch/2, pitch/2, -pitch/2, pitch/2])
    plt.xlabel('x [cm]')
    plt.ylabel('y [cm]')
    cb = plt.colorbar()
    cb.set_label(r'Heat generation rate [W/cm$^3$]')
    plt.savefig('openmc_nek_heat_.pdf')
    openmc.lib.finalize()
