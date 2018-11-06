#!/usr/bin/env python3

import glob
from multiprocessing import Pool
import os
from pathlib import Path
import tarfile
import tempfile
from urllib.parse import urljoin
import warnings
import zipfile

import openmc.data
from openmc._utils import download


def process_neutron(filename, output_dir):
    """Process ENDF neutron sublibrary file into HDF5 and write into a
    specified output directory."""
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            data = openmc.data.IncidentNeutron.from_njoy(
                filename, temperatures=temperatures
            )
    except Exception as e:
        print(filename, e)
        raise
    data.export_to_hdf5(output_dir / f'{data.name}.h5', 'w')
    print(f'Finished {filename}')


base_endf = 'http://www.nndc.bnl.gov/endf/b7.1/zips/'
files = [
    (base_endf, 'ENDF-B-VII.1-neutrons.zip', 'e5d7f441fc4c92893322c24d1725e29c'),
    (base_endf, 'ENDF-B-VII.1-photoat.zip', '5192f94e61f0b385cf536f448ffab4a4'),
    (base_endf, 'ENDF-B-VII.1-atomic_relax.zip', 'fddb6035e7f2b6931e51a58fc754bd10'),
    (base_endf, 'ENDF-B-VII.1-thermal_scatt.zip', 'fe590109dde63b2ec5dc228c7b8cab02')
]
wmp_version = '1.1'
wmp_base = f'https://github.com/mit-crpg/WMP_Library/releases/download/v{wmp_version}/'
wmp_filename = f'WMP_Library_v{wmp_version}.tar.gz'

temperatures = [293.6, 500.0, 750.0, 1000.0, 1250.0]
pwd = Path.cwd()
output_dir = pwd / 'endf71_multitemp'

os.makedirs(output_dir / 'photon', exist_ok=True)

with tempfile.TemporaryDirectory() as tmpdir:
    # Save current working directory and temporarily change dir
    os.chdir(tmpdir)
    library = openmc.data.DataLibrary()

    # =========================================================================
    # Download files from NNDC server
    for base, fname, checksum in files:
        download(urljoin(base, fname), checksum)

    # =========================================================================
    # EXTRACT FROM ZIP FILES

    for _, f, _ in files:
        print(f'Extracting {f}...')
        zipfile.ZipFile(f).extractall()

    # =========================================================================
    # PROCESS INCIDENT NEUTRON DATA IN PARALLEL

    with Pool() as pool:
        neutron_files = sorted(glob.glob('neutrons/*.endf'))
        results = []
        for f in neutron_files:
            r = pool.apply_async(process_neutron, (f, output_dir))
            results.append(r)
        for r in results:
            r.wait()

    for f in sorted(glob.glob(f'{output_dir}/*.h5')):
        library.register_file(f)

    # =========================================================================
    # PROCESS THERMAL SCATTERING FOR WATER

    print('Processing H2O...')
    h2o = openmc.data.ThermalScattering.from_njoy(
        'neutrons/n-001_H_001.endf',
        'thermal_scatt/tsl-HinH2O.endf',
        stdout=True
    )
    filename = f'{output_dir}/c_H_in_H2O.h5'
    h2o.export_to_hdf5(filename, 'w')
    library.register_file(filename)
    print('Finished processing H2O')

    # =========================================================================
    # INCIDENT PHOTON DATA

    for z in range(1, 101):
        element = openmc.data.ATOMIC_SYMBOL[z]
        print('Generating HDF5 file for Z={} ({})...'.format(z, element))

        # Generate instance of IncidentPhoton
        photo_file = Path('photoat') / f'photoat-{z:03}_{element}_000.endf'
        atom_file = Path('atomic_relax') / f'atom-{z:03}_{element}_000.endf'
        data = openmc.data.IncidentPhoton.from_endf(photo_file, atom_file)

        # Write HDF5 file and register it
        outfile = output_dir / 'photon' / f'{element}.h5'
        data.export_to_hdf5(outfile, 'w')
        library.register_file(outfile)

    # =========================================================================
    # WINDOWED MULTIPOLE DATA

    # Download and extract data
    download(urljoin(wmp_base, wmp_filename))
    with tarfile.open(wmp_filename, 'r') as tgz:
        tgz.extractall(output_dir)
    os.rename(output_dir / 'WMP_Library', output_dir / 'wmp')

    # Add multipole data to library
    for f in sorted(glob.glob(f'{output_dir}/wmp/*.h5')):
        library.register_file(f)

    library.export_to_xml(output_dir / 'cross_sections.xml')

    # Change back to original directory
    os.chdir(pwd)
