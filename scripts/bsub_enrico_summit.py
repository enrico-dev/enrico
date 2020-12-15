#!/usr/bin/env python3

# Based on vendor/nekRS/scripts/nrsqsub_summit
import os
import shutil
import subprocess
import datetime
from configparser import ConfigParser
from xml.etree import ElementTree

class BsubSummitDriver:

    def __init__(self, n_nodes, time, proj_id=None, backend="CUDA"):
        self.casename = None
        self.n_nodes = int(n_nodes)
        self.time = time
        self.proj_id = proj_id
        self.backend = backend.upper()

        self.occa_cache_dir = os.getcwd() + "/.cache/occa"

        # Specific to Summit
        self.nvme_home = "/mnt/bb/" + os.environ["USER"]
        self.xl_home = "/sw/summit/xl/16.1.1-3/xlC/16.1.1"
        self.gpus_per_node = 6
        self.cores_per_socket = 21

        if not self.proj_id:
            self.proj_id = os.environ["PROJ_ID"]

        if self.backend == "CUDA":
            self.n_tasks = self.n_nodes * self.gpus_per_node
            self.device_number = "0"
        elif self.backend == "SERIAL":
            self.n_tasks = self.n_nodes * self.cores_per_socket * 2
            self.device_number = "LOCAL-RANK"
        else:
            raise ValueError("Incorrect value for NekRS backend.  Allowed values are: CUDA, SERIAL")

        self.setup_env()
        self.setup_files()

    def setup_env(self):
        self.nekrs_home = os.environ.get('NEKRS_HOME')
        if not self.nekrs_home:
            raise RuntimeError("NEKRS_HOME was not defined in environment")
        if not os.path.isdir(self.nekrs_home):
            raise NotADirectoryError("NEKRS_HOME was set to {self.nekrs_home}, but that directory doesn't exist")
        nekrs_env = { "OCCA_CACHE_DIR" : self.occa_cache_dir, 
                "NEKRS_HYPRE_NUM_THREADS" : "1", 
                "OGS_MPI_SUPPORT" : "1", 
                "OCCA_CXX" : f"{self.xl_home}/bin/xlc", 
                "OCCA_CXXFLAGS" : "-O3 -qarch=pwr9 -qhot -DUSE_OCCA_MEM_BYTE_ALIGN=64", 
                "OCCA_LDFLAGS" : f"{self.xl_home}/lib/libibmc++.a" }
        os.environ.update(nekrs_env)

    def setup_files(self):
        # Get info from enrico.xml
        enrico_xml = ElementTree.parse("enrico.xml")
        root = enrico_xml.getroot()
        if root.find('heat_fluids').find('driver').text != 'nekrs':
            raise ValueError("enrico.xml was not configured to use nekrs as the heat_fluids driver")
        self.casename = root.find('heat_fluids').find('casename').text

        files = [f"{self.casename}.{ext}" for ext in ["par", "co2", "udf", "oudf", "re2"]]
        for f in files:
            if not os.path.isfile(f):
                raise FileNotFoundError(f"Could not find {f}, needed for NekRS")

        os.makedirs(self.occa_cache_dir, exist_ok=True)

        # Write correct info in parfile
        par = ConfigParser()
        parfile = self.casename + ".par"
        par.read(parfile)
        rewrite = False
        if par.get('OCCA', 'backend', fallback='') != self.backend:
            par.set('OCCA', 'backend', self.backend)
            rewrite = True
        if par.get('OCCA', 'devicenumber', fallback='') != self.device_number:
            par.set('OCCA', 'devicenumber', self.device_number)
            rewrite = True
        if rewrite:
            backup_file = "{}.{}.bk".format(parfile, datetime.datetime.now().isoformat().replace(':', '.'))
            print(f"Rewriting {parfile} with necessary runtime params.  Original will be moved to {backup_file}")
            shutil.copy2(parfile, backup_file)
            with open(parfile, "w") as f:
                par.write(f)

    def precompile(self):
        cmd = f"mpirun -pami_noib -np 1 {self.nekrs_home}/bin/nekrs --setup {self.casename} --build-only {self.n_tasks} --backend {self.backend}"
        print("Running precompile command:", cmd)
        # Timeout is 2 hours.  This won't be run in bsub, so I'm using an explicit timeout
        subprocess.run(cmd, shell=True, check=True, timeout=7200)

    def bsub(self):
        enrico_bin = f"{self.nekrs_home}/bin/enrico"
        if self.backend == "SERIAL":
            jsrun_cmd = f"jsrun -X 1 -n{self.n_nodes} -r1 -a1 -c1 -g0 -b packed:1 -d packed cp -a {self.occa_cache_dir}/* {self.nvme_home}; export OCCA_CACHE_DIR={self.nvme_home}; jsrun -X 1 -n{self.n_tasks} -a{self.cores_per_socket} -c{self.cores_per_socket} -g0 -b packed:1 -d packed {enrico_bin}" 
        else: # backend == "CUDA"
            jsrun_cmd = f"jsrun -X 1 -n{self.n_nodes} -r1 -a1 -c1 -g0 -b packed:1 -d packed cp -a {self.occa_cache_dir}/* {self.nvme_home}; export OCCA_CACHE_DIR={self.nvme_home}; jsrun --smpiargs='-gpu' -X 1 -n{self.n_tasks} -r{self.gpus_per_node} -a1 -c2 -g1 -b rs -d packed {enrico_bin}" 
        bsub_cmd = f'bsub -nnodes {self.n_nodes} -alloc_flags NVME -W {self.time} -P {self.proj_id} -J nekRS_{self.casename} "{jsrun_cmd}"'

        print("Running bsub command:", bsub_cmd)
        subprocess.run(bsub_cmd, shell=True, check=True)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--nodes", type=int, required=True, help="Specify the number of NODES (not processes) to use")
    parser.add_argument("-t", "--time", type=str, required=True, help="Sets the runtime limit of the job")
    parser.add_argument("-b", "--backend", type=str, help="Sets the OCCA kernel backend [default: CUDA]", choices=['CUDA', 'SERIAL'], default='CUDA')
    parser.add_argument("--no-precompile", action="store_true", help="Skips pre-compile step for NekRS kernels")
    args = parser.parse_args()

    driver = BsubSummitDriver(n_nodes=args.nodes, time=args.time, backend=args.backend)
    if not args.no_precompile:
        driver.precompile()
    driver.bsub()

