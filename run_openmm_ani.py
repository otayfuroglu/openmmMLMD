#
import os
import pandas as pd
import csv
from openmm import app
from openmm import unit
from openmmml import MLPotential
from openmm_funcs import (equilibrate,
                          production)


# set device
<<<<<<< HEAD
device = "CPU"
#  device = "CUDA"
#  device = "OpenCL"
=======
#  device = "CUDA"
device = "OpenCl"
>>>>>>> 70d9c7386982845255101dfc2fc9a5a5012687c0

# Load from gromac files
#  gro = app.GromacsGroFile("mol.gro")
#  forcefield = app.GromacsTopFile("mol.top",
#                                  periodicBoxVectors=gro.getPeriodicBoxVectors())
#file1 = open("./list.txt", 'r')
#file2 = file1.readlines()
#for mole in file2:
#    print(mole)
#    mole=mole.replace("\n","")
#    print(mole)

    # Load from pdb
pdb_path = "/truba_scratch/otayfuroglu/openmmMLMD/test/mudu_works/test_1mole/mobley_5631798.pdb"
pdb = app.PDBFile(pdb_path)
temp = 300
step_size = 2 * unit.femtoseconds
total_time = 1*unit.nanoseconds

# set ani potential as a forcefield
forcefield = MLPotential('ani2x')

# Load pdb into modeller and add solvent
modeller = app.Modeller(pdb.topology, pdb.positions)
# modeller.addExtraParticles(forcefield)
#  modeller.addHydrogens(forcefield)
#  modeller.addSolvent(forcefield, model='tip3p',
#                      padding=1*nanometer, neutralize=False)

# Equilibrate
simulation = equilibrate(
    modeller,
    forcefield,
    temp_range=range(100, temp+1, 100),
    time_per_temp_increment=0.001*unit.nanoseconds,
    time_final_stage=0.01*unit.nanoseconds,
    step_size=step_size,
    minimize=False,
    platform=device,
)


# Production
simulation.positions = simulation.context\
        .getState(getPositions=True).getPositions()
production(
    simulation,
    forcefield,
    temperature=temp*unit.kelvin,
    step_size=step_size,
    duration=total_time,
    output_pdb_filename="traj_production.pdb",
    output_state_data_filename="production_state_data.csv",
    platform=device,
)
