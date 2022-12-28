#
from openmm import app
from openmm import unit
from openmmml import MLPotential
from openmm_funcs import (equilibrate,
                          production)


# set device
#  device = "CUDA"
device = "OpenCl"

# Load from gromac files
#  gro = app.GromacsGroFile("mol.gro")
#  forcefield = app.GromacsTopFile("mol.top",
#                                  periodicBoxVectors=gro.getPeriodicBoxVectors())

# Load from pdb
pdb = app.PDBFile("minimized_struct.pdb")
temp = 300
step_size = 2 * unit.femtoseconds
total_time = 1*unit.nanoseconds,

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
    temp_range=range(100, temp, 50),
    time_per_temp_increment=0.001*unit.nanoseconds,
    time_final_stage=0.01*unit.nanoseconds,
    step_size=step_size,
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
