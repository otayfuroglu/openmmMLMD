#
from openmm import app
from openmm import unit
from openmmml import MLPotential
from openmm_funcs import (equilibrate,
                          production)


# set device
device = "CPU"
#  device = "CUDA"
#  device = "OpenCL"

# Load from gromac files
#  gro = app.GromacsGroFile("mol.gro")
#  forcefield = app.GromacsTopFile("mol.top",
#                                  periodicBoxVectors=gro.getPeriodicBoxVectors())

# Load from pdb
#  pdb_path = "solv.pdb"
pdb_path = "/truba_scratch/otayfuroglu/openmmMLMD/test/akocak_test/test.pdb"
pdb = app.PDBFile(pdb_path)
box_vectors = pdb.topology.getPeriodicBoxVectors()

temp = 300
step_size = 2 * unit.femtoseconds
total_time = 0.5*unit.nanoseconds

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
    box_vectors,
    temp_range=range(300, temp, 100), # initial temp, ultimate temp, temp step size for increasing
    time_per_temp_increment=0.001*unit.nanoseconds,
    time_final_stage=0.01*unit.nanoseconds,
    step_size=step_size,
    minimize=True,
    platform=device,
)

# Production
simulation.positions = simulation.context\
        .getState(getPositions=True).getPositions()
production(
    simulation,
    forcefield,
    box_vectors,
    temperature=temp*unit.kelvin,
    step_size=step_size,
    duration=total_time,
    output_pdb_filename="traj_production.pdb",
    output_state_data_filename="production_state_data.csv",
    platform=device,
)
