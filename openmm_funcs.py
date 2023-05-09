#
#  from openmm.app import *
from openmm import app
#  from openmm import *
import openmm
from openmm import unit


# Minimization function

def minimize_f(simulation, iters=500):
    simulation.minimizeEnergy(tolerance=0.001, maxIterations=iters)
    position = simulation.context.getState(getPositions=True).getPositions()
    outFile = "minimized_struct.pdb"
    with open(outFile, 'w') as fl:
        app.PDBFile.writeFile(simulation.topology, position, fl)
    return simulation


# Equilibration function - Constant pressure & temp

def equilibrate(
    coords: app.Topology,
    forcefield: app.ForceField,
    box_vectors: unit.quantity.Quantity,
    final_pressure: unit.Quantity = 1*unit.atmosphere,
    temp_range: range = range(0, 300, 25),
    output_state_data_filename="equilibration_state_data.csv",
    output_pdb_filename="traj_equilibration.pdb",
    friction_coeff: unit.Quantity = 1/unit.femtosecond,
    step_size: unit.Quantity = 4*unit.femtoseconds,
    time_per_temp_increment: unit.Quantity = 0.005*unit.nanoseconds,
    time_final_stage: unit.Quantity = 0.05*unit.nanoseconds,
    minimize: bool = True,
    platform="CPU"
):
    print("Initialising equilibration run...")
    # adjust the range to include the highest temp (stop value)
    inclusive_temp_range = range(
        temp_range.start,
        temp_range.stop + temp_range.step,
        temp_range.step
    )
    temperatures = unit.Quantity(inclusive_temp_range, unit.kelvin)
    steps_per_temp_increment = int(time_per_temp_increment / step_size)
    steps_final_stage = int(time_final_stage / step_size)

    # Create system
    system = forcefield.createSystem(
        coords.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1*unit.nanometer,
        constraints=app.AllBonds,
        hydrogenMass=4*unit.amu,
    )
    # Create constant temp integrator
    integrator = openmm.LangevinMiddleIntegrator(
        temperatures.min(),
        friction_coeff,
        step_size
    )
    # Create simulation and set initial positions
    simulation = app.Simulation(
        coords.topology,
        system,
        integrator,
        # Platform.getPlatformByName("OpenCL" if FORCEFIELD=="ani" else "CUDA")
        openmm.Platform.getPlatformByName(platform)
    )
    simulation.context.setPositions(coords.positions)

    # seet PBC for simulataion if there is
    if box_vectors is not None:
        simulation.context.setPeriodicBoxVectors(*box_vectors)

    simulation.reporters.append(app.PDBReporter(output_pdb_filename, 500))
    state_reporter = app.StateDataReporter(
        output_state_data_filename,
        steps_per_temp_increment//10,
        temperature=True,
        density=True,
        potentialEnergy=True,
        kineticEnergy=True,
        totalEnergy=True,
    )
    simulation.reporters.append(state_reporter)

    # Local energy minimisation
    if minimize:
        print("Local energy minimisation...")
        simulation = minimize_f(simulation, 500)

    # Heating to final temp
    print(f"Equilibrating {temperatures.min()}\
          to {temperatures.max()} in {len(temperatures)}\
          stages, {time_per_temp_increment} per stage")
    for stage, temperature in enumerate(temperatures):
        print(f"Heating stage {stage+1}/{len(temperatures)} at {temperature}")
        integrator.setTemperature(temperature)
        simulation.step(steps_per_temp_increment)
    # Final equilibration, constant pressure
    print(f"Final equilibration at {final_pressure} for {time_final_stage}")
    barostat = openmm.MonteCarloBarostat(
        final_pressure,
        temperatures.max()
    )
    system.addForce(barostat)
    simulation.step(steps_final_stage)
    print("Done")
    return simulation


# Production function - Constant pressure & volume

def production(
    coords: app.Topology,
    forcefield: app.ForceField,
    box_vectors: unit.quantity.Quantity,
    output_state_data_filename="production_state_data.csv",
    output_pdb_filename="traj_production.pdb",
    temperature: unit.Quantity = 300*unit.kelvin,
    pressure: unit.Quantity = 1*unit.atmosphere,
    friction_coeff: unit.Quantity = 1/unit.femtosecond,
    step_size: unit.Quantity = 4*unit.femtoseconds,
    duration: unit.Quantity = 1*unit.nanoseconds,
    steps_per_saved_frame: int = 500,
    platform="CPU"
):
    print("Initialising production run...")

    # Create system
    system = forcefield.createSystem(
        coords.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1*unit.nanometer,
        constraints=app.AllBonds,
        hydrogenMass=4*unit.amu,
    )
    # Create constant temp integrator
    integrator = openmm.LangevinMiddleIntegrator(
        temperature,
        friction_coeff,
        step_size
    )
    # Create simulation and set initial positions
    simulation = app.Simulation(
        coords.topology,
        system,
        integrator,
        # Platform.getPlatformByName("OpenCL" if FORCEFIELD=="ani" else "CUDA")
        openmm.Platform.getPlatformByName(platform)
    )
    simulation.context.setPositions(coords.positions)

    # seet PBC for simulataion if there is
    if box_vectors is not None:
        simulation.context.setPeriodicBoxVectors(*box_vectors)

    state_reporter = app.StateDataReporter(
        output_state_data_filename,
        steps_per_saved_frame,
        temperature=True,
        density=True,
        potentialEnergy=True,
        kineticEnergy=True,
        totalEnergy=True,
    )
    simulation.reporters.append(state_reporter)
    simulation.reporters.append(app.PDBReporter(output_pdb_filename,
                                                steps_per_saved_frame))
    #  simulation.reporters.append(app.PDBReporter('traj_production.pdb', 10))

    barostat = openmm.MonteCarloBarostat(
        pressure,
        temperature,
    )
    system.addForce(barostat)

    # Production run
    print("Running production...")
    simulation.step(int(duration / step_size))
    print("Done")
    return simulation
