"""Simulate MD trajectory of protein in water.

Usage:
  simulate_trajectory.py [options] <input.ext>

Options :
  -h --help             Show this screen.
  --keep-water          Do not remove water molecules from PDBx/mmCIF.
  --preset=<preset>     Use pre-defined simulation settings, listed below.
  --force-field=<ff>    (preset) Force field, "amber99-implicit", "amber14-implicit", or "amber14-explicit". [default: amber14-implicit]
  --waterbox-pad=<pad>  (preset) Waterbox padding width in nm [default: 1.0].
  --temperature=<T>     (preset) System temperature in Kelvin [default: 310].
  --timestep=<ts>       (preset) Integration time step in femtoseconds [default: 1.0].
  --friction=<f>        (preset) Langevin friction in 1.0/ps [default: 0.3].
  --old-integrator      (preset) Use LangevinIntegrator, not LangevinMiddleIntegrator.
  --do-not-minimize     Do not perform energy minimization.
  --min-tol=<mintol>    Energy minimization tolerance in kJ/mol [default: 2.0].
  --log=<logsteps>      Number steps between stdout report [default: 10000].
  --sampling=<steps>    Integration time in ps [default: 100000].
  --burnin=<burnin>     Integration time in ps for burn-in steps [default: 1000].
  --spacing=<spacing>   Spacing time in ps [default: 5].
  --no-checkpointing    Do not use checkpointing.
  --cpmins=<cpmins>     Set checkpoint frequency in minutes [default: 5].
  --name=<name>         Number of simulations to run [default: ].

    
The extension of <input.ext> needs to be '.pdb' for processed PDB files that
will be directly simulated.

The following pre-defined simulation settings are available.  All presets use
a temperature of 310K, friction of 0.3/ps, and timestep of 0.5fs.  If you use
the --preset option then the following parameters, marked with "(preset)" are
determined from the preset: --force-field, --waterbox-pad, --explicit-water,
--temperature, --timestep, --friction, --old-integrator.
All other parameters are determined from their respective option values.
The available presets are as follows:

* 'amber99-implicit': AMBER99 with implicit OBC water.
* 'amber14-implicit': AMBER14 with implicit OBC1 water.
* 'amber14-explicit': AMBER14 with explicit TIP3 PFB water, 1nm waterbox.

"""

import os
import sys

import openmm as mm
import openmm.unit as u
from docopt import docopt

from md import timewarp_md
from md.npzreporter import NPZReporter, RegularSpacing

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))


def main():
    args = docopt(__doc__)
    name = args["--name"]

    input = args["<input.ext>"]
    output = input.replace(".pdb", f"{name}.npz")

    model = timewarp_md.get_openmm_model(input)
    model.addHydrogens()

    parameters = get_parameters(args)
    print(
        f"\nRunning {args['--sampling']}ps of simulation using {args['--force-field']} forcefield and integration timestep of {args['--timestep']}fs integration step and a {args['--burnin']}ps of burnin, saving frames every {args['--spacing']}ps starting from {input} saving to {output}, total of {parameters['total_steps']} steps"
    )

    simulation = setup_environment(parameters, model)
    minimize_energy(simulation, parameters["min_tol"])

    spacing = RegularSpacing(parameters["spacing"])

    print("\nBurnin ...")
    simulate(simulation, parameters["burnin_steps"], parameters["temperature"])

    protein_atoms = len(model.getPositions())

    simulation.reporters.append(
        NPZReporter(output, spacing, atom_indices=range(protein_atoms))
    )
    print("\nSampling ...")
    simulate(simulation, parameters["sampling_steps"], parameters["temperature"])


def simulate(simulation, steps, temperature):
    print("Initializing VELOCITIES to %s" % temperature)
    simulation.context.setVelocitiesToTemperature(temperature)
    if steps > 0:
        simulation.step(steps)


def minimize_energy(simulation, tolerance):
    print(f"Performing ENERGY MINIMIZATION to {tolerance}")
    simulation.minimizeEnergy(tolerance=tolerance)
    print("Completed ENERGY MINIMIZATION")


def get_parameters(args):
    burnin_steps = int(float(args["--burnin"]) * 1000 / float(args["--timestep"]))
    sampling_steps = int(float(args["--sampling"]) * 1000 / float(args["--timestep"]))
    total_steps = burnin_steps + sampling_steps
    spacing = int(float(args["--spacing"]) * 1000 / float(args["--timestep"]))

    parameters = {
        "forcefield": args["--force-field"],
        "temperature": float(args["--temperature"]) * u.kelvin,
        "friction": float(args["--friction"]) / u.picosecond,
        "timestep": float(args["--timestep"]) * u.femtosecond,
        "waterbox_pad": float(args["--waterbox-pad"]) * u.nanometers,
        "integrator": "LangevinMiddleIntegrator",
        "log_steps": int(args["--log"]),
        "min_tol": float(args["--min-tol"]),  # * u.kilojoules_per_mole,
        "total_steps": total_steps,
        "burnin_steps": burnin_steps,
        "sampling_steps": sampling_steps,
        "spacing": spacing,
    }

    return parameters


def setup_environment(parameters, model):
    simulation = timewarp_md.get_simulation_environment_from_model(model, parameters)

    if parameters["log_steps"] > 0:
        simulation.reporters.append(
            mm.app.StateDataReporter(
                sys.stdout,
                parameters["log_steps"],
                step=True,
                potentialEnergy=True,
                kineticEnergy=True,
                speed=True,
                temperature=True,
                progress=True,
                totalSteps=parameters["total_steps"],
            )
        )
    simulation.context.setPositions(model.positions)
    return simulation


if __name__ == "__main__":
    main()
