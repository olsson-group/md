This repo is mostly lifted from https://github.com/microsoft/timewarp 

I recommend setting up a local conda environment, python virtual environment doesn't work because we need packages from conda-forge

```bash
conda create --prefix ./venv python=3.9
conda activate ./venv
```

Now run the install script, it will install the environment in the venv.

Now we can start a MD simulating using the following command:

```bash
python scripts/simulate.py --timestep {timestep in fs} --sampling {sampling time in ps} --burnin {burnin time in ps} --spacing {spacing betewen samples in ps} path/to/pdbfile.pdb
```

This will save a trajectory at `path/to/pdbfile.npz` with positions, velocities, and forces at every {spacing} ps.

it is also possible to generate pdb files for amino acid sequences using 

```bash
python scripts/generate_peptide.py {sequence}
```

It is also possible to generate random sequences using amino acids in the same frequency as in vertabrates with with length n


```bash
python scripts/generate_peptide --random n
```




# Environment Setup
To get started, I recommend using a local conda environment. A Python virtual environment won't work as we need some packages from `conda-forge`.
Set up and activate the environment with the following commands:

``` bash
conda create --prefix ./venv python=3.9
conda activate ./venv
```

Next, run the provided install script. This will install all necessary dependencies into the venv environment.

``` bash
source install.sh
```

# Running Molecular Dynamics Simulations
You can start a molecular dynamics simulation using the following command:

``` bash
python scripts/simulate.py --timestep {timestep in fs} --sampling {sampling time in ps} --burnin {burnin time in ps} --spacing {spacing between samples in ps} path/to/pdbfile.pdb
```

--timestep: The integration time step in femtoseconds (fs).
--sampling: The sampling time in picoseconds (ps).
--burnin: The burn-in time in picoseconds (ps).
--spacing: The spacing between samples in picoseconds (ps).

This command will save a trajectory file at path/to/pdbfile.npz containing:

Positions
Velocities
Forces
The trajectory will have data at every {spacing} ps.

# Generating PDB Files for peptides
From Amino Acid Sequences
To generate PDB files for specific amino acid sequences, use the following command:

``` bash
python scripts/generate_peptide.py {sequence}
```
For Random Sequences

To generate random amino acid sequence frequencies and a specific length n, use:

``` bash
python scripts/generate_peptide.py --random n
```

Both of these scripts will save peptides to `results/{mono/di/tri/...}/{sequence}/traj.pdb`.

