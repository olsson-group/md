import argparse
import os
import random

import PeptideBuilder
from Bio.PDB import PDBIO
from openmm.app import Modeller, PDBFile
from pdbfixer import PDBFixer

# Observed frequencies of amino acids in vertabrates
aa_frequencies = {
    "A": 7.4,  # Alanine
    "R": 4.2,  # Arginine
    "N": 4.4,  # Asparagine
    "D": 5.9,  # Aspartic Acid
    "C": 3.3,  # Cysteine
    "E": 5.8,  # Glutamic Acid
    "Q": 3.7,  # Glutamine
    "G": 7.4,  # Glycine
    "H": 2.9,  # Histidine
    "I": 3.8,  # Isoleucine
    "L": 7.6,  # Leucine
    "K": 7.2,  # Lysine
    "M": 1.8,  # Methionine
    "F": 4.0,  # Phenylalanine
    "P": 5.0,  # Proline
    "S": 8.1,  # Serine
    "T": 6.2,  # Threonine
    "W": 1.3,  # Tryptophan
    "Y": 3.3,  # Tyrosine
    "V": 6.8,  # Valine
}
random.seed(42)

names = {
    1: "monopeptides",
    2: "dipeptides",
    3: "tripeptides",
    4: "tetrapeptides",
    5: "pentapeptides",
    6: "hexapeptides",
    7: "heptapeptides",
    8: "octapeptides",
}


def main(args):
    if args.sequence:
        sequence = args.sequence
    else:
        sequence = get_random_sequence(args.random)

    tmp_pdb_path = generate_peptide(sequence)
    sequence_length = len(sequence)
    path = f"results/{names[sequence_length]}/{sequence}"
    os.makedirs(path, exist_ok=True)
    tmp_pdb_path = fix_peptide(tmp_pdb_path)
    modeller = add_hydrogen_and_termini(tmp_pdb_path)

    output_pdb_path = os.path.join(path, f"traj.pdb")
    with open(output_pdb_path, "w") as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f)


def get_random_sequence(n):
    amino_acids = list(aa_frequencies.keys())
    weights = list(aa_frequencies.values())
    total = sum(weights)
    normalized_weights = [w / total for w in weights]

    sequence = random.choices(amino_acids, weights=normalized_weights, k=n)
    sequence = "".join(sequence)
    print(f"Generated sequence: {sequence}")
    return sequence


def generate_peptide(sequence):
    phi = [-120] * len(sequence)
    psi_im1 = [140] * len(sequence)
    structure = PeptideBuilder.make_structure(sequence, phi, psi_im1)

    io = PDBIO()
    io.set_structure(structure)
    io.save("/tmp/peptide.pdb")
    return "/tmp/peptide.pdb"


def fix_peptide(pdb):
    fixer = PDBFixer(filename=pdb)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.0)

    temp_pdb_path = f"/tmp/peptide.pdb"

    with open(temp_pdb_path, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)

    return temp_pdb_path


def add_hydrogen_and_termini(pdb_file_path):
    pdb = PDBFile(pdb_file_path)
    modeller = Modeller(pdb.topology, pdb.positions)

    return modeller


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("sequence", nargs="?")
    parser.add_argument("--random", type=int)

    main(parser.parse_args())
