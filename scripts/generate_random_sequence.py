import argparse
import random

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


def main(args):
    for _ in range(args.n):
        sequence = get_sequence(args.l)
        print(sequence)


def get_sequence(n):
    amino_acids = list(aa_frequencies.keys())
    weights = list(aa_frequencies.values())
    total = sum(weights)
    normalized_weights = [w / total for w in weights]

    sequence = random.choices(amino_acids, weights=normalized_weights, k=n)
    return "".join(sequence)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("l", type=int, help="Length of the sequence")
    parser.add_argument("n", type=int, help="Number of sequences to generate")
    main(parser.parse_args())
