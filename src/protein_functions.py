AA_SET = {'V', 'I', 'L', 'E', 'Q', 'D', 'N', 'H', 'W', 'F', 'Y', 'R', 'K', 'S', 'T', 'M', 'A', 'G', 'P', 'C',
          'v', 'i', 'l', 'e', 'q', 'd', 'n', 'h', 'w', 'f', 'y', 'r', 'k', 's', 't', 'm', 'a', 'g', 'p', 'c'}
HYDROPHOBIC_AA = ['A', 'V', 'L', 'I', 'P', 'F', 'W', 'M']
HYDROPHILIC_AA = ['R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'K', 'S', 'T', 'Y']
AMINO_ACIDS = {'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
               'K': 'Lys', 'L': 'Leu', 'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg', 'S': 'Ser',
               'T': 'Thr', 'V': 'Val',
               'W': 'Trp', 'Y': 'Tyr'}


def is_aa(seq: str) -> bool:
    """
    Check if a sequence contains only amino acids.

        Args:
        seq (str): The input sequence to be checked.

    Returns:
        bool: True if the sequence contains only amino acids, False otherwise.
    """
    unique_chars = set(seq)
    return unique_chars <= AA_SET
