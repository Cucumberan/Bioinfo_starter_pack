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


def choose_weight(weight: str) -> dict:
    """
    Choose the weight type of amino acids - average or monoisotopic.

    Args:
        weight (str): The type of weight to choose, either 'average' or 'monoisotopic'.

    Returns:
        dict: A dictionary mapping amino acids to their weights based on the chosen type.
    """
    if weight == 'average':
        weights = {
            'A': 71.0788, 'R': 156.1875, 'N': 114.1038, 'D': 115.0886, 'C': 103.1388,
            'E': 129.1155, 'Q': 128.1307, 'G': 57.0519, 'H': 137.1411, 'I': 113.1594,
            'L': 113.1594, 'K': 128.1741, 'M': 131.1926, 'F': 147.1766, 'P': 97.1167,
            'S': 87.0782, 'T': 101.1051, 'W': 186.2132, 'Y': 163.1760, 'V': 99.1326
            }
    elif weight == 'monoisotopic':
        weights = {
            'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694, 'C': 103.00919,
            'E': 129.04259, 'Q': 128.05858, 'G': 57.02146, 'H': 137.05891, 'I': 113.08406,
            'L': 113.08406, 'K': 128.09496, 'M': 131.04049, 'F': 147.06841, 'P': 97.05276,
            'S': 87.03203, 'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841
            }
    else:
        raise ValueError(f"I do not know what '{weight}' is :( \n Read help or just do not write anything except your "
                         f"sequence")

    return weights


def aa_weight(seq: str, weight: str = 'average') -> float:
    """
    Calculate the amino acids weight in a protein sequence.

    Args:
        seq (str): The amino acid sequence to calculate the weight for.
        weight (str, optional): The type of weight to use, either 'average' or 'monoisotopic'. Default is 'average'.

    Returns:
        float: The calculated weight of the amino acid sequence.
    """
    weights_aa = choose_weight(weight)
    final_weight = 0
    for aa in seq.upper():
        final_weight += weights_aa[aa]
    return round(final_weight, 3)


def count_hydroaffinity(seq: str) -> list:
    """
    Count the quantity of hydrophobic and hydrophilic amino acids in a protein sequence.

    Args:
        seq (str): The protein sequence for which to count hydrophobic and hydrophilic amino acids.

    Returns:
        tuple: A tuple containing the count of hydrophobic and hydrophilic amino acids, respectively.
    """
    hydrophobic_count = 0
    hydrophilic_count = 0
    seq = seq.upper()

    for aa in seq:
        if aa in HYDROPHOBIC_AA:
            hydrophobic_count += 1
        elif aa in HYDROPHILIC_AA:
            hydrophilic_count += 1

    return [hydrophobic_count, hydrophilic_count]


def peptide_cutter(sequence: str, enzyme: str = "trypsin") -> str:
    """
    This function identifies cleavage sites in a given peptide sequence using a specified enzyme.
    
    Args: sequence (str): The input peptide sequence. enzyme (str): The enzyme to be used for cleavage. Choose
    between "trypsin" and "chymotrypsin". Default is "trypsin".
        
    Returns: str: A message indicating the number and positions of cleavage sites, or an error message if an invalid
    enzyme is provided.
    """
    cleavage_sites = []
    if enzyme not in ("trypsin", "chymotrypsin"):
        return "You have chosen an enzyme that is not provided. Please choose between trypsin and chymotrypsin."

    if enzyme == "trypsin":  # Trypsin cuts peptide chains mainly at the carboxyl side of the amino acids lysine or
        # arginine.
        for aa in range(len(sequence) - 1):
            if sequence[aa] in ['K', 'R', 'k', 'r'] and sequence[aa + 1] not in ['P', 'p']:
                cleavage_sites.append(aa + 1)

    if enzyme == "chymotrypsin":  # Chymotrypsin preferentially cleaves at Trp, Tyr and Phe in position P1(high
        # specificity)
        for aa in range(len(sequence) - 1):
            if sequence[aa] in ['W', 'Y', 'F', 'w', 'y', 'f'] and sequence[aa + 1] not in ['P', 'p']:
                cleavage_sites.append(aa + 1)

    if cleavage_sites:
        return f"Found {len(cleavage_sites)} {enzyme} cleavage sites at positions {', '.join(map(str, cleavage_sites))}"
    else:
        return f"No {enzyme} cleavage sites were found."



