def check_valid_sequence(sequences: tuple) -> bool:
    """
    Checks if the input sequences consist only of valid DNA or RNA characters.

    Args:
        sequences (iterable of str): List of sequences to be validated.
    Returns:
        bool: True if all sequences are valid, False otherwise.
    """
    allowed_characters = set('ATGCUatgcu')
    for sequence in sequences:
        for nucleotide in sequence:
            if nucleotide not in allowed_characters:
                return False
    return True
