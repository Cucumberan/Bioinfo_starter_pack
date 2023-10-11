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


def contains_T_and_U_at_the_same_time(sequences: tuple) -> bool:
    """
    Checks if any sequence in the input list contains both 'T' and 'U' nucleotides simultaneously.
    Args:
        sequences (iterable of str): List of sequences to be checked.
    Returns:
        bool: True if any sequence contains both 'T' and 'U', False otherwise.
    """
    for sequence in sequences:
        sequence = sequence.upper()
        if sequence.count("T") and sequences.count("U"):
            return True
    return False


def get_first_sequence(my_tuple: tuple or list[str]) -> str:
    """
    Extracts the first sequence from the input tuple, if applicable.

    Args:
        my_tuple (str or tuple): Input tuple of sequences.

    Returns:
        str or None: The first sequence if available, otherwise None.

    """
    if isinstance(my_tuple, str):
        return my_tuple
    elif len(my_tuple) == 1:
        return str(my_tuple[0])
