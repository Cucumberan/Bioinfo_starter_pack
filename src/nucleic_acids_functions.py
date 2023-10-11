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


def is_dna(sequences: tuple) -> bool:
    """
    Checks if all sequences in the input list consist only of valid DNA characters.

    Args:
        sequences (iterable of str): List of sequences to be validated.

    Returns:
        bool: True if all sequences are valid DNA, False otherwise.
    """
    allowed_characters = set('ATGCatgc')
    for sequence in sequences:
        for nucleotide in sequence:
            if nucleotide not in allowed_characters:
                return False
    return True


def transcribe(sequences: tuple) -> str or list[str]:
    """
    Transcribes DNA sequences to RNA sequences.

    Args:
        sequences (iterable of str): List of DNA sequences to be transcribed.

    Returns:
        str or list of str: Transcribed RNA sequence(s).
    """
    for sequence in sequences:
        if not is_dna(sequence):
            raise ValueError("At least one of your sequences is RNA instead of DNA, and RNA can not be transcribed")
    first_sequence = get_first_sequence(sequences)
    if first_sequence:
        return first_sequence.replace("T", "U").replace('t', 'u')
    else:
        return [sequence.replace("T", "U").replace('t', 'u') for sequence in sequences]

