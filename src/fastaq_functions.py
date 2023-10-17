def check_gc_content(sequence: str, gc_bounds: tuple or float) -> bool:
    """
    Checks if the GC content of a sequence is within the specified bounds.
    Args:
        sequence (str): The input DNA sequence.
        gc_bounds (tuple or float): Tuple with lower and upper bounds or a single float representing the upper bound.
    Returns:
        bool: True if GC content is within bounds, False otherwise.
    """
    gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence) * 100
    if isinstance(gc_bounds, tuple):
        return gc_bounds[0] <= gc_content <= gc_bounds[1]
    else:
        return gc_content <= gc_bounds


def check_length(sequence: str, length_bounds: tuple or int) -> bool:
    """
    Checks if the length of a sequence is within the specified bounds.

    Args:
        sequence (str): The input DNA sequence.
        length_bounds (tuple or int): Tuple with lower and upper bounds or a single integer representing the upper bound

    Returns:
        bool: True if length is within bounds, False otherwise.
    """
    seq_length = len(sequence)
    if isinstance(length_bounds, tuple):
        return length_bounds[0] <= seq_length <= length_bounds[1]
    else:
        return seq_length <= length_bounds


def check_quality(quality_scores: str, quality_threshold: int) -> bool:
    """
    Checks if the average quality of a sequence is above a specified threshold.

    Args:
        quality_scores (str):  The quality (Phred) score for each base, represented as an ASCII character.
        quality_threshold (int): The threshold for average quality (scale Phred33).

    Returns:
        bool: True if average quality is above threshold, False otherwise.
    """
    avg_quality = sum(ord(score) - 33 for score in quality_scores) / len(quality_scores)
    return avg_quality >= quality_threshold


def read_fastaq(input_path: str) -> dict:
    """
    Reads a FASTQ file and returns a dictionary.

    Args:
        input_path (str): The path to the FASTQ file.

    Returns:
        dict: A dictionary where keys are sequence names (starting with "@")
              and values are tuples of (sequence (line 2), quality (line 4)).

    Example:
        Given a FASTQ file like this:

        @Sequence1
        AGCTAGCTAGCTAGCT
        +
        !@#$!@#$!@#$!@#$
        @Sequence2
        CGATCGATCGATCGAT
        +
        !@#$!@#$!@#$!@#$

        The function will return:
        {'@Sequence1': ('AGCTAGCTAGCTAGCT', '!@#$!@#$!@#$!@#$'),
         '@Sequence2': ('CGATCGATCGATCGAT', '!@#$!@#$!@#$!@#$')}

    """
    seqs = {}
    with open(input_path) as fastaq:
        lines = fastaq.readlines()
        number_of_line = 0
        while number_of_line < len(lines):
            if lines[number_of_line].startswith("@"):
                name = lines[number_of_line].strip()
                sequence = lines[number_of_line + 1].strip()
                quality = lines[number_of_line + 3].strip()
                seqs[name] = (sequence, quality)
                number_of_line += 4
            else:
                number_of_line += 1
    return seqs
