from src.fastaq_functions import check_gc_content, check_length, check_quality


def filter_fastq(seqs: dict, gc_bounds: tuple = (0, 100), length_bounds: tuple = (0, 2 ** 32),
                 quality_threshold: int = 0) -> dict:
    """
    Filters a dictionary of FASTQ sequences based on specified criteria.

    Args:
        1. seqs (dict): Dictionary of sequences. Keys are sequence names, values are tuples of (0) sequence and (1)
        quality scores.
        2. gc_bounds (tuple or float, optional): Tuple with lower and upper bounds or a single float
        representing the upper bound for GC content (default is (0, 100)).
        3. length_bounds (tuple or int, optional): Tuple with lower and upper bounds or a single integer representing
        the upper bound for sequence length (default is (0, 2**32)).
        4. quality_threshold (int, optional): The threshold for average quality (default is 0).

    Returns:
        dict: Filtered dictionary of sequences.
    """
    filtered_seqs = {}

    for seq_name, (sequence, quality) in seqs.items():
        if not check_gc_content(sequence, gc_bounds):
            continue

        if not check_length(sequence, length_bounds):
            continue

        if not check_quality(quality, quality_threshold):
            continue

        filtered_seqs[seq_name] = (sequence, quality)

    return filtered_seqs



