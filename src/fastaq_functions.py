import os


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


def check_quality(quality_scores, quality_threshold: int) -> bool:
    """
    Checks the average quality of a sequence, accepting both preprocessed numerical quality scores
    and raw ASCII character quality scores.

    This function allows for flexible handling of quality scores, whether they come directly from FASTQ files
    as ASCII characters or have been preprocessed into numerical scores. It calculates the average quality
    and compares it to a specified threshold to determine if the sequence meets the quality criteria.

    Args:
        quality_scores: Numerical list of quality scores or a string of ASCII quality characters.
        quality_threshold (int): The threshold for average quality.

    Returns:
        bool: True if the average quality is above the threshold, False otherwise.

    Raises:
        ValueError: If `quality_scores` is neither a string nor a list/tuple.
    """
    # If quality_scores is a string, assume these are ASCII characters, raw data from a FASTQ file
    if isinstance(quality_scores, str):
        avg_quality = sum(ord(score) - 33 for score in quality_scores) / len(quality_scores)
    # If quality_scores is a list or tuple, assume these are numerical quality scores
    elif isinstance(quality_scores, (list, tuple)):
        avg_quality = sum(quality_scores) / len(quality_scores)
    else:
        raise ValueError("quality_scores must be either a string or a list/tuple")

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


def save_fastq_from_dict(filtered_seqs: dict, output_filename=None) -> None:
    """
    Save sequences from a dictionary to a FASTQ file.

    Args:
        filtered_seqs (dict): A dictionary where keys are sequence names and values
                              are tuples of (sequence, quality).
        output_filename (str, optional): The output filename. If not provided,
                                         the default is 'fastq_filtrator_results/filtered_data.fastq'.

    Returns:
        None

    """
    if not output_filename:
        output_filename = 'fastq_filtrator_results/filtered_data.fastq'
    else:
        output_filename = f'fastq_filtrator_results/{output_filename}.fastq'

    output_folder = os.path.dirname(output_filename)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    with open(output_filename, 'w') as output_file:
        for name, (sequence, quality) in filtered_seqs.items():
            output_file.write(f'{name}\n{sequence}\n+\n{quality}\n')
