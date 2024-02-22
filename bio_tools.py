from src.fastaq_functions import check_gc_content, check_length, check_quality, read_fastaq, save_fastq_from_dict
from src.nucleic_acids_functions import check_valid_sequence, contains_T_and_U_at_the_same_time, transcribe, \
    complement, reverse_complement, reverse
from src.protein_functions import is_aa, aa_weight, count_hydroaffinity, peptide_cutter, one_to_three_letter_code, \
    sulphur_containing_aa_counter

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.Seq import Seq


def filter_fastq(input_path: str, output_filename: str = None, gc_bounds: tuple = (0, 100), length_bounds: tuple = (0, 2 ** 32),
                 quality_threshold: int = 0) -> dict:
    """
        Filters a dictionary of FASTQ sequences based on specified criteria. Saves the output FASTAQ file.

        Args:
            input_path (str): Path to the input FASTQ file.
            output_filename (str, optional): Name of the output FASTQ file. If not provided,
                the filtered sequences will not be saved to a file (default is None).
            gc_bounds (tuple or float, optional): Tuple with lower and upper bounds or a single float
                representing the upper bound for GC content (default is (0, 100)).
            length_bounds (tuple or int, optional): Tuple with lower and upper bounds or a single integer
                representing the upper bound for sequence length (default is (0, 2**32)).
            quality_threshold (int, optional): The threshold for average quality (default is 0).

        Returns:
            dict: Filtered dictionary of sequences. Keys are sequence names, values are tuples of
            (0) sequence and (1) quality scores.
        """
    seqs = read_fastaq(input_path)
    filtered_seqs = {}

    for seq_name, (sequence, quality) in seqs.items():
        if not check_gc_content(sequence, gc_bounds):
            continue

        if not check_length(sequence, length_bounds):
            continue

        if not check_quality(quality, quality_threshold):
            continue

        filtered_seqs[seq_name] = (sequence, quality)

    save_fastq_from_dict(filtered_seqs, output_filename)

    return filtered_seqs


def filter_fastq_with_Bio(input_path: str, output_filename: str = None, gc_bounds: tuple = (0, 100), length_bounds: tuple = (0, 2 ** 32), quality_threshold: int = 0) -> dict:
    """
    Filters FASTQ sequences from fastq format file based on specified criteria. Saves the output FASTAQ file. Uses Biopython libraries
    """
    filtered_seqs = {}
    
    for record in SeqIO.parse(input_path, "fastq"):
        sequence = str(record.seq)
        quality_scores = record.letter_annotations["phred_quality"]
        gc_content =  gc_fraction(record.seq)*100

        if not (gc_bounds[0] <= gc_content <= gc_bounds[1]):
            continue
        
        if not (length_bounds[0] <= len(sequence) <= length_bounds[1]):
            continue
        
        if not check_quality(quality_scores, quality_threshold):
            continue
        
        filtered_seqs[record.id] = (sequence, quality_scores)
    
    if output_filename:
        with open(output_filename, "w") as output_handle:
            SeqIO.write((SeqIO.SeqRecord(Seq(seq), id=seq_id, description="", letter_annotations={"phred_quality": quality}) for seq_id, (seq, quality) in filtered_seqs.items()), output_handle, "fastq")
    
    return filtered_seqs


def run_dna_rna_tools(*arguments):
    """
    Executes DNA/RNA sequence manipulation procedures.

    Args:
        *arguments (tuple): Variable-length argument list containing sequences and procedure.

    Returns:
        str or list of str: Result of the selected procedure.
    """
    procedure = arguments[-1]
    sequences = arguments[:-1]
    if not check_valid_sequence(sequences):
        raise ValueError("At least one of your sequences does not correspond to either DNA or RNA")
    if contains_T_and_U_at_the_same_time(sequences):
        raise ValueError(
            "One of your sequences contains both thymine and uracil at the same time, which is not possible((((")
    if procedure == "transcribe":
        return transcribe(sequences)
    elif procedure == "reverse":
        return reverse(sequences)
    elif procedure == "complement":
        return complement(sequences)
    elif procedure == "reverse_complement":
        return reverse_complement(sequences)
    else:
        return "Something went wrong, please, verify the chosen procedure is written correctly"


def run_amino_analyzer(sequence: str, procedure: str, *, weight_type: str = 'average', enzyme: str = 'trypsin'):
    """
    This is the main function to run the amino-analyzer.py tool.
    
    Args:
        sequence (str): The input protein sequence in one-letter code.
        procedure (str): amino-analyzer.py tool has 5 functions at all:
            1. aa_weight - Calculate the amino acids weight in a protein sequence. Return float weight
            weight_type = 'average': default argument for 'aa_weight' function. weight_type = 'monoisotopic' can be
            used as a second option.
            2. count_hydroaffinity - Count the quantity of hydrophobic and hydrophilic amino acids in a protein
            sequence. Return list in order: hydrophobic, hydrophilic
            3. peptide_cutter - This function identifies cleavage sites in a given peptide sequence using a specified
            enzyme. Return list of cleavage sites enzyme = 'trypsin': default argument for 'peptide_cutter' function.
            enzyme = 'chymotrypsin' can be used as a second option.
            4. one_to_three_letter_code - This function converts a protein sequence from one-letter amino acid code
            to three-letter code. Return string of amino acids in three-letter code
            5. sulphur_containing_aa_counter - This function counts sulphur-containing amino acids in a protein
            sequence. Return quantity of sulphur-containing amino acids.

    Returns:
        The result of the specified procedure.

    Raises:
        ValueError: If the procedure is not recognized or if the input sequence contains non-amino acid characters.

    Note: - Supported amino acid characters: V, I, L, E, Q, D, N, H, W, F, Y, R, K, S, T, M, A, G, P, C, v, i, l, e,
    q, d, n, h, w, f, y, r, k, s, t, m, a, g, p, c. - Make sure to provide a valid procedure name and sequence for
    analysis. :param enzyme: :param sequence: :param procedure: :param weight_type:
    """
    procedures = ['aa_weight', 'count_hydroaffinity', 'peptide_cutter', 'one_to_three_letter_code',
                  'sulphur_containing_aa_counter']
    if procedure not in procedures:
        raise ValueError(f"Incorrect procedure. Acceptable procedures: {', '.join(procedures)}")

    if not is_aa(sequence):
        raise ValueError("Incorrect sequence. Only amino acids are allowed (V, I, L, E, Q, D, N, H, W, F, Y, R, K, S, "
                         "T, M, A, G, P, C, v, i, l, e, q, d, n, h, w, f, y, r, k, s, t, m, a, g, p, c).")
    result = ''
    if procedure == 'aa_weight':
        result = aa_weight(sequence, weight_type)
    elif procedure == 'count_hydroaffinity':
        result = count_hydroaffinity(sequence)
    elif procedure == 'peptide_cutter':
        result = peptide_cutter(sequence, enzyme)
    elif procedure == 'one_to_three_letter_code':
        result = one_to_three_letter_code(sequence)
    elif procedure == 'sulphur_containing_aa_counter':
        result = sulphur_containing_aa_counter(sequence)
    return result



