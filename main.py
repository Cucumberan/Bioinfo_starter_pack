#HW_14
from abc import ABC, abstractmethod
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.Seq import Seq

class BiologicalSequence(ABC):
    """
    Abstract base class representing a biological sequence.

    Attributes:
        sequence (str): The biological sequence as a string.

    Methods:
        __len__(): Returns the length of the sequence.
        is_valid(): Abstract method to check if the sequence is valid based on the specified alphabet.
        __getitem__(key): Allows for indexing the sequence.
        __str__(): Returns the string representation of the sequence.
        __repr__(): Returns a string representation that can be used to recreate the object (official representation).
    """

    def __init__(self, sequence: str):
        self.sequence = sequence
 
    def __len__(self) -> int:
        return len(self._sequence)

    @abstractmethod
    def is_valid(self) -> bool:
        """
        Checks if the sequence matches the specified alphabet.
        Returns:
            bool: True if the sequence is correct, otherwise False.
        """
        pass

    def __getitem__(self, key):
        return self.sequence[key]
    
    def __str__(self) -> str:
        return self.sequence
    
    def __repr__(self) -> str:
        return f"{self.__class__.__name__}('{self.sequence}')"



class NucleicAcidSequence(BiologicalSequence):
    """
    Abstract base class for nucleic acid sequences, derived from BiologicalSequence.

    Implements generic nucleic acid sequence functionalities, which can be extended by specific types of nucleic acids.

    Attributes:
        complement_map (dict): A dictionary mapping each nucleotide to its complement.

    Methods:
        is_valid(): Checks if all nucleotides in the sequence are valid according to complement_map.
        complement(): Returns the complementary sequence.
        gc_content(): Calculates the GC content of the sequence.
    """
    complement_map = {}

    def is_valid(self) -> bool:
        return all(nucleotide in self.complement_map for nucleotide in self.sequence)

    def complement(self):
        return ''.join(self.complement_map[nucleotide] for nucleotide in self.sequence)

    def gc_content(self):
        gc_content = (self.sequence.count('G') + self.sequence.count('C')) / len(self.sequence) if self.sequence else 0
        return gc_content

class DNASequence(NucleicAcidSequence):
    """
    Represents a DNA sequence, derived from NucleicAcidSequence.

    Contains specific methods and attributes for DNA, including transcription.

    Attributes:
        complement_map (dict): Maps each DNA nucleotide to its complement.

    Methods:
        transcribe(): Converts the DNA sequence to an RNA sequence.
    """
    complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    def transcribe(self):
        return RNASequence(self.sequence.replace('T', 'U'))

class RNASequence(NucleicAcidSequence):
    """
    Represents an RNA sequence, derived from NucleicAcidSequence.

    Contains specific methods and attributes for RNA.

    Attributes:
        complement_map (dict): Maps each RNA nucleotide to its complement.
    """
    complement_map = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}


class AminoAcidSequence(BiologicalSequence):
    """
    Represents an amino acid sequence, derived from BiologicalSequence.

    Contains specific methods and attributes for amino acid sequences.

    Methods:
        is_valid(): Checks if all characters in the sequence are valid amino acids.
        one_to_three_letter_code(): Converts the sequence from one-letter code to three-letter code.
    """

    def is_valid(self) -> bool:
        amino_acids = "ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy"
        return all(aa in amino_acids for aa in self.sequence)
    

    def one_to_three_letter_code(self) -> str:
        """
        This function converts a protein sequence from one-letter amino acid code to three-letter code.
    
        Args:
            sequence (str): The input protein sequence in one-letter code.
        
        Returns:
            str: The converted protein sequence in three-letter code.
        """
        AMINO_ACIDS = {'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
               'K': 'Lys', 'L': 'Leu', 'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg', 'S': 'Ser',
               'T': 'Thr', 'V': 'Val',
               'W': 'Trp', 'Y': 'Tyr'}
        three_letter_code = [AMINO_ACIDS.get(aa.upper()) for aa in self.sequence]
        return '-'.join(three_letter_code)
    
    
def filter_fastq(input_path: str, output_filename: str = None, gc_bounds: tuple = (0, 100), length_bounds: tuple = (0, 2 ** 32), quality_threshold: int = 0) -> dict:
    """
    Filters FASTQ sequences from fastq format file based on specified criteria. 
    Saves the output FASTAQ file. Uses Biopython libraries.
    """
    filtered_seqs = {}
    
    for record in SeqIO.parse(input_path, "fastq"):
        sequence = str(record.seq)
        quality_scores = record.letter_annotations["phred_quality"]
        gc_content =  GC(record.seq)*100

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

#HW_17
import os
import sys
import io
import time
from datetime import timedelta
import requests
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Retrieving the API token and chat ID from environment variables
API_TOKEN = os.getenv("API_TOKEN")

def telegram_logger(chat_id):
    """
    A decorator for logging the execution of a function and sending the results to Telegram.

    Args:
        chat_id (str): Telegram chat ID where messages will be sent.

    Returns:
        Callable: The decorated function.
    """
    def decorator(func):
        def wrapper(*args, **kwargs):
            start_time = time.time()
            stdout_backup = sys.stdout
            stderr_backup = sys.stderr
            sys.stdout = io.StringIO()  # Redirecting stdout
            sys.stderr = io.StringIO()  # Redirecting stderr

            try:
                result = func(*args, **kwargs)
                duration = time.time() - start_time
                stdout_content = sys.stdout.getvalue()
                stderr_content = sys.stderr.getvalue()
                duration = time.time() - start_time

                if duration < 86400:
                    duration_str = str(timedelta(seconds=duration))
                else:
                    days = duration // 86400
                    duration_str = f"{days} days, {str(timedelta(seconds=duration % 86400))}"

                message = f"Function `{func.__name__}` completed in {duration_str}"
                send_message(chat_id, message)
                return result
            except Exception as e:
                message = f"Function `{func.__name__}` failed with error: {type(e).__name__}: {e}"
                send_message(chat_id, message)
                raise e
            finally:
                stdout_content = sys.stdout.getvalue()
                stderr_content = sys.stderr.getvalue()

                if stdout_content or stderr_content:
                    send_document(chat_id, stdout_content + stderr_content, f"{func.__name__}.log")

                sys.stdout = stdout_backup
                sys.stderr = stderr_backup

        return wrapper
    return decorator

def send_message(chat_id, text):
    """
    Sends a text message to the specified Telegram chat.

    Args:
        chat_id (str): Telegram chat ID where the message will be sent.
        text (str): The text of the message to be sent.
    """
    url = f"https://api.telegram.org/bot{API_TOKEN}/sendMessage"
    data = {'chat_id': chat_id, 'text': text, 'parse_mode': 'Markdown'}
    requests.post(url, data=data)

def send_document(chat_id, document_content, document_name):
    """
    Sends a document to the specified Telegram chat.

    Args:
        chat_id (str): Telegram chat ID where the document will be sent.
        document_content (str): The content of the document to be sent.
        document_name (str): The file name of the document to be sent.

    Notes:
        The document is sent as a virtual file created from a string.
    """
    url = f"https://api.telegram.org/bot{API_TOKEN}/sendDocument"
    files = {'document': (document_name, document_content)}
    data = {'chat_id': chat_id}
    requests.post(url, data=data, files=files)