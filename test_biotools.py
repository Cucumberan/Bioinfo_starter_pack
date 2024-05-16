import unittest
from biotools import DNASequence, RNASequence, check_quality
from bio_files_processor import  OpenFasta, FastaRecord
import os

class TestBiotools(unittest.TestCase):
    def test_sequence_length(self):
        """
        Test the __len__ method to ensure it returns the correct sequence length.
        """
        seq = DNASequence("AGTC")
        self.assertEqual(len(seq), 4)

class TestDNA(unittest.TestCase):
    def test_dna_validation(self):
        """
        Test that `is_valid` raises ValueError when invalid nucleotides are included in the DNA sequence.
        """
        seq = DNASequence("AGTX")  # Assuming 'X' is an invalid nucleotide.
        with self.assertRaises(ValueError):
            seq.is_valid()

    def test_gc_content(self):
        """
        Test the gc_content method in DNASequence to verify correct GC percentage calculation.
        """
        seq = DNASequence("GCGC")
        self.assertEqual(seq.gc_content(), 1)

    def test_transcription(self):
        '''
        Test the transcribe method in DNASequence to ensure DNA is correctly transcribed to RNA.
        '''
        dna = DNASequence("ATGC")
        rna = dna.transcribe()
        self.assertIsInstance(rna, RNASequence)
        self.assertEqual(str(rna), "AUGC")


class TestComplement(unittest.TestCase):
    def test_dna_complement(self):
        """
        Test the complement method in DNASequence to ensure it returns the correct complementary sequence.
        """
        seq = DNASequence("ATGC")
        self.assertEqual(seq.complement(), "TACG")


class TesReadFasta(unittest.TestCase):
    def test_read_fasta_file(self):
        """
        Test reading from a FASTA file using OpenFasta to ensure sequences are correctly parsed and returned as FastaRecord objects.
        """
        with OpenFasta("data/example.fasta") as fasta_file:
            records = list(fasta_file.read_records())
            self.assertTrue(len(records) > 0)
            self.assertIsInstance(records[0], FastaRecord)


class TestQualityCheck(unittest.TestCase):
    def test_invalid_quality_score_type(self):
        """
        Test that check_quality raises ValueError when quality_scores is not a string or list/tuple.
        """
        quality_scores_dict = {'score': 30}  # Incorrect datatype for quality_scores
        with self.assertRaises(ValueError):
            check_quality(quality_scores_dict, 20)


class TestConvertFasta(unittest.TestCase):
    """
    A unit test class for testing the convert_multiline_fasta_to_oneline function from the bio_files_processor module.

    This test class ensures that the conversion from a multiline FASTA format to a single-line FASTA format is performed correctly. The class uses the unittest framework to manage setup, execution, and teardown of tests.

    Methods:
        setUp(self): Prepares the environment before each test function is executed. It creates a sample multiline FASTA file to be used in the tests.
        test_conversion(self): Executes the test to verify that the convert_multiline_fasta_to_oneline function correctly converts a given multiline FASTA file to the expected single-line format. The test checks that the content of the converted file matches the expected results.
        tearDown(self): Cleans up the environment after each test function has been executed. It removes any files created during the test setup to ensure a clean state for subsequent tests.

    Usage:
        Running this test class checks the functionality of the convert_multiline_fasta_to_oneline function to ensure it accurately handles the conversion process without data loss or format errors.
    """
    def setUp(self):
        self.input_fasta = "test.fasta"
        self.output_fasta = "test_one_line.fasta"
        with open(self.input_fasta, 'w') as f:
            f.write(">seq1\nATGCA\nTGACG\n>seq2\nTTGAC\nGTACA\n")

    def test_conversion(self):
        # Import the function to test
        from bio_files_processor import convert_multiline_fasta_to_oneline
        convert_multiline_fasta_to_oneline(self.input_fasta, self.output_fasta)
        expected_results = ">seq1\nATGCATGACG\n>seq2\nTTGACGTACA\n"
        with open(self.output_fasta, 'r') as f:
            results = f.read()
        self.assertEqual(results, expected_results)

    def tearDown(self):
        # Clean up: Remove files created for the test
        os.remove(self.input_fasta)
        os.remove(self.output_fasta)

if __name__ == '__main__':
    unittest.main()
