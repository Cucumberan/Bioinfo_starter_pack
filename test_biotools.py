import unittest
from biotools import DNASequence, RNASequence, check_quality
from bio_files_processor import  OpenFasta, FastaRecord

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


if __name__ == '__main__':
    unittest.main()
