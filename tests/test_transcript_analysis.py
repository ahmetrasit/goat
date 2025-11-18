"""
Unit tests for TranscriptAnalysis.py
Tests SAM parsing and multi-mapper categorization
"""

import unittest
from TranscriptAnalysis import TranscriptAnalysis


class TestTranscriptAnalysis(unittest.TestCase):
    """Test TranscriptAnalysis functionality"""

    def test_flag2strand_forward(self):
        """Test SAM flag to strand conversion for forward"""
        ta = TranscriptAnalysis(None)
        result = ta.flag2strand(0)
        self.assertEqual(result, 's')

    def test_flag2strand_reverse(self):
        """Test SAM flag to strand conversion for reverse"""
        ta = TranscriptAnalysis(None)
        result = ta.flag2strand(16)
        self.assertEqual(result, 'a')

    def test_getTranscriptGeneName_valid(self):
        """Test transcript name extraction"""
        ta = TranscriptAnalysis(None)
        result = ta.getTranscriptGeneName('TRANSCR:Y74C9A.6.t1:WBGene12345')
        # Should extract the transcript ID
        self.assertIsNotNone(result)
        self.assertNotEqual(result, '-')

    def test_getTranscriptGeneName_invalid(self):
        """Test transcript name extraction with invalid input"""
        ta = TranscriptAnalysis(None)
        result = ta.getTranscriptGeneName('invalid')
        self.assertEqual(result, '-')

    def test_findMappers_unique(self):
        """Test identification of unique mappers"""
        ta = TranscriptAnalysis(None)
        seq2genes = {
            'ATCGATCG': {'s': {'gene1'}},
            'GCTAGCTA': {'s': {'gene2'}}
        }
        unique, multi = ta.findMappers(seq2genes)
        self.assertEqual(len(unique), 2)
        self.assertEqual(len(multi), 0)

    def test_findMappers_multi(self):
        """Test identification of multi-mappers"""
        ta = TranscriptAnalysis(None)
        seq2genes = {
            'ATCGATCGATCG': {'s': {'gene1', 'gene2'}},  # Multi-mapper
            'GCTAGCTAGCTA': {'s': {'gene3'}}  # Unique
        }
        unique, multi = ta.findMappers(seq2genes)
        self.assertEqual(len(multi), 1)
        self.assertIn('ATCGATCGATCG', multi)

    def test_filterSeqBySpecies_correct_logic(self):
        """Test sequence filtering uses correct length logic"""
        ta = TranscriptAnalysis(None)
        seq_set = {
            'ATCGATCGATCGATCGATCGA',  # 21nt
            'ATCGATCGATCGATCGATCGAT',  # 22nt
            'ATCGATCGATCGATCGATCGATC',  # 23nt
            'ATCGATCGATCGATCGATCGATCG'  # 24nt
        }
        norm_seq2ppm = {seq: 100.0 for seq in seq_set}

        # Filter for 21-23nt G-starting sequences
        filtered, has_results = ta.filterSeqBySpecies(seq_set, 21, 23, 'A', norm_seq2ppm)

        # Should have 21, 22, 23nt sequences (3 total)
        self.assertEqual(len(filtered), 3)
        self.assertTrue(has_results)

    def test_filterSeqBySpecies_nucleotide_filter(self):
        """Test sequence filtering by nucleotide"""
        ta = TranscriptAnalysis(None)
        seq_set = {
            'ATCGATCGATCGATCGATCGA',  # Starts with A
            'GTCGATCGATCGATCGATCGA',  # Starts with G
            'CTCGATCGATCGATCGATCGA'  # Starts with C
        }
        norm_seq2ppm = {seq: 100.0 for seq in seq_set}

        # Filter for G-starting sequences
        filtered, _ = ta.filterSeqBySpecies(seq_set, 21, 21, 'G', norm_seq2ppm)

        # Should have only G-starting sequence
        self.assertEqual(len(filtered), 1)
        self.assertIn('GTCGATCGATCGATCGATCGA', filtered)


class TestNormalization(unittest.TestCase):
    """Test normalization methods"""

    def test_normAll_identity(self):
        """Test 'all' normalization returns unchanged"""
        ta = TranscriptAnalysis(None)
        gene2total_ppm = {'gene1': 100.0, 'gene2': 200.0}
        seq2ppm = {'seq1': 50.0, 'seq2': 75.0}

        norm_gene, norm_seq = ta.normAll(gene2total_ppm, seq2ppm)

        self.assertEqual(norm_gene, gene2total_ppm)
        self.assertEqual(norm_seq, seq2ppm)

    def test_sumOfType_calculation(self):
        """Test sum calculation by gene type"""
        ta = TranscriptAnalysis(None)
        ta.gene2element = {'gene1': 'MIRNA', 'gene2': 'MIRNA', 'gene3': 'PC'}

        gene2total_ppm = {
            'gene1': {'a': {'u': 100, 'm': 0}, 's': {'u': 50, 'm': 0}},
            'gene2': {'a': {'u': 200, 'm': 0}, 's': {'u': 100, 'm': 0}},
            'gene3': {'a': {'u': 300, 'm': 0}, 's': {'u': 150, 'm': 0}}
        }

        # Sum of MIRNA type
        result = ta.sumOfType(gene2total_ppm, 'MIRNA')

        # gene1: 100+50 = 150, gene2: 200+100 = 300, total = 450
        self.assertEqual(result, 450)


if __name__ == '__main__':
    unittest.main()
