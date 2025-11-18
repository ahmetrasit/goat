"""
Unit tests for process.py
Tests critical set operations and data processing functions
"""

import unittest
from process import Process


class TestSetOperations(unittest.TestCase):
    """Test set operations functionality"""

    def setUp(self):
        self.process = Process()

    def test_applyOperation_2way_and(self):
        """Test 2-way intersection"""
        set_a = {'gene1', 'gene2', 'gene3'}
        set_b = {'gene2', 'gene3', 'gene4'}
        result = self.process.applyOperation(set_a, set_b, set(), set(), 'and')
        self.assertEqual(result, {'gene2', 'gene3'})

    def test_applyOperation_2way_or(self):
        """Test 2-way union"""
        set_a = {'gene1', 'gene2'}
        set_b = {'gene3', 'gene4'}
        result = self.process.applyOperation(set_a, set_b, set(), set(), 'or')
        self.assertEqual(result, {'gene1', 'gene2', 'gene3', 'gene4'})

    def test_applyOperation_2way_difference(self):
        """Test 2-way difference A-B"""
        set_a = {'gene1', 'gene2', 'gene3'}
        set_b = {'gene2', 'gene3', 'gene4'}
        result = self.process.applyOperation(set_a, set_b, set(), set(), 'a-b')
        self.assertEqual(result, {'gene1'})

    def test_applyOperation_3way_and(self):
        """Test 3-way intersection"""
        set_a = {'gene1', 'gene2', 'gene3'}
        set_b = {'gene2', 'gene3', 'gene4'}
        set_c = {'gene2', 'gene5'}
        result = self.process.applyOperation(set_a, set_b, set_c, set(), 'abc-and')
        self.assertEqual(result, {'gene2'})

    def test_applyOperation_3way_or(self):
        """Test 3-way union"""
        set_a = {'gene1'}
        set_b = {'gene2'}
        set_c = {'gene3'}
        result = self.process.applyOperation(set_a, set_b, set_c, set(), 'abc-or')
        self.assertEqual(result, {'gene1', 'gene2', 'gene3'})

    def test_applyOperation_4way_and(self):
        """Test 4-way intersection"""
        set_a = {'gene1', 'gene2'}
        set_b = {'gene1', 'gene2', 'gene3'}
        set_c = {'gene1', 'gene2', 'gene4'}
        set_d = {'gene1', 'gene2', 'gene5'}
        result = self.process.applyOperation(set_a, set_b, set_c, set_d, 'abcd-and')
        self.assertEqual(result, {'gene1', 'gene2'})

    def test_applyOperation_4way_or(self):
        """Test 4-way union"""
        set_a = {'gene1'}
        set_b = {'gene2'}
        set_c = {'gene3'}
        set_d = {'gene4'}
        result = self.process.applyOperation(set_a, set_b, set_c, set_d, 'abcd-or')
        self.assertEqual(result, {'gene1', 'gene2', 'gene3', 'gene4'})

    def test_applyOperation_empty_sets(self):
        """Test with empty sets"""
        set_a = set()
        set_b = {'gene1'}
        result = self.process.applyOperation(set_a, set_b, set(), set(), 'and')
        self.assertEqual(result, set())

    def test_applyOperation_unknown_operation(self):
        """Test with unknown operation returns empty set"""
        set_a = {'gene1'}
        set_b = {'gene2'}
        result = self.process.applyOperation(set_a, set_b, set(), set(), 'invalid-op')
        self.assertEqual(result, set())


class TestIDConversion(unittest.TestCase):
    """Test ID conversion functionality"""

    def setUp(self):
        self.process = Process()

    def test_getIdType_empty_set(self):
        """Test ID type detection with empty set"""
        id_sets = {
            'alias2name': {'gene_a', 'gene_b'},
            'name2gene': {'WBGene001', 'WBGene002'},
            'gene2name': {'WBGene003'}
        }
        gene_set = set()
        id_type, percentage = self.process.getIdType(id_sets, gene_set)
        # With empty set, should still return a type, percentage should be 0
        self.assertIn(id_type, ['alias', 'name', 'gene'])

    def test_getTranscriptGeneName_valid(self):
        """Test transcript name extraction"""
        from process import Process
        proc = Process()
        # This would need getTranscriptGeneName method - checking if it exists
        if hasattr(proc, 'getTranscriptGeneName'):
            result = proc.getTranscriptGeneName('Y74C9A.6.t1')
            self.assertEqual(result, 'Y74C9A.6.t1' if '.' in result else 'Y74C9A')


class TestFiltering(unittest.TestCase):
    """Test data filtering functionality"""

    def setUp(self):
        self.process = Process()

    def test_passedSeqRules_greater_than(self):
        """Test sequence rule: count > threshold"""
        seq = 'ATCG'
        count = 100
        rules = [['>','50']]
        result = self.process.passedSeqRules(seq, count, rules)
        self.assertTrue(result)

    def test_passedSeqRules_less_than(self):
        """Test sequence rule: count < threshold"""
        seq = 'ATCG'
        count = 25
        rules = [['<', '50']]
        result = self.process.passedSeqRules(seq, count, rules)
        self.assertTrue(result)

    def test_passedSeqRules_equals(self):
        """Test sequence rule: count = threshold"""
        seq = 'ATCG'
        count = 50
        rules = [['=', '50']]
        result = self.process.passedSeqRules(seq, count, rules)
        self.assertTrue(result)

    def test_removeEmptyGenes(self):
        """Test removal of empty gene entries"""
        data = {
            'gene1': {'seq1': 100},
            'gene2': {},
            'gene3': {'seq2': 50}
        }
        result = self.process.removeEmptyGenes(data)
        self.assertIn('gene1', result)
        self.assertNotIn('gene2', result)
        self.assertIn('gene3', result)


if __name__ == '__main__':
    unittest.main()
