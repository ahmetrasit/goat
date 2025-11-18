"""
Unit tests for utils.py
Tests validation and security functions
"""

import unittest
import tempfile
import os
from utils import (
    validate_path, sanitize_filename, validate_sequence_length,
    validate_nucleotides, safe_load_json, safe_save_json,
    get_unique_filename, validate_length_range, safe_divide,
    calculate_percentage
)


class TestPathValidation(unittest.TestCase):
    """Test path validation and security"""

    def test_validate_path_safe(self):
        """Test validation of safe path"""
        result = validate_path('data/binned/sample.json', ['data'])
        self.assertTrue(result)

    def test_validate_path_traversal_blocked(self):
        """Test path traversal attack is blocked"""
        result = validate_path('../../../etc/passwd', ['data'])
        self.assertFalse(result)

    def test_validate_path_absolute_blocked(self):
        """Test absolute path outside allowed dirs is blocked"""
        result = validate_path('/etc/passwd', ['data'])
        self.assertFalse(result)

    def test_sanitize_filename_removes_special_chars(self):
        """Test filename sanitization removes special characters"""
        result = sanitize_filename('test<>file?.json')
        self.assertNotIn('<', result)
        self.assertNotIn('>', result)
        self.assertNotIn('?', result)

    def test_sanitize_filename_handles_spaces(self):
        """Test filename sanitization handles spaces"""
        result = sanitize_filename('my test file.json')
        self.assertNotIn(' ', result)

    def test_sanitize_filename_max_length(self):
        """Test filename length limiting"""
        long_name = 'a' * 300 + '.json'
        result = sanitize_filename(long_name, max_length=50)
        self.assertLessEqual(len(result), 50)


class TestSequenceValidation(unittest.TestCase):
    """Test sequence validation functions"""

    def test_validate_sequence_length_valid(self):
        """Test valid sequence length"""
        valid, msg = validate_sequence_length(21)
        self.assertTrue(valid)
        self.assertEqual(msg, "")

    def test_validate_sequence_length_too_short(self):
        """Test sequence length below minimum"""
        valid, msg = validate_sequence_length(5)
        self.assertFalse(valid)
        self.assertIn('below minimum', msg)

    def test_validate_sequence_length_too_long(self):
        """Test sequence length above maximum"""
        valid, msg = validate_sequence_length(50)
        self.assertFalse(valid)
        self.assertIn('exceeds maximum', msg)

    def test_validate_sequence_length_wrong_type(self):
        """Test sequence length with wrong type"""
        valid, msg = validate_sequence_length("21")
        self.assertFalse(valid)
        self.assertIn('must be an integer', msg)

    def test_validate_nucleotides_valid(self):
        """Test valid nucleotide string"""
        valid, msg = validate_nucleotides('ATGC')
        self.assertTrue(valid)

    def test_validate_nucleotides_invalid_chars(self):
        """Test nucleotide string with invalid characters"""
        valid, msg = validate_nucleotides('ATGCX')
        self.assertFalse(valid)
        self.assertIn('X', msg)

    def test_validate_nucleotides_empty(self):
        """Test empty nucleotide string"""
        valid, msg = validate_nucleotides('')
        self.assertFalse(valid)

    def test_validate_length_range_valid(self):
        """Test valid length range"""
        valid, msg = validate_length_range(21, 23)
        self.assertTrue(valid)

    def test_validate_length_range_reversed(self):
        """Test reversed length range"""
        valid, msg = validate_length_range(23, 21)
        self.assertFalse(valid)
        self.assertIn('cannot exceed', msg)


class TestFileOperations(unittest.TestCase):
    """Test file operation functions"""

    def setUp(self):
        """Create temporary directory for tests"""
        self.test_dir = tempfile.mkdtemp()

    def tearDown(self):
        """Clean up temporary directory"""
        import shutil
        shutil.rmtree(self.test_dir, ignore_errors=True)

    def test_safe_save_and_load_json(self):
        """Test safe JSON save and load"""
        test_file = os.path.join(self.test_dir, 'test.json')
        test_data = {'gene1': 100, 'gene2': 200}

        # Save
        result = safe_save_json(test_file, test_data)
        self.assertTrue(result)
        self.assertTrue(os.path.exists(test_file))

        # Load
        loaded = safe_load_json(test_file)
        self.assertEqual(loaded, test_data)

    def test_safe_load_json_nonexistent(self):
        """Test loading non-existent file returns default"""
        result = safe_load_json('/nonexistent/file.json', default={})
        self.assertEqual(result, {})

    def test_get_unique_filename_no_conflict(self):
        """Test getting unique filename when no conflict"""
        base_path = os.path.join(self.test_dir, 'test')
        result = get_unique_filename(base_path, '.json')
        self.assertEqual(result, base_path + '.json')

    def test_get_unique_filename_with_conflict(self):
        """Test getting unique filename when file exists"""
        base_path = os.path.join(self.test_dir, 'test')

        # Create existing file
        with open(base_path + '.json', 'w') as f:
            f.write('{}')

        # Should get _1 suffix
        result = get_unique_filename(base_path, '.json')
        self.assertEqual(result, base_path + '_1.json')


class TestMathHelpers(unittest.TestCase):
    """Test math and statistics helpers"""

    def test_safe_divide_normal(self):
        """Test safe division with normal values"""
        result = safe_divide(10, 2)
        self.assertEqual(result, 5.0)

    def test_safe_divide_by_zero(self):
        """Test safe division by zero returns default"""
        result = safe_divide(10, 0, default=0.0)
        self.assertEqual(result, 0.0)

    def test_safe_divide_custom_default(self):
        """Test safe division with custom default"""
        result = safe_divide(10, 0, default=-1.0)
        self.assertEqual(result, -1.0)

    def test_calculate_percentage(self):
        """Test percentage calculation"""
        result = calculate_percentage(25, 100)
        self.assertEqual(result, 25.0)

    def test_calculate_percentage_zero_total(self):
        """Test percentage with zero total"""
        result = calculate_percentage(25, 0)
        self.assertEqual(result, 0.0)

    def test_calculate_percentage_rounding(self):
        """Test percentage rounding"""
        result = calculate_percentage(1, 3, decimals=2)
        self.assertEqual(result, 33.33)


if __name__ == '__main__':
    unittest.main()
