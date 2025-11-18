"""
Utility functions for G.O.A.T
Includes validation, security, and helper functions
"""

import os
import re
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Any
from markupsafe import escape
from constants import (
    MIN_SEQUENCE_LENGTH, MAX_SEQUENCE_LENGTH,
    VALID_NUCLEOTIDES, MAX_DATASET_NAME_LENGTH,
    VALID_FILE_EXTENSIONS, FILENAME_COUNTER_START,
    MAX_FILENAME_RENAME_ATTEMPTS
)

logger = logging.getLogger(__name__)


# ============================================================================
# Security and Validation
# ============================================================================

def validate_path(filepath: str, allowed_base_dirs: List[str]) -> bool:
    """
    Validate that a filepath is within allowed directories (prevents path traversal)

    Args:
        filepath: Path to validate
        allowed_base_dirs: List of allowed base directory paths

    Returns:
        True if path is safe, False otherwise

    Example:
        >>> validate_path('data/binned/sample.json', ['data'])
        True
        >>> validate_path('../etc/passwd', ['data'])
        False
    """
    try:
        # Resolve to absolute path
        abs_path = Path(filepath).resolve()

        # Check if within any allowed directory
        for base_dir in allowed_base_dirs:
            abs_base = Path(base_dir).resolve()
            try:
                abs_path.relative_to(abs_base)
                return True
            except ValueError:
                continue

        logger.warning(f"Path traversal attempt blocked: {filepath}")
        return False
    except Exception as e:
        logger.error(f"Error validating path {filepath}: {e}")
        return False


def sanitize_filename(filename: str, max_length: int = MAX_DATASET_NAME_LENGTH) -> str:
    """
    Sanitize filename to prevent security issues

    Args:
        filename: Original filename
        max_length: Maximum allowed length

    Returns:
        Sanitized filename
    """
    # Remove path components
    filename = os.path.basename(filename)

    # Escape HTML
    filename = str(escape(filename))

    # Remove special characters, keep alphanumeric, dash, underscore, dot
    filename = re.sub(r'[^\w\s.-]', '', filename.lower())

    # Replace spaces and multiple dashes/underscores
    filename = re.sub(r'[-\s]+', '-', filename).strip('-_')

    # Limit length
    if len(filename) > max_length:
        name, ext = os.path.splitext(filename)
        filename = name[:max_length - len(ext)] + ext

    return filename


def validate_sequence_length(length: int) -> Tuple[bool, str]:
    """
    Validate sequence length is within acceptable range

    Args:
        length: Sequence length to validate

    Returns:
        Tuple of (is_valid, error_message)
    """
    if not isinstance(length, int):
        return False, f"Length must be an integer, got {type(length)}"

    if length < MIN_SEQUENCE_LENGTH:
        return False, f"Length {length} is below minimum {MIN_SEQUENCE_LENGTH}"

    if length > MAX_SEQUENCE_LENGTH:
        return False, f"Length {length} exceeds maximum {MAX_SEQUENCE_LENGTH}"

    return True, ""


def validate_nucleotides(nucleotides: str) -> Tuple[bool, str]:
    """
    Validate nucleotide string contains only valid bases

    Args:
        nucleotides: String of nucleotides (e.g., 'ATGC')

    Returns:
        Tuple of (is_valid, error_message)
    """
    if not nucleotides:
        return False, "Nucleotide string cannot be empty"

    invalid_chars = set(nucleotides.upper()) - VALID_NUCLEOTIDES
    if invalid_chars:
        return False, f"Invalid nucleotides: {invalid_chars}"

    return True, ""


def validate_file_extension(filename: str, allowed_extensions: Optional[Set[str]] = None) -> bool:
    """
    Validate file has an allowed extension

    Args:
        filename: Filename to check
        allowed_extensions: Set of allowed extensions (uses default if None)

    Returns:
        True if extension is allowed
    """
    if allowed_extensions is None:
        allowed_extensions = VALID_FILE_EXTENSIONS

    ext = os.path.splitext(filename)[1].lower()

    # Handle .fastq.gz, .fa.gz etc.
    if filename.endswith('.gz'):
        base_ext = os.path.splitext(filename[:-3])[1].lower()
        return f"{base_ext}.gz" in allowed_extensions or ext in allowed_extensions

    return ext in allowed_extensions


# ============================================================================
# File Operations with Error Handling
# ============================================================================

def safe_load_json(filepath: str, default: Any = None) -> Any:
    """
    Safely load JSON file with comprehensive error handling

    Args:
        filepath: Path to JSON file
        default: Default value to return on error

    Returns:
        Loaded JSON data or default value
    """
    try:
        if not os.path.exists(filepath):
            logger.warning(f"File not found: {filepath}")
            return default

        with open(filepath, 'r') as f:
            return json.load(f)

    except json.JSONDecodeError as e:
        logger.error(f"Invalid JSON in {filepath}: {e}")
        return default
    except PermissionError as e:
        logger.error(f"Permission denied reading {filepath}: {e}")
        return default
    except Exception as e:
        logger.error(f"Error loading {filepath}: {e}")
        return default


def safe_save_json(filepath: str, data: Any, indent: Optional[int] = None) -> bool:
    """
    Safely save JSON file with error handling

    Args:
        filepath: Path to save JSON file
        data: Data to serialize
        indent: JSON indentation (None for compact)

    Returns:
        True if successful, False otherwise
    """
    try:
        # Ensure directory exists
        os.makedirs(os.path.dirname(filepath), exist_ok=True)

        # Write atomically using temporary file
        temp_path = f"{filepath}.tmp"
        with open(temp_path, 'w') as f:
            json.dump(data, f, indent=indent)

        # Rename to final location (atomic on POSIX)
        os.replace(temp_path, filepath)

        logger.info(f"Saved JSON to {filepath}")
        return True

    except PermissionError as e:
        logger.error(f"Permission denied writing {filepath}: {e}")
        return False
    except Exception as e:
        logger.error(f"Error saving {filepath}: {e}")
        return False


def get_unique_filename(base_path: str, extension: str = '.json') -> str:
    """
    Generate unique filename with proper counter (fixes race condition)

    Args:
        base_path: Base path without extension
        extension: File extension

    Returns:
        Unique filepath that doesn't exist

    Raises:
        RuntimeError: If cannot find unique name after max attempts
    """
    # Clean base path
    base_path = re.sub(r'\.json$', '', base_path)
    base_path = sanitize_filename(base_path)

    # Try without counter first
    filepath = f"{base_path}{extension}"
    if not os.path.exists(filepath):
        return filepath

    # Try with counter
    for i in range(FILENAME_COUNTER_START, MAX_FILENAME_RENAME_ATTEMPTS):
        filepath = f"{base_path}_{i}{extension}"
        if not os.path.exists(filepath):
            return filepath

    raise RuntimeError(f"Could not generate unique filename after {MAX_FILENAME_RENAME_ATTEMPTS} attempts")


# ============================================================================
# Data Validation
# ============================================================================

def validate_form_data(formdata: Dict, required_fields: List[str]) -> Tuple[bool, List[str]]:
    """
    Validate that required fields are present in form data

    Args:
        formdata: Form data dictionary
        required_fields: List of required field names

    Returns:
        Tuple of (is_valid, missing_fields)
    """
    missing = [field for field in required_fields if field not in formdata or not formdata[field]]
    return len(missing) == 0, missing


def validate_length_range(len_start: int, len_end: int) -> Tuple[bool, str]:
    """
    Validate sequence length range

    Args:
        len_start: Start length (minimum)
        len_end: End length (maximum)

    Returns:
        Tuple of (is_valid, error_message)
    """
    # Validate individual lengths
    valid_start, msg_start = validate_sequence_length(len_start)
    if not valid_start:
        return False, f"Start length: {msg_start}"

    valid_end, msg_end = validate_sequence_length(len_end)
    if not valid_end:
        return False, f"End length: {msg_end}"

    # Validate range
    if len_start > len_end:
        return False, f"Start length ({len_start}) cannot exceed end length ({len_end})"

    return True, ""


# ============================================================================
# String and Data Processing
# ============================================================================

def parse_rule_string(rule_str: str) -> Optional[List[str]]:
    """
    Parse rule string from form data

    Args:
        rule_str: Rule string (e.g., "Top;25" or ">;100")

    Returns:
        List of rule components or None if invalid
    """
    try:
        parts = rule_str.split(';')
        if len(parts) >= 2:
            return parts
        return None
    except Exception as e:
        logger.error(f"Error parsing rule string '{rule_str}': {e}")
        return None


def format_gene_count(count: int) -> str:
    """
    Format gene count for display

    Args:
        count: Number of genes

    Returns:
        Formatted string
    """
    if count == 0:
        return "No genes"
    elif count == 1:
        return "1 gene"
    else:
        return f"{count:,} genes"


def truncate_string(s: str, max_length: int = 50, suffix: str = "...") -> str:
    """
    Truncate string to maximum length

    Args:
        s: String to truncate
        max_length: Maximum length
        suffix: Suffix to add when truncated

    Returns:
        Truncated string
    """
    if len(s) <= max_length:
        return s
    return s[:max_length - len(suffix)] + suffix


# ============================================================================
# Math and Statistics Helpers
# ============================================================================

def safe_divide(numerator: float, denominator: float, default: float = 0.0) -> float:
    """
    Safely divide with default for division by zero

    Args:
        numerator: Numerator value
        denominator: Denominator value
        default: Value to return if denominator is zero

    Returns:
        Result of division or default
    """
    try:
        if denominator == 0:
            return default
        return numerator / denominator
    except (TypeError, ZeroDivisionError):
        return default


def calculate_percentage(part: float, total: float, decimals: int = 1) -> float:
    """
    Calculate percentage with rounding

    Args:
        part: Part value
        total: Total value
        decimals: Number of decimal places

    Returns:
        Percentage rounded to specified decimals
    """
    if total == 0:
        return 0.0
    return round((part / total) * 100, decimals)


# ============================================================================
# Logging Helpers
# ============================================================================

def setup_logging(log_level: str = 'INFO', log_file: Optional[str] = None) -> None:
    """
    Configure logging for the application

    Args:
        log_level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        log_file: Optional log file path
    """
    from constants import LOG_FORMAT_DETAILED, LOG_DATE_FORMAT

    level = getattr(logging, log_level.upper(), logging.INFO)

    handlers = [logging.StreamHandler()]
    if log_file:
        handlers.append(logging.FileHandler(log_file))

    logging.basicConfig(
        level=level,
        format=LOG_FORMAT_DETAILED,
        datefmt=LOG_DATE_FORMAT,
        handlers=handlers
    )

    logger.info(f"Logging configured: level={log_level}, file={log_file}")


# ============================================================================
# Performance Helpers
# ============================================================================

def chunked_iterable(iterable, chunk_size: int = 1000):
    """
    Yield successive chunks from iterable

    Args:
        iterable: Any iterable
        chunk_size: Size of each chunk

    Yields:
        Chunks of the iterable
    """
    chunk = []
    for item in iterable:
        chunk.append(item)
        if len(chunk) >= chunk_size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk
