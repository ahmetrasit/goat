"""
Security extensions and middleware for G.O.A.T
Includes CSRF protection, rate limiting, and security headers
"""

from flask import request
from flask_wtf.csrf import CSRFProtect
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address
from flask_cors import CORS
import logging
from constants import RATE_LIMIT_PREPROCESS, RATE_LIMIT_GENERAL
from utils import validate_path

logger = logging.getLogger(__name__)

# CSRF Protection
csrf = CSRFProtect()

# Rate Limiting
limiter = Limiter(
    key_func=get_remote_address,
    default_limits=[RATE_LIMIT_GENERAL],
    storage_uri="memory://",  # Use Redis in production: "redis://localhost:6379"
    strategy="fixed-window"
)


def init_security(app):
    """
    Initialize security extensions for Flask app

    Args:
        app: Flask application instance
    """
    # CSRF Protection
    csrf.init_app(app)
    logger.info("CSRF protection enabled")

    # Rate Limiting
    limiter.init_app(app)
    logger.info("Rate limiting enabled")

    # CORS (if needed for API access)
    if app.config.get('ENABLE_CORS', False):
        CORS(app, resources={
            r"/api/*": {"origins": app.config.get('CORS_ORIGINS', '*')}
        })
        logger.info("CORS enabled")

    # Security Headers
    @app.after_request
    def add_security_headers(response):
        """Add security headers to all responses"""
        response.headers['X-Content-Type-Options'] = 'nosniff'
        response.headers['X-Frame-Options'] = 'SAMEORIGIN'
        response.headers['X-XSS-Protection'] = '1; mode=block'

        # Strict-Transport-Security (only in production with HTTPS)
        if not app.debug:
            response.headers['Strict-Transport-Security'] = 'max-age=31536000; includeSubDomains'

        # Content Security Policy
        csp = "default-src 'self'; script-src 'self' 'unsafe-inline' https://d3js.org; style-src 'self' 'unsafe-inline' https://cdn.jsdelivr.net; img-src 'self' data:"
        response.headers['Content-Security-Policy'] = csp

        return response

    # Path Validation Middleware
    @app.before_request
    def validate_file_paths():
        """Validate file paths in requests to prevent path traversal"""
        # Skip for static files
        if request.path.startswith('/static/'):
            return

        # Check for suspicious patterns
        if '..' in request.path or request.path.startswith('/etc/'):
            logger.warning(f"Suspicious path detected: {request.path} from {get_remote_address()}")
            return "Access Denied", 403

        # Validate file parameters
        file_params = ['dataset', 'file', 'fileA', 'fileB', 'path']
        for param in file_params:
            value = request.form.get(param) or request.args.get(param)
            if value and not validate_path(value, ['data', 'mappers', 'templates']):
                logger.warning(f"Path traversal attempt: {param}={value} from {get_remote_address()}")
                return "Invalid file path", 400

    logger.info("Security middleware initialized")


def apply_rate_limits(app):
    """
    Apply specific rate limits to routes

    Args:
        app: Flask application instance
    """
    # Preprocessing endpoint - very limited
    limiter.limit(RATE_LIMIT_PREPROCESS)(app.view_functions.get('preprocess', lambda: None))

    # Data retrieval endpoints - moderate limits
    data_endpoints = ['getDataPair', 'getgenelist', 'preview']
    for endpoint in data_endpoints:
        if endpoint in app.view_functions:
            limiter.limit("50 per minute")(app.view_functions[endpoint])

    logger.info("Route-specific rate limits applied")


def exempt_from_csrf(*routes):
    """
    Exempt specific routes from CSRF protection

    Args:
        *routes: Route names to exempt

    Example:
        exempt_from_csrf('api_endpoint', 'webhook')
    """
    for route in routes:
        csrf.exempt(route)


# Custom decorators for enhanced security

def require_safe_filename(f):
    """Decorator to ensure filenames are sanitized"""
    from functools import wraps
    from utils import sanitize_filename

    @wraps(f)
    def decorated_function(*args, **kwargs):
        # Sanitize any filename parameters
        if request.form:
            for key in request.form:
                if 'filename' in key.lower() or 'file' in key.lower():
                    request.form = request.form.copy()
                    request.form[key] = sanitize_filename(request.form[key])
        return f(*args, **kwargs)
    return decorated_function


def log_security_event(event_type, details):
    """
    Log security events for monitoring

    Args:
        event_type: Type of security event (e.g., 'path_traversal', 'rate_limit')
        details: Details about the event
    """
    logger.warning(f"SECURITY EVENT [{event_type}]: {details} | IP: {get_remote_address()} | Path: {request.path}")
