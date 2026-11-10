"""Shared defaults; search sensitivity never changes the result display limit."""
DEFAULT_HITS = 100
SEARCH_PRESETS = {
    'standard': {'evalue': 1e-5, 'sensitivity': 'sensitive'},
    'sensible': {'evalue': 1e-2, 'sensitivity': 'very-sensitive'},
    'strict': {'evalue': 1e-10, 'sensitivity': 'sensitive'},
}
