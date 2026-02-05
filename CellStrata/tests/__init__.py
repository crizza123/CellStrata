[pytest]
# Test discovery
testpaths = tests
python_files = test_*.py
python_classes = Test*
python_functions = test_*

# Markers
markers =
    integration: marks tests as integration tests (require network access)
    slow: marks tests as slow running

# Output settings
addopts = -v --tb=short

# Ignore warnings from dependencies
filterwarnings =
    ignore::DeprecationWarning
    ignore::PendingDeprecationWarning