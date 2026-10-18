import sys
from pathlib import Path

# Make the repository root importable when pytest is run from a checkout
# without installing the package.
sys.path.insert(0, str(Path(__file__).resolve().parent))
