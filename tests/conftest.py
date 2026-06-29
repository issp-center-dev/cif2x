import os
import sys

# When CIF2X_TEST_INSTALLED is set (by the CI workflow, after `pip install .`),
# test the installed distribution. Otherwise always use the in-repo src/ tree so
# the tests exercise the current checkout regardless of any globally-installed
# cif2x.
if not os.environ.get("CIF2X_TEST_INSTALLED"):
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
