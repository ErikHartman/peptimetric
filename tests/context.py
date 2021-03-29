import sys
from pathlib import Path


module_= Path(__file__).parent.parent / "src/pept-compare/"
sys.path.insert(0, str(module_path.absolute()))


# TODO: import other modules and make more test_module.py files
import methods
from playground import merge_dataframes
