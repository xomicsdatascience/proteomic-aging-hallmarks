"""Compatibility entry point for current paper figures."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent / "reanalysis"))
from run_A08 import main
if __name__=='__main__':main()
