"""Run the frozen analysis in dependency order, stopping on any error."""
import argparse
import os
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--figures-only',action='store_true',help='Regenerate figures, prose and Word exports from existing A01-A07 tables')
    args=parser.parse_args()
    steps=[] if args.figures_only else ['run_A01_A03.py','run_A02.py','run_A04.py','run_A05.py','run_A06.py','run_A07.py']
    steps+=['run_A08.py','build_manuscript.py','A08/build_docx.py','validate_outputs.py']
    env=os.environ.copy();env['OPENBLAS_NUM_THREADS']='1';env['OMP_NUM_THREADS']='1'
    for step in steps:
        print('Running '+step,flush=True)
        subprocess.run([sys.executable,str(ROOT/'reanalysis'/step)],cwd=ROOT,env=env,check=True)

if __name__=='__main__':main()
