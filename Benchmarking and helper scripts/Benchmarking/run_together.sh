#!/bin/bash
#SBATCH --job-name=mosaicprot
#SBATCH --partition=scc-gpu
#SBATCH --time=2-0:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=250G
#SBATCH --exclusive


#python run_code_module1.py
#python run_code_module2.py
python run_code_candidate_list.py
python run_code_module3.py

echo "Completed"

