#!/bin/bash
#Write output to JOB_ID.out
#SBATCH -o out-%j.out
#Write error output to JOB_ID.err
#SBATCH -e error-%j.err
#submit to elsc.q partition
#SBATCH --partition=elsc.q
#ask for 50g for the job
#SBATCH --mem=100g
#ask for 8 cpus
#SBATCH -c 16
#mail address
#SBATCH --mail-user=guy.lutsker@mail.huji.ac.il
#send mail for errors or finish
#SBATCH --mail-type=ALL
CUDA_LAUNCH_BLOCKING=1 python -u ex1.py
