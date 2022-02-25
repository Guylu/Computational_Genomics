#!/bin/bash
#Write output to JOB_ID.out
#SBATCH -o out-%j.out
#Write error output to JOB_ID.err
#SBATCH -e error-%j.err
#submit to elsc.q partition
#SBATCH --partition=elsc.q
#ask for 50g for the job
#SBATCH --mem=150g
#ask for 8 cpus
#SBATCH -c 1
#mail address
#SBATCH --mail-user=guy.lutsker@mail.huji.ac.il
#send mail for errors or finish
#SBATCH --mail-type=ALL

python -u ex1_rc.py 20


