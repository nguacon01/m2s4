#!/bin/bash

#SBATCH -p fast ###la partition sur laquelle on soumet le job
#SBATCH -n 2 ##le nombre de coeurs alloués
#SBATCH --mem-per-cpu=2000 ###la mémoire requise par CPU (pas sûr que ce paramètre soit réllement utile
#SBATCH -o /home/mddo/stage/M2S4/logs/slurmPools.%N.%j.out   ###position du log d'output
#SBATCH -e /home/mddo/stage/M2S4/logs/slurmPools.%N.%j.err   ###position du log d'erreur

python /home/mddo/stage/M2S4/code/main.py

