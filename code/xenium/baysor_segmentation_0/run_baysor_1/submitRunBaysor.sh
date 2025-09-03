#!/bin/bash
#BSUB -J submitrun_baysor
#BSUB -n 10
#BSUB -q cpuqueue
#BSUB -W 01:00
#BSUB -R "rusage[mem=15G]"
#BSUB -o /lila/data/deyk/harry/spatial/RA_Xenium/job_logs/submit_run_baysor.log
#BSUB -e /lila/data/deyk/harry/spatial/RA_Xenium/job_logs/submit_run_baysor.log

python /lila/data/deyk/harry/spatial/RA_Xenium/job_scripts/runBaysor/runBaysor.py