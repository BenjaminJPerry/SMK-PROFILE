#!/usr/bin/env python3
# Forked from: https://github.com/rusalkaguy/snakemake-slurm-module/blob/master/slurm-status.py

import re
import subprocess
import shlex
import sys
import time


jobid = sys.argv[1] # original from https://github.com/Snakemake-Profiles/slurm
#jobid = jobid.split()[3] # needed when --parsable not indicated with sbatch

output = str(subprocess.check_output("sacct -j %s --format State --noheader | head -1 | awk '{print $1}'" % jobid, shell=True).strip())

running_status=["PENDING", "CONFIGURING", "COMPLETING", "RUNNING", "SUSPENDED"]
if "COMPLETED" in output:
  print("success")
elif any(r in output for r in running_status):
  print("running")
else:
  print("failed")

