#!/usr/bin/env python3
# Forked from: https://github.com/rusalkaguy/snakemake-slurm-module/blob/master/slurm-status.py

import re
import subprocess
import shlex
import sys
import time


jobid = sys.argv[1] # original from https://github.com/Snakemake-Profiles/slurm
jobid = jobid.split()[3]

status = str(subprocess.check_output("sacct -j %s --format State --noheader | head -1 | awk '{print $1}'" % jobid, shell=True).strip())

if status == "BOOT_FAIL":
    print("failed")
elif status == "OUT_OF_MEMORY":
    print("failed")
elif status.startswith("CANCELLED"):
    print("failed")
elif status == "COMPLETED":
    print("success")
elif status == "DEADLINE":
    print("failed")
elif status == "FAILED":
    print("failed")
elif status == "NODE_FAIL":
    print("failed")
elif status == "PREEMPTED":
    print("failed")
elif status == "TIMEOUT":
    print("failed")
# Unclear whether SUSPENDED should be treated as running or failed
elif status == "SUSPENDED":
    print("failed")
else:
    print("running")

