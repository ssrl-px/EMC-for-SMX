from subprocess import *
import numpy as np
import sys
import os

fp = open(sys.argv[1], "r")
lines = fp.readlines()
fp.close()

length = len(lines)
for t in range(length):
    if ("home_dir" in lines[t]):
        home_dir = lines[t].split(" = ")[1].strip()
    if ("nproc" in lines[t]):
        nproc = np.int(lines[t].split(" = ")[1])

data_dir = os.path.join(home_dir, "Data")

if ("A" in sys.argv[1]):
    local_r2peak_file = os.path.join(data_dir, "mpi-r2peak-local-A.bin")
    local_peak2r_file = os.path.join(data_dir, "mpi-peak2r-local-A.bin")
elif ("B" in sys.argv[1]):
    local_r2peak_file = os.path.join(data_dir, "mpi-r2peak-local-B.bin")
    local_peak2r_file = os.path.join(data_dir, "mpi-peak2r-local-B.bin")
else:
    local_r2peak_file = os.path.join(data_dir, "mpi-r2peak-local.bin")
    local_peak2r_file = os.path.join(data_dir, "mpi-peak2r-local.bin")

fp = open(sys.argv[1], "w")
for t in range(length):
    if ("local_r2peak_file" in lines[t]):
        lines[t] = "local_r2peak_file = {0:s}\n".format(local_r2peak_file)
    if ("local_peak2r_file" in lines[t]):
        lines[t] = "local_peak2r_file = {0:s}\n".format(local_peak2r_file)
    fp.write(lines[t])
fp.close()

cmd = "mpicc ../make-Ematrix/mpi-make-Emat.c -O3 -lm -o emat"
p = Popen(cmd, shell=True)
p.wait()

cmd = "mpirun -np {0:d} ./emat {1:s}".format(nproc, sys.argv[1])
p = Popen(cmd, shell=True)
p.wait()
