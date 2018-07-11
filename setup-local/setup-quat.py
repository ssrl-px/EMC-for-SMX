from subprocess import *
import numpy as np
import sys
import os

nproc = 10
fp = open(sys.argv[1], "r")
lines = fp.readlines()
fp.close()

length = len(lines)
for t in range(length):
    entry_name = lines[t].split(" = ")[0]
    if ("res_cutoff" == entry_name):
        res_cutoff = np.double(lines[t].split(" = ")[1])
    if ("high_res_cutoff" == entry_name):
        high_res_cutoff = np.double(lines[t].split(" = ")[1])
    if ("quat_file" == entry_name):
        num_coarse_div = np.int(lines[t].split("quaternion")[1].split(".")[0])
    if ("home_dir" == entry_name):
        home_dir = lines[t].split(" = ")[1].strip()
    if ("prob_dir" == entry_name):
        prob_dir = lines[t].split(" = ")[1].strip()

num_fine_div = np.int(np.ceil(res_cutoff*num_coarse_div/high_res_cutoff))
quat_map_file = "quaternion-{0:d}-{1:d}.dat".format(num_coarse_div, num_fine_div)

if ("A" in sys.argv[1]):
    dir_hdr = "low-res-recon-A"
elif ("B" in sys.argv[1]):
    dir_hdr = "low-res-recon-B"
else:
    dir_hdr = "low-res-recon"

    cmd = "gcc ../make-quaternion/make-quaternion.c -O3 -lm -o quat"
    p = Popen(cmd, shell=True)
    p.wait()
    
    cmd = "./quat -bin {0:d}".format(num_coarse_div)
    p = Popen(cmd, shell=True)
    p.wait()
    
    cmd = "./quat -bin {0:d}".format(num_fine_div)
    p = Popen(cmd, shell=True)
    p.wait()
    
    if (os.path.exists(quat_map_file) == False):
        cmd = "mpicc mpi-embed-quaternions.c -O3 -lm -Wall -o embed"
        p = Popen(cmd, shell=True)
        p.wait()
        cmd = "mpirun -np {0:d} ./embed {1:d} {2:d}".format(nproc, num_coarse_div, num_fine_div)
        p = Popen(cmd, shell=True)
        p.wait()

data_id = 0
dirname = "{0:s}-{1:d}".format(dir_hdr, data_id)
data_dir = os.path.join(prob_dir, dirname)
while (os.path.exists(data_dir) == True):
    data_id += 1
    dirname = "{0:s}-{1:d}".format(dir_hdr, data_id)
    data_dir = os.path.join(prob_dir, dirname)

data_id -= 1
dirname = "{0:s}-{1:d}".format(dir_hdr, data_id)
data_dir = os.path.join(prob_dir, dirname)

iter_max = 0
for root, dirs, files in os.walk(data_dir, topdown=True):
    for name in dirs:
        if ("iter_flag" not in name):
            continue
        flag_id = np.int(name.split('-')[1])
        if (flag_id % 2 == 0 and iter_max < flag_id/2):
            iter_max = flag_id/2
            prob_dir = os.path.join(root, name)

aux_dir = os.path.join(home_dir, "aux")

if ("A" in sys.argv[1]):
    start_prob_dir = os.path.join(aux_dir, "start_prob-A")
elif ("B" in sys.argv[1]):
    start_prob_dir = os.path.join(aux_dir, "start_prob-B")
else:
    start_prob_dir = os.path.join(aux_dir, "start_prob")

cmd = "mkdir -p {0:s}".format(start_prob_dir)
p = Popen(cmd, shell=True)
p.wait()

cmd = "cp {0:s}/high_p-*.dat {1:s}".format(prob_dir, start_prob_dir)
p = Popen(cmd, shell=True)
p.wait()

fp = open(sys.argv[1], "w")
for t in range(length):
    if ("start_prob_dir" in lines[t]):
        lines[t] = "start_prob_dir = {0:s}\n".format(start_prob_dir)
    if ("fine_quat_file" in lines[t]):
        if ("A" in sys.argv[1]):
            fine_quat_file = "c-reduced-quaternion{0:d}-A.bin".format(num_fine_div)
        elif ("B" in sys.argv[1]):
            fine_quat_file = "c-reduced-quaternion{0:d}-B.bin".format(num_fine_div)
        else:
            fine_quat_file = "c-reduced-quaternion{0:d}.bin".format(num_fine_div)

        lines[t] = "fine_quat_file = {0:s}\n".format(os.path.join(aux_dir, fine_quat_file))
    if ("quat_table_file" in lines[t]):
        if ("A" in sys.argv[1]):
            quat_table_file = "reduced-quaternion-{0:d}-{1:d}-A.dat".format(num_coarse_div, num_fine_div)
        elif ("B" in sys.argv[1]):
            quat_table_file = "reduced-quaternion-{0:d}-{1:d}-B.dat".format(num_coarse_div, num_fine_div)
        else:
            quat_table_file = "reduced-quaternion-{0:d}-{1:d}.dat".format(num_coarse_div, num_fine_div)
        lines[t] = "quat_table_file = {0:s}\n".format(os.path.join(aux_dir, quat_table_file))
    fp.write(lines[t])
 
fp.close()

cmd = "gcc reduce-quaternions.c -O3 -lm -o reduce"
p = Popen(cmd, shell=True)
p.wait()

cmd = "./reduce {0:s} c-quaternion{1:d}.bin".format(sys.argv[1], num_fine_div)
p = Popen(cmd, shell=True)
p.wait()

cmd = "mv c-reduced-quaternion{0:d}.bin {1:s}/{2:s}".format(num_fine_div, aux_dir, fine_quat_file)
p = Popen(cmd, shell=True)
p.wait()

cmd = "mv reduced-quaternion-{0:d}-{1:d}.dat {2:s}/{3:s}".format(num_coarse_div, num_fine_div, aux_dir, quat_table_file)
p = Popen(cmd, shell=True)
p.wait()
