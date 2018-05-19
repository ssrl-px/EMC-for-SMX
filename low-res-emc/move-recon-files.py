from subprocess import *
import numpy as np
import os

fp = open("../config.ini", "r")
lines = fp.readlines()
fp.close()

length = len(lines)
for t in range(length):
    if ("nproc" in lines[t]):
        nproc = np.int(lines[t].split(" = ")[1])
    if ("prob_dir" == lines[t].split(" = ")[0]):
        prob_dir = lines[t].split(" = ")[1].strip()
    if ("start_phi_file" in lines[t]):
        start_phi_file = lines[t].split(" = ")[1].strip()
    if ("start_intens_file" in lines[t]):
        start_intens_file = lines[t].split(" = ")[1].strip()

dir_list = []
for root, dirs, files in os.walk(prob_dir, topdown=True):
    for name in dirs:
        if ("iter_flag" not in name):
            continue
        dir_name = os.path.join(root, name)
        if (prob_dir != dir_name.split('/iter_flag')[-2]):
            continue
        dir_list.append(dir_name)

out_id = 0
dirname = "low-res-recon-{0:d}".format(out_id)
out_dir = os.path.join(prob_dir, dirname)
while (os.path.exists(out_dir) == True):
    out_id += 1
    dirname = "low-res-recon-{0:d}".format(out_id)
    out_dir = os.path.join(prob_dir, dirname)

cmd = "mkdir -p {0:s}".format(out_dir)
p = Popen(cmd, shell=True)
p.wait()

for dir_name in dir_list:
    cmd = "mv {0:s} {1:s}".format(dir_name, out_dir)
    p = Popen(cmd, shell=True)
    p.wait()

cmd = "mv {0:s} {1:s}".format(start_phi_file, out_dir)
p = Popen(cmd, shell=True)
p.wait()

cmd = "mv {0:s} {1:s}".format(start_intens_file, out_dir)
p = Popen(cmd, shell=True)
p.wait()

cmd = "mv run.log {0:s}".format(out_dir)
p = Popen(cmd, shell=True)
p.wait()

cmd = "mv EMC.log {0:s}".format(out_dir)
p = Popen(cmd, shell=True)
p.wait()
