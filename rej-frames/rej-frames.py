from subprocess import *
import numpy as np
import os

fp = open(sys.argv[1], "r")
lines = fp.readlines()
fp.close()

length = len(lines)
for t in range(length):
    if ("prob_dir" == lines[t].split(" = ")[0]):
        prob_dir = lines[t].split(" = ")[1].strip()
    if ("start_phi_file" in lines[t]):
        start_phi_file = lines[t].split(" = ")[1].strip()

if ("A" in sys.argv[1]):
    dir_hdr = "low-res-recon-A"
elif ("B" in sys.argv[1]):
    dir_hdr = "low-res-recon-B"
else:
    dir_hdr = "low-res-recon"

data_id = 0
dirname = "{0:s}-{1:d}".format(dir_hdr, out_id)
data_dir = os.path.join(prob_dir, dirname)
while (os.path.exists(data_dir) == True):
    data_id += 1
    dirname = "{0:s}-{1:d}".format(dir_hdr, out_id)
    data_dir = os.path.join(prob_dir, dirname)

data_id -= 1
dirname = "{0:s}-{1:d}".format(dir_hdr, out_id)
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

cmd = "gcc reject-frames.c -O3 -lm -o rej"
p = Popen(cmd, shell=True)
p.wait()

cmd = "./rej {0:s} {1:s} > run.log".format(sys.argv[1], prob_dir)
p = Popen(cmd, shell=True)
p.wait()

cmd = "mv out-phi.dat {0:s}".format(start_phi_file)
p = Popen(cmd, shell=True)
p.wait()
