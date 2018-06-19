"""
make filelists

usage:
python make-filelists.py [path to dir that stores raw data] > run.log &

makes:
cbflist.dat, radiallist.dat, orienlist.dat
outlierlist.dat, peaklist.dat

"""

from subprocess import *
import sys
import os

raw_data_dir = os.path.abspath(sys.argv[1])
new_data_dir = "{0:s}/Data".format(os.path.abspath("../"))

dir_list = []
for root, dirs, files in os.walk(raw_data_dir, topdown=True):
    for name in dirs:
        dir_name = os.path.join(root, name)
        contents = os.listdir(dir_name)
        if_cbf = 0
        for f in contents:
            if (f[-4:] == ".cbf"):
                if_cbf = 1
        if (if_cbf == 0):
            continue
        dir_list.append(dir_name)

dir_list = sorted(dir_list)

cbf_files = []
radial_files = []
orien_files = []
outlier_files = []
peak_files = []

for src_dir in dir_list:
    dir_name = src_dir.split(raw_data_dir)[1]
    radial_dir = "{0:s}{1:s}/radial-bg".format(new_data_dir, dir_name)
    cmd = "mkdir -p {0:s}".format(radial_dir)
    p = Popen(cmd, shell=True)
    p.wait()
    
    outlier_dir = "{0:s}{1:s}/outlier".format(new_data_dir, dir_name)
    cmd = "mkdir -p {0:s}".format(outlier_dir)
    p = Popen(cmd, shell=True)
    p.wait()
    
    files = sorted(os.listdir(src_dir))
    for f in files:
        if (f[-4:] == ".cbf"):
            cur_str = "{0:s}/{1:s}\n".format(src_dir, f)
            cbf_files.append(cur_str)

            fileid = f.split('_')[-1].split('.')[0]
            cur_str = "{0:s}/ave_bg-{1:s}.bin\n".format(radial_dir, fileid)
            radial_files.append(cur_str)

            cur_str = "{0:s}/prob_orien-{1:s}.dat\n".format(radial_dir, fileid)
            orien_files.append(cur_str)

            cur_str = "{0:s}/outlier-{1:s}.dat\n".format(outlier_dir, fileid)
            outlier_files.append(cur_str)

            cur_str = "{0:s}/peak-{1:s}.dat\n".format(outlier_dir, fileid)
            peak_files.append(cur_str)

fp = open("cbflist.dat", "w")
for line in cbf_files:
    fp.write(line)
fp.close()

fp = open("radiallist.dat", "w")
for line in radial_files:
    fp.write(line)
fp.close()

fp = open("orienlist.dat", "w")
for line in orien_files:
    fp.write(line)
fp.close()

fp = open("outlierlist.dat", "w")
for line in outlier_files:
    fp.write(line)
fp.close()

fp = open("peaklist.dat", "w")
for line in peak_files:
    fp.write(line)
fp.close()

num_raw_data = len(cbf_files)
fp = open("../config.ini", "r")
lines = fp.readlines()
fp.close()

for t in range(len(lines)):
    if ("num_raw_data" in lines[t]):
        lines[t] = "num_raw_data = {0:d}\n".format(num_raw_data)

fp = open("../config.ini", "w")
for t in range(len(lines)):
    fp.write(lines[t])
fp.close()
