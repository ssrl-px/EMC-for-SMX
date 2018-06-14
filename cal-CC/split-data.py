from subprocess import *
import numpy as np
import sys
import os

fp = open("../config.ini", "r")
lines = fp.readlines()
fp.close()

linesA = []
linesB = []

length = len(lines)
for t in range(length):
    linesA.append(lines[t])
    linesB.append(lines[t])
    entry_name = lines[t].split(" = ")[0]
    if (" = " in lines[t]):
        file_name = lines[t].split(" = ")[1].strip()
    if ("mpi_bgfile" == entry_name):
        linesA[t] = "{0:s} = {1:s}-A.bin\n".format(entry_name, file_name.split('.')[0])
        linesB[t] = "{0:s} = {1:s}-B.bin\n".format(entry_name, file_name.split('.')[0])
    if ("mpi_datafile" == entry_name):
        linesA[t] = "{0:s} = {1:s}-A.bin\n".format(entry_name, file_name.split('.')[0])
        linesB[t] = "{0:s} = {1:s}-B.bin\n".format(entry_name, file_name.split('.')[0])
    if ("prob_orien_file" == entry_name):
        linesA[t] = "{0:s} = {1:s}-A.bin\n".format(entry_name, file_name.split('.')[0])
        linesB[t] = "{0:s} = {1:s}-B.bin\n".format(entry_name, file_name.split('.')[0])
    if ("reduced_data_id_file" == entry_name):
        linesA[t] = "{0:s} = {1:s}-A.dat\n".format(entry_name, file_name.split('.')[0])
        linesB[t] = "{0:s} = {1:s}-B.dat\n".format(entry_name, file_name.split('.')[0])
        reduced_data_id_file = file_name
    if ("start_phi_file" == entry_name):
        linesA[t] = "{0:s} = {1:s}-A.dat\n".format(entry_name, file_name.split('.')[0])
        linesB[t] = "{0:s} = {1:s}-B.dat\n".format(entry_name, file_name.split('.')[0])
    if ("start_intens_file" == entry_name):
        linesA[t] = "{0:s} = {1:s}-A.bin\n".format(entry_name, file_name.split('.')[0])
        linesB[t] = "{0:s} = {1:s}-B.bin\n".format(entry_name, file_name.split('.')[0])
    if ("start_prob_dir" == entry_name):
        linesA[t] = "{0:s} = {1:s}-A\n".format(entry_name, file_name.split('.')[0])
        linesB[t] = "{0:s} = {1:s}-B\n".format(entry_name, file_name.split('.')[0])
    if ("local_r2peak_file" == entry_name):
        linesA[t] = "{0:s} = {1:s}-A.bin\n".format(entry_name, file_name.split('.')[0])
        linesB[t] = "{0:s} = {1:s}-B.bin\n".format(entry_name, file_name.split('.')[0])
    if ("local_peak2r_file" == entry_name):
        linesA[t] = "{0:s} = {1:s}-A.bin\n".format(entry_name, file_name.split('.')[0])
        linesB[t] = "{0:s} = {1:s}-B.bin\n".format(entry_name, file_name.split('.')[0])

fp = open(reduced_data_id_file, "r")
lines = fp.readlines()
fp.close()

num_data = np.int(lines[0])
num_dataB = num_data/2
num_dataA = num_data - num_dataB

fp1 = open("../reduce-data/reduced-data_id-A.dat", "w")
fp2 = open("../reduce-data/reduced-data_id-B.dat", "w")
tmp = "{0:d}\n".format(num_dataA)
fp1.write(tmp)
tmp = "{0:d}\n".format(num_dataB)
fp2.write(tmp)

for i in xrange(num_data):
    if (i % 2 == 0):
        fp1.write(lines[i+1])
    else:
        fp2.write(lines[i+1])
fp1.close()
fp2.close()

fp1 = open("../config-A.ini", "w")
fp2 = open("../config-B.ini", "w")
for t in xrange(length):
    fp1.write(linesA[t])
    fp2.write(linesB[t])
fp1.close()
fp2.close()
