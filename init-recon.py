from subprocess import *
import sys
import os

cur_dir = os.getcwd()
data_dir = "{0:s}/Data".format(sys.argv[1])
prob_dir = "{0:s}/high-prob".format(data_dir)
if not os.path.exists(prob_dir):
    os.makedirs(prob_dir)

cmd = "ln -s {0:s} Data".format(data_dir)
p = Popen(cmd, shell=True)
p.wait()

fp = open("config.ini", "r")
lines = fp.readlines()
fp.close()

length = len(lines)
for t in range(length):
    entry_name = lines[t].split(" = ")[0]
    if ("home_dir" == entry_name):
        lines[t] = "home_dir = {0:s}\n".format(cur_dir)
    if ("prob_dir" == entry_name):
        lines[t] = "prob_dir = {0:s}\n".format(prob_dir)
    if ("quat_file" == entry_name):
        lines[t] = "quat_file = {0:s}/aux/c-quaternion70.bin\n".format(cur_dir)
    if ("start_phi_file" == entry_name):
        lines[t] = "start_phi_file = {0:s}/aux/start-phi.dat\n".format(cur_dir)
    if ("start_intens_file" == entry_name):
        lines[t] = "start_intens_file = {0:s}/aux/start_intensity.bin\n".format(cur_dir)
    if ("reduced_data_id_file" == entry_name):
        lines[t] = "reduced_data_id_file = {0:s}/reduce-data/reduced-data_id.dat\n".format(cur_dir)
    if ("prob_orien_file" == entry_name):
        lines[t] = "prob_orien_file = {0:s}/aux/prob-orien.bin\n".format(cur_dir)
    if ("mpi_bgfile" == entry_name):
        lines[t] = "mpi_bgfile = {0:s}/mpi-bg_model.bin\n".format(data_dir)
    if ("mpi_datafile" == entry_name):
        lines[t] = "mpi_datafile = {0:s}/mpi-datafile.bin\n".format(data_dir)
    if ("r2peak_file" == entry_name):
        lines[t] = "r2peak_file = {0:s}/mpi-r2peak-low.bin\n".format(data_dir)
    if ("peak2r_file" == entry_name):
        lines[t] = "peak2r_file = {0:s}/mpi-peak2r-low.bin\n".format(data_dir)
    
fp = open("config.ini", "w")
for line in lines:
    fp.write(line)
fp.close()
