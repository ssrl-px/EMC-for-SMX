from subprocess import *
import sys
import os

num_div = int(sys.argv[1])

cmd = "gcc make-quaternion.c -O3 -lm -o quat"
p = Popen(cmd, shell=True)
p.wait()

cmd = "./quat -bin {0:d}".format(num_div)
p = Popen(cmd, shell=True)
p.wait()

filename = "c-quaternion{0:d}.bin".format(num_div)
cmd = "mv {0:s} ../aux".format(filename)
p = Popen(cmd, shell=True)
p.wait()

filepath = "{0:s}/aux/{1:s}".format(os.path.abspath("../"), filename)
fp = open("../config.ini", "r")
lines = fp.readlines()
fp.close()

for t in range(len(lines)):
    if ("quat_file" == lines[t].split(" = ")[0]):
        lines[t] = "quat_file = {0:s}\n".format(filepath)

fp = open("../config.ini", "w")
for t in range(len(lines)):
    fp.write(lines[t])
fp.close()
