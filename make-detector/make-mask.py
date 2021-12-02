"""
make detector mask

usage:
python make-mask.py [path to a data frame in cbf format] > run.log &

need:
config.ini

makes:
mask.dat

"""

import numpy as np
import matplotlib.pyplot as plt
import fabio
import sys

short_max = 2**15 - 1
fp = open("../config.ini", "r")
lines = fp.readlines()
fp.close()

print("run make-mask.py:\n")
for line in lines:
    words = line.split(" = ")
    if (words[0] == "num_row"):
        num_row = np.int(words[1].strip())
    if (words[0] == "num_col"):
        num_col = np.int(words[1].strip())
    if (words[0] == "cx"):
        c_row = np.double(words[1].strip())
    if (words[0] == "cy"):
        c_col = np.double(words[1].strip())

cbf_file = sys.argv[1]
print("cbf_file = {0:s}".format(cbf_file))
print("num_row = {0:d}, num_col = {1:d}, c_row = {2:.1f}, \
c_col = {3:.1f}\n".format(num_row, num_col, c_row, c_col))

img = fabio.open(cbf_file)
img = np.array(img.data)
mask = np.zeros_like(img, dtype=np.int)
# mask detector gaps
mask[img < 0] = 1
# mask hot pixels
mask[img > short_max] = 1

fp = open("../aux/mask.dat", "w")
MASK = []
for i in range(num_row):
    for j in range(num_col):
        x = i - c_row
        y = j - c_col
        # mask shadow of beamstop holder
        if (y < 0 and np.abs(x) < 110):
            mask[i][j] = 1
        tmp = "{0:d}\n".format(mask[i][j])
        MASK.append(mask[i][j])
        fp.write(tmp)
fp.close()
MASK = np.array(MASK)
MASK.tofile("MASK.bin", format="%d") 
img[mask == 1] = 0.
img = np.log10(img + 1.)
fig, ax = plt.subplots()
im = ax.imshow(img, interpolation='nearest')
cb = fig.colorbar(im, ax=ax)
plt.show()
