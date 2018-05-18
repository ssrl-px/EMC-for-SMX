import matplotlib.pyplot as plt
import numpy as np

fp = open("2d-pseudo-powder.dat", "r")
det = np.array(fp.readlines()).astype(np.int)
fp.close()

fp = open("../config.ini", "r")
lines = fp.readlines()
fp.close()
for line in lines:
    if "num_row" in line:
        num_row = np.int(line.split(" = ")[1])
    if "num_col" in line:
        num_col = np.int(line.split(" = ")[1])

det = np.reshape(det, (num_row, num_col))
img = np.log10(det+1.0)
fig, ax = plt.subplots()
im = ax.imshow(img, interpolation='nearest')
plt.show()
