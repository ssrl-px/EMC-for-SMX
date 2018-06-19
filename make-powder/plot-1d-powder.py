import matplotlib.pyplot as plt
import numpy as np
import sys

fin = open("1d-pseudo-powder.dat", "r")
lines = fin.readlines()
fin.close()

length = len(lines)
qval = np.zeros(length)
radial_hist = np.zeros((2, length))
for i in range(length):
    words = lines[i].split()
    qval[i] = np.double(words[0])
    radial_hist[0][i] = np.double(words[1])
    radial_hist[1][i] = np.double(words[2])

a = 79.1
b = 79.1
c = 38.4
idx_max = 10
qmax = 0.075
epsilon = 1.e-10

# reciprocal basis vectors
vec_a = np.array([1./a, 0., 0.])
vec_b = np.array([0., 1./b, 0.])
vec_c = np.array([0., 0., 1./c])
rvec = np.zeros(3)

for p in range(2):

    if (p == 0):
        print("1D histogram of inter-peak distances in reciprocal space")
    else:
        print("1D histogram of spatial frequency magnitudes of peaks")

    ymax = np.max(radial_hist[p])*1.1
    fig, ax = plt.subplots(figsize=(15, 6))
    plt.plot(qval, radial_hist[p], 'b')
    
    for h in range(-idx_max, idx_max+1):
        for k in range(-idx_max, idx_max+1):
            for l in range(-idx_max, idx_max+1):
                q = 0.
                for i in range(3):
                    rvec[i] = h*vec_a[i] + k*vec_b[i] + l*vec_c[i]
                    q += rvec[i]*rvec[i]
                q = np.sqrt(q)
                if (q > qmax):
                    continue
    
                y = np.arange(0, ymax, 100)
                x = np.ones_like(y)*q
                plt.plot(x, y, 'r--')
                
    plt.xlim([0, qmax])
    plt.ylim([0, ymax])
    plt.show()
