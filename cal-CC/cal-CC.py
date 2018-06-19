import matplotlib.pyplot as plt
import numpy as np
import sys

fp = open("../aux/sym-op.dat", "r")
num_sym_op = np.int(fp.readline().strip())
sym_op = np.zeros((num_sym_op, 9))
for i in range(num_sym_op):
    sym_op[i] = np.copy(np.array(fp.readline().split()).astype(np.double))
fp.close()

fp = open("../aux/basis-vec.dat", "r")
lines = fp.readlines()
fp.close()

basis = np.zeros((3, 3))
for i in range(3):
    basis[i] = np.copy(np.array(lines[i].split()).astype(np.double))

rec_basis = np.zeros((3, 3))
rec_basis[0] = np.copy(np.cross(basis[1], basis[2]))
rec_basis[1] = np.copy(np.cross(basis[2], basis[0]))
rec_basis[2] = np.copy(np.cross(basis[0], basis[1]))

vol = np.abs(np.dot(rec_basis[0], basis[0]))
rec_basis[0] /= vol
rec_basis[1] /= vol
rec_basis[2] /= vol
proj = np.linalg.inv(rec_basis)

out_id = 0
prob_dir = "../Data/high-prob"
dir_A = "{0:s}/high-res-recon-A-{1:d}".format(prob_dir, out_id)
dir_B = "{0:s}/high-res-recon-B-{1:d}".format(prob_dir, out_id)
while (os.path.exists(dir_A) == True and os.path.exists(dir_B) == True):
    out_id += 1
    dir_A = "{0:s}/high-res-recon-A-{1:d}".format(prob_dir, out_id)
    dir_B = "{0:s}/high-res-recon-B-{1:d}".format(prob_dir, out_id)
out_id -= 1

iter_flag = 1
dir_A = "{0:s}/high-res-recon-A-{1:d}/iter_flag-{2:d}".format(prob_dir, out_id, iter_flag)
dir_B = "{0:s}/high-res-recon-B-{1:d}/iter_flag-{2:d}".format(prob_dir, out_id, iter_flag)
while (os.path.exists(dir_A) == True and os.path.exists(dir_B) == True):
    iter_flag += 1
    dir_A = "{0:s}/high-res-recon-A-{1:d}/iter_flag-{2:d}".format(prob_dir, out_id, iter_flag)
    dir_B = "{0:s}/high-res-recon-B-{1:d}/iter_flag-{2:d}".format(prob_dir, out_id, iter_flag)
iter_flag -= 1

intens_file_A = "{0:s}/intens-mean-std.dat".format(dir_A)
intens_file_B = "{0:s}/intens-mean-std.dat".format(dir_B)
fp1 = open(intens_file_A, "r")
fp2 = open(intens_file_B, "r")
linesA = fp1.readlines()
linesB = fp2.readlines()
fp1.close()
fp2.close()

intensA = {}
intensB = {}
hkl_indices = []
num_hkl = len(linesA)
for t in range(num_hkl):
    words = linesA[t].split()
    hkl = "{0:s} {1:s} {2:s}".format(words[1], words[2], words[3])
    intensA[hkl] = [np.double(words[4]), 0]
    words = linesB[t].split()
    hkl = "{0:s} {1:s} {2:s}".format(words[1], words[2], words[3])
    intensB[hkl] = [np.double(words[4]), 0]
    hkl_indices.append([np.int(words[1]), np.int(words[2]), np.int(words[3])])

unique_peaks = []
for t in range(num_hkl):
    ave_intensA = 0.
    ave_intensB = 0.
    hkl_ct = 0
    hval = hkl_indices[t][0]
    kval = hkl_indices[t][1]
    lval = hkl_indices[t][2]
    hkl = [hval, kval, lval]
    rvec = np.zeros(3)
    for i in range(3):
        for j in range(3):
            rvec[i] += rec_basis[i][j] * hkl[j]
    for s in range(num_sym):
        rot_rvec = np.zeros(3)
        for i in range(3):
            for j in range(3):
                rot_rvec[i] += sym_op[s][3*i+j]*rvec[j]
        new_hkl = np.zeros(3)
        for i in range(3):
            for j in range(3):
                new_hkl[i] += proj[i][j]*rot_rvec[j]

        new_hval = np.int(new_hkl[0])
        new_kval = np.int(new_hkl[1])
        new_lval = np.int(new_hkl[2])
        
        hkl_id = "{0:d} {1:d} {2:d}".format(new_hval, new_kval, new_lval)
        if (intensA[hkl_id][1] == 0):
            ave_intensA += intensA[hkl_id][0]
            ave_intensB += intensB[hkl_id][0]
            hkl_ct += 1
            intensA[hkl_id][1] = 1
            intensB[hkl_id][1] = 1
    
    ave_intensA /= hkl_ct
    ave_intensB /= hkl_ct
    qval = np.sqrt(rvec[0]**2 + rvec[1]**2 + rvec[2]**2)
    unique_peaks.append([hval, kval, lval, qval, ave_intensA, ave_intensB])

qval_max = 0.
length = len(unique_peaks)
for t in range(length):
    if (unique_peaks[t][3] > qval_max):
        qval_max = unique_peaks[t][3]

qval_max *= 1.01
len_qval = 20
qval = np.linspace(0, qval_max, len_qval)
inv_dq = 1. / (qval[1] - qval[0])
ave_intensA = np.zeros(len_qval)
ave_intensB = np.zeros(len_qval)
count = np.zeros(len_qval, dtype=np.int)
for t in range(length):
    idx = np.int(np.round(unique_peaks[t][3]*inv_dq - 0.5))
    ave_intensA[idx] += unique_peaks[t][4]
    ave_intensB[idx] += unique_peaks[t][5]
    count[idx] += 1

for i in range(len_qval):
    if (count[i] > 0):
        ave_intensA[i] /= count[i]
        ave_intensB[i] /= count[i]

var_AB = np.zeros(len_qval)
var_AA = np.zeros(len_qval)
var_BB = np.zeros(len_qval)
for t in range(length):
    idx = np.int(np.round(unique_peaks[t][3]*inv_dq - 0.5))
    dA = unique_peaks[t][4] - ave_intensA[idx]
    dB = unique_peaks[t][5] - ave_intensB[idx]
    var_AB[idx] += dA*dB
    var_AA[idx] += dA*dA
    var_BB[idx] += dB*dB

CC_half = np.zeros(len_qval)
CC_true = np.zeros(len_qval)
for i in range(len_qval):
    CC_half[i] = var_AB[i]/np.sqrt(var_AA[i]*var_BB[i])
    CC_true[i] = np.sqrt(2*CC_half[i]/(1+CC_half[i]))

plt.plot(ave_qval, CC_true, 'b')
plt.show()
