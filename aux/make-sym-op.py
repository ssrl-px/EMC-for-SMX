import numpy as np

num_op = 16
sym_op = []

sym_op.append([ [1, 0, 0], [0, 1, 0], [0, 0, 1] ])
sym_op.append([ [0, -1, 0], [1, 0, 0], [0, 0, 1] ])
sym_op.append([ [-1, 0, 0], [0, -1, 0], [0, 0, 1] ])
sym_op.append([ [0, 1, 0], [-1, 0, 0], [0, 0, 1] ])
sym_op.append([ [-1, 0, 0], [0, 1, 0], [0, 0, -1] ])
sym_op.append([ [1, 0, 0], [0, -1, 0], [0, 0, -1] ])
sym_op.append([ [0, 1, 0], [1, 0, 0], [0, 0, -1] ])
sym_op.append([ [0, -1, 0], [-1, 0, 0], [0, 0, -1] ])

fp = open("sym-op.dat", "w")
line = "%d\n" % num_op
fp.write(line)
for k in range(int(num_op/2.)):
    for i in range(3):
        for j in range(3):
            line = "%1.3e " % np.double(sym_op[k][i][j])
            fp.write(line)
    fp.write("\n")

for k in range(int(num_op/2.)):
    for i in range(3):
        for j in range(3):
            if (i == j):
                sym_op[k][i][j] *= -1
            line = "%1.3e " % np.double(sym_op[k][i][j])
            fp.write(line)
    fp.write("\n")

fp.close()
