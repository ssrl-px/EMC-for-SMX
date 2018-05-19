import numpy as np

num_op = 8
sym_op = []

sym_op.append([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
sym_op.append([[-1, 0, 0], [0, -1, 0], [0, 0, 1]])
sym_op.append([[1, 0, 0], [0, -1, 0], [0, 0, -1]])
sym_op.append([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])

fp = open("sym-op.dat", "w")
line = "%d\n" % num_op
fp.write(line)
for k in xrange(num_op/2):
    for i in xrange(3):
        for j in xrange(3):
            line = "%1.3e " % np.double(sym_op[k][i][j])
            fp.write(line)
    fp.write("\n")

for k in xrange(num_op/2):
    for i in xrange(3):
        for j in xrange(3):
            if (i == j):
                sym_op[k][i][j] *= -1
            line = "%1.3e " % np.double(sym_op[k][i][j])
            fp.write(line)
    fp.write("\n")

fp.close()
