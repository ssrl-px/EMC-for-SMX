import numpy as np

fp = open("reduced-peakfiles.dat", "r")
infiles = fp.readlines()
fp.close()

fp = open("../config.ini", "r")
for line in fp:
    if ("max_patch_sz" in line):
        max_patch_sz = np.int(line.split("=")[1].strip())
fp.close()

num_data = len(infiles)
phi = np.zeros(num_data)
phi_sum = 0.
for d in xrange(num_data):
    infile = infiles[d].strip()
    fp = open(infile, "r")
    content = fp.readlines()
    fp.close()

    num_peak = len(content)/2
    num_diffuse_peak = 0
    for i in xrange(num_peak):
        patch_sz = int(content[2*i].strip())
        if (patch_sz > max_patch_sz):
            num_diffuse_peak += 1
            continue
        words = content[2*i+1].split()
        for j in xrange(patch_sz):
            phi[d] += np.double(words[2*j+1])
    phi[d] /= (num_peak - num_diffuse_peak)
    phi_sum += phi[d]

rescale = num_data/phi_sum
fout = open("start-phi.dat", "w")
for d in xrange(num_data):
    line = "%1.15e\n" % (phi[d]*rescale)
    fout.write(line)
fout.close()

fa = open("start-phi-A.dat", "w")
fb = open("start-phi-B.dat", "w")
for d in xrange(num_data):
    line = "%1.15e\n" % (phi[d]*rescale)
    if (d % 2 == 0):
        fa.write(line)
    else:
        fb.write(line)
fa.close()
fb.close()
