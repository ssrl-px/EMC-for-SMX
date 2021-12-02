"""
make filelists

usage:
python make-filelists.py [glob pointing to raw data cbfs (surrounded in quotes)] > run.log &

e.g.

python make-filelists.py "/path//data/data[1-8]/*cbf"

makes:
cbflist.dat, radiallist.dat, orienlist.dat
outlierlist.dat, peaklist.dat

"""

import numpy as np
import sys
import os
import glob

raw_cbffiles = glob.glob(sys.argv[1])

print("Found %d cbf files in glob" % len(raw_cbffiles))
new_data_dir = "{0:s}/Data".format(os.path.abspath("../"))

dirnames = set([os.path.dirname(f) for f in raw_cbffiles])

outfolders = "radial-bg", "outlier"
for d in dirnames:
    new_d = os.path.join(new_data_dir, os.path.basename(d))
    for d in outfolders:
        sub_d = os.path.join(new_d, d)
        if not os.path.exists(sub_d):
            os.makedirs(sub_d)
        print(sub_d)


new_file_templates = {"ave_bg-%s.bin": (outfolders[0], [], "radiallist.dat"),
                      "prob_orien-%s.dat": (outfolders[0], [], "orienlist.dat"),
                      "outlier-%s.dat": (outfolders[1], [], "outlierlist.dat"),
                      "peak-%s.dat": (outfolders[1], [], "peaklist.dat")}

for i_f, f in enumerate(raw_cbffiles):
    raw_dirname = os.path.dirname(f)
    new_dirname = os.path.join(new_data_dir, os.path.basename(raw_dirname))

    tag = os.path.basename(f).split("_")[-1].split(".cbf")[0]
    for f_template, (subdir, store, _) in new_file_templates.items():
        new_f = os.path.join(new_dirname, subdir, f_template % tag)
        store.append(new_f)
    print(i_f, len(raw_cbffiles))


np.savetxt("cbflist.dat", raw_cbffiles, "%s")
for _, store, outname in new_file_templates.values():
    np.savetxt(outname, store, "%s")
