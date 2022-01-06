import numpy as np
import time
from scipy.spatial.transform import Rotation
import glob
from dials.array_family import flex
from mpi4py import MPI
COMM = MPI.COMM_WORLD
from dxtbx.model import ExperimentList
from xfel.merging.application.utils.memory_usage import get_memory_usage
from collections import Counter


def print0(*args, **kwargs):
    if COMM.rank==0:
        print(*args, **kwargs)

RES_CUT = 4
HKL_CUT = 0.2
MIN_SPOTS_FOR_PROBABLE_ROT = 3  # NOTE: might want to increase this considering DIALS pickle up many more spots than the original EMC code

refGeo = "../make-detector/geom_+x_-y.expt"
El = ExperimentList.from_file(refGeo)
detector = El[0].detector
beam = El[0].beam

# Load the data produced by make-quaternion code
num_quat = np.fromfile('../aux/c-quaternion70.bin', np.int32, 1)[0]
quat_data = np.fromfile('../aux/c-quaternion70.bin', np.float64, offset=4)
quat_data = quat_data.reshape((num_quat, 5))

# Convert these quats to rotation matrices using scipy
rotMats = Rotation.from_quat(quat_data[:, :4]).as_matrix()
# TODO: Consider making rotMats a mpi window

a = np.array([79.1, 0, 0])
b = np.array([0, 79.1, 0])
c = np.array([0, 0, 38.4])

Bmat = np.vstack((a, b, c)).T

# all of the strong spot reflections

fnames = glob.glob("../../all_data[1-8]_spots4/*strong.refl")
num_xtals = len(fnames)
print0("Found data on %d xtals, and loaded %d orientations" % (num_xtals, num_quat))

# chunk over the rotMats to save memory
#num_chunks = 10
#chunk_inds = np.array_split(np.arange(num_quat), num_chunks)
# NOTE: sacrifice the chunk code for readability purposes

for i_f, f in enumerate(fnames):
    if i_f % COMM.size != COMM.rank:
        continue

    tstart = time.time()
    print0("Loading refl table from %s" % f)
    R = flex.reflection_table.from_file(f)

    # add in the rlp (Q-vector) column
    print0("Converting strong spot centroids to rlp (inverse Angstrom)")
    R.centroid_px_to_mm(El)
    R.map_centroids_to_reciprocal_space(El)

    # resolution of each refl; only use refls within the resolution cutoff
    print0("Applying resolution cut (%f Angstrom)" % RES_CUT)
    R['res'] = flex.double(1 / np.linalg.norm(R['rlp'].as_numpy_array(), axis=1))
    nR = len(R)
    R = R.select(R['res'] < RES_CUT)
    print0("-- %d / %d refls remain after res cut" % (len(R), nR))

    # the fractional HKL, for each rotMat and for each Qxyz (strong spot centroid) is given by
    # Hfrac = rotMat*Bmat*Qxyz
    # We can vectorize this using numpy
    qvecs = R['rlp'].as_numpy_array()
    possible_rot_ids = []
    for i_q, q in enumerate(qvecs):
        print0("Iterating over rlp %d / %d" % (i_q+1,qvecs.shape[0]))
        Bq = np.dot(Bmat, q)
        #for inds in chunk_inds:
        #    istart, istop = inds[0], inds[-1]+1
        #    Hfracs = np.dot(rotMats[istart: istop], Bq)
        #    hkl = np.ceil(Hfracs - 0.5)
        Hfracs = np.dot(rotMats, Bq)
        Hkl = np.ceil(Hfracs - 0.5)
        Hkl_dist = np.linalg.norm(Hfracs - Hkl , axis=1)
        within_hkl = Hkl_dist < HKL_CUT
        possible_rot_ids.append(np.where(within_hkl)[0])

    # combined them all and check for duplicates. If a rotation was flagged for more then 3 refls, then we consider
    # is a probably rotation
    print0("Finding probable orientations...")
    possible_rot_ids = np.hstack(possible_rot_ids)
    C = Counter(possible_rot_ids)  # creates dictionary of rot_id -> number of occurences
    probable_rot_ids = [rot_id for rot_id, num_spots in C.items() if num_spots >= MIN_SPOTS_FOR_PROBABLE_ROT]
    print0("Crystal corresponding to %s (%d / %d) has %d probable orientations (took %.4f sec)"
          % (f, i_f, num_xtals, len(probable_rot_ids), time.time()-tstart))


