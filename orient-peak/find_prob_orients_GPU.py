import numpy as np
import time
import h5py
from scipy.spatial.transform import Rotation
import orient
from mpi4py import MPI
COMM = MPI.COMM_WORLD
O = orient.probable_orients(1)
maxNumQ = 200

RES_CUT = 4  # Angstrom, only use strong spots out to this res
HKL_CUT = 0.1  # fractional HKL, measures how close a spot is to the prediction
MIN_SPOTS_FOR_PROBABLE_ROT = 3  # NOTE: might want to increase this considering DIALS pickle up many more spots than the original EMC code

# aligned unit cell matrix (upper diagonal)
a = np.array([79.1, 0, 0])  # real space a
b = np.array([0, 79.1, 0])  # real space b
c = np.array([0, 0, 38.4])  # real space c
Bmat = np.vstack((a, b, c)).T

# path to reference geometry
h5_refl_file = "ALL_qvecs.h5"


def print0(*args, **kwargs):
    if COMM.rank==0:
        print(*args, **kwargs)


F = h5py.File(h5_refl_file, "r", driver="mpio", comm=COMM)


def main():

    # Load the data produced by make-quaternion code
    num_quat = np.fromfile('../aux/c-quaternion70.bin', np.int32, 1)[0]
    quat_data = np.fromfile('../aux/c-quaternion70.bin', np.float64, offset=4)
    quat_data = quat_data.reshape((num_quat, 5))

    # Convert these quats to rotation matrices using scipy
    rotMats = Rotation.from_quat(quat_data[:, :4]).as_matrix()
    #rotMats = rotMats[:8000000]
    # ALLOCATE / COPY ORIENTATIONS FOR GPU
    print("allocating")
    O.allocate_orientations(0, rotMats.ravel(), maxNumQ)
    print("Done allocating")
    # TODO: Consider making rotMats a mpi window

    # all of the strong spot reflections

    #with h5py.File(h5_refl_file, "r", driver="mpio", comm=COMM ) as F:
    with h5py.File("ALL_rotIDX_%d.h5" % COMM.rank, "w") as OUT:

        Nref_per_xtal = F["num_refls"]
        all_Qvecs = F["qvecs"]
        all_Reso = F["resolution"]

        num_xtals = all_Qvecs.shape[0]
        dt = h5py.vlen_dtype(int)
        idx_dset = OUT.create_dataset("probable_rot_inds", (num_xtals,), dtype=dt)
        xtal_inds = []

        print0("Found data on %d xtals, and loaded %d orientations" % (num_xtals, num_quat))
        for i_f in range(num_xtals):
            if i_f % COMM.size != COMM.rank:
                continue

            tstart = time.time()

            print0("Loading data")
            num_refl = Nref_per_xtal[i_f]
            qvecs = all_Qvecs[i_f,:num_refl]
            reso = all_Reso[i_f, :num_refl]
            qvecs = qvecs[reso < RES_CUT]

            # the fractional HKL, for each rotMat and for each Qxyz (strong spot centroid) is given by
            # Hfrac = rotMat*Bmat*Qxyz
            # We can vectorize this using numpy
            Bqs = np.dot(Bmat, qvecs.T).T
            O.orient_peaks(Bqs.ravel(), HKL_CUT, MIN_SPOTS_FOR_PROBABLE_ROT, False)
            probable_rot_ids = O.get_probable_orients()
            idx_dset[i_f] = probable_rot_ids
            #OUT.create_dataset("probable_rots/shot%d"% i_f, data=probable_rot_ids)


            #probable_rot_ids = [rot_id for rot_id, num_spots in C.items() if num_spots >= MIN_SPOTS_FOR_PROBABLE_ROT]
            print0("Crystal %d / %d has %d probable orientations (took %.4f sec)"
                  % ( i_f+1, num_xtals, len(probable_rot_ids), time.time()-tstart))
            xtal_inds.append(i_f)
        OUT.create_dataset("xtal_inds", data=xtal_inds)


if __name__ == "__main__":
    main()
