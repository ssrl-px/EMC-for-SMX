## EMC-for-SMX

This program is designed to reconstruct the Bragg intensities from serial mircocrystallography (SMX) data collected at synchrotron storage ring sources using the [expand-maximize-compress (EMC) algorithm](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.80.026705).

Our program reconstructs the 3D crystal intensity distribution using all the data frames simultaneously without requiring
individual data frames to be able to be oriented by normal indexing methods.

The program was written in *C* and *Python*, with *cbf* as the default input data format.
It is executed on *Linux* using the *MPI* parallelization framework.
The required packages include:

- Requirements for *C*: *gcc*, *OpenMPI* and *OpenSSL*.
- Requirements for *Python*: *NumPy*, *Matplotlib* and *FabIO*.
- *Git*, *X Window System*.

See [the wiki page](https://github.com/tl578/EMC-for-SMX/wiki) for the execution details.

Send comments, bug reports, to: tl578@cornell.edu
