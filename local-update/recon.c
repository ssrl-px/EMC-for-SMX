/*
high-resolution intensity reconstruction using 
the local update scheme of the EMC algorithm

usage:
mpirun -np [number of processors] ./emc [path of the config_file] [number of iterations]

needs:
config_file, mpi_datafile, mpi_bgfile, r2peak_file, peak2r_file,
rec-vectors.bin, rec2pix-map.dat, pix-map.dat, basis-vec.dat,
quat_file, quat_table_file, start-phi.dat, start_intensity.bin

makes:
probability distribution, finish_intensity.bin, hkl-val.dat

*/

#include "emc.h"

const long long int_max = 2147483647 ;
const double golden_ratio = 1.618034 ;
const double CP_THRES = 0.999 ;
const double P_MIN = 0.0001 ;
const int mem_step_sz = 128 ;

char config_file[256], data_file[256], peak2r_file[256] ;
char start_intens_file[256], home_dir[256], prob_dir[256] ;
double **ave_bg, **log_bg, *phi, **intens1, *logP_offset, **pval_d2r ;
double rlatt_vec[3][3], proj[3][3], rot[3][3], logP_MIN ;
double qval_max, rmin2, rmax2, wl_D, gw, min_rcell_sz ;
int *r_block_id, *r_block_sz_mpi, **r_block_id_mpi ;
int *rid_map, *pix_map, **my_pval_id, *my_pval_ct ;
int **quat_table_idx, *num_quat_neighbor, *tomo_idx ;
int **most_probable_quat_id, *most_probable_num_p ;
int *num_r2hkl, **r2hkl_mpi, **r2hkl, *num_hkl2vox_id ;
int *if_visit_r, *num_hkl2r, *hkl2r_mpi, *relevant_pix_id ;
int *len_d2r, *count_d2r, **idx_d2r, **box_vox_id ;
int *my_Acount, *my_Alen, *my_pval_len, *qid_map, *hkl_table ;
int num_row, num_col, num_data, s_num_data, iter_data_block ;
int VN, qlen, data_block_sz, box_half_len, box_len, box_len2 ;
int nproc, myid, num_coarse_rot, num_rot, s_num_rot ;
int num_pix, full_num_pix, iter_flag, Nx, Ny, Nz ;
int hmax, kmax, lmax, hlen, klen, llen, num_hkl ;
int max_num_r2hkl, max_num_hkl2r, max_r_block_sz ;
float (*quat)[5], (*pix)[3], *polar, *inv_polar, *log_polar ;
float *tomo_val, **my_Aval, **my_Bval ;
long long NyNz, **hkl2vox_id ;
short **data_frame ;
E_matrix *Emat ;

int main(int argc, char* argv[]){

    MPI_Init(&argc, &argv) ;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc) ;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid) ;

    FILE *fp ;
    int i, iter ;
    double t1, t2 ;

    t1 = MPI_Wtime() ;
    logP_MIN = log(P_MIN) ;
    if (argc == 3){
        sprintf(config_file, "%s", argv[1]) ;
        iter = atoi(argv[2]) ;
    }
    else{
        if (myid == 0)
            printf("expected two argument: config_file & iter\n") ;
        return 0 ;
    }

    if (!setup())
        return 0 ;

    initialize() ;
    MPI_Barrier(MPI_COMM_WORLD) ;
    t2 = MPI_Wtime() ;
    if (myid == 0)
        printf("setup takes %.0f sec\n", t2-t1) ;

    for (i = 1 ; i <= iter ; i++){
        if (myid == 0 && i == 1){
            fp = fopen("EMC.log", "w") ;
            fprintf(fp, "nproc = %d    num_coarse_rot = %d    ", nproc, num_coarse_rot) ;
            fprintf(fp, "num_rot = %d    num_data = %d    ", num_rot, num_data) ;
            fprintf(fp, "num_pix = %d\n", num_pix) ;
            fclose(fp) ;
        }

        iter_flag = i ;
        if (myid == 0)
            printf("\n --- %d-th iteration: --- \n\n", i) ;

        emc() ;
    }

    free_mem() ;
    MPI_Barrier(MPI_COMM_WORLD) ;
    MPI_Finalize() ;
    return 0 ;
}
