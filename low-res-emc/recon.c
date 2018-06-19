/*
low-resolution intensity reconstruction using the EMC algorithm

usage:
mpirun -np [number of processors] ./emc (-intens) [path of the config_file] [number of iterations]
"-intens": only update intensity; otherwise, update intensity & phi alternately

needs:
config_file, mpi_datafile, mpi_bgfile, r2peak_file, peak2r_file,
quat_file, rec-vectors.bin, rec2pix-map.dat, pix-map.dat, basis-vec.dat 

makes:
probability distribution, finish_intensity.bin, hkl-val.dat, total-phi.dat

*/

#include "emc.h"

const long long int_max = 2147483647 ;
const double golden_ratio = 1.618034 ;
const double CP_THRES = 0.999 ;
const double P_MIN = 0.0001 ;
const int mem_step_sz = 128 ;

char config_file[256], data_file[256], peak2r_file[256] ;
char start_intens_file[256], home_dir[256], prob_dir[256] ;
double **ave_bg, **log_bg, *phi, *updated_phi ;
double **intens1, *logP_offset, **pval_d2r ;
double rlatt_vec[3][3], proj[3][3], rot[3][3], logP_MIN ;
double qval_max, rmin2, rmax2, wl_D, gw, min_rcell_sz ;
int *r_block_id, *r_block_sz_mpi, **r_block_id_mpi ;
int *rid_map, *pix_map, **prob_orien, *num_prob_orien ;
int *num_r2hkl, **r2hkl_mpi, **r2hkl, *num_hkl2vox_id ;
int *if_visit_r, *num_hkl2r, *hkl2r_mpi, *relevant_pix_id ;
int *len_d2r, *count_d2r, **idx_d2r, **box_vox_id, *tomo_idx ;
int **my_Aidx, **my_pix_ct, **my_pval_id, *my_pval_ct ;
int *my_Acount, *my_Alen, *my_pval_len, *qid_map, *hkl_table ;
int num_row, num_col, num_data, s_num_data, iter_data_block ;
int VN, qlen, data_block_sz, box_half_len, box_len, box_len2 ;
int nproc, myid, num_rot, s_num_rot, num_pix, full_num_pix ;
int Nx, Ny, Nz, hmax, kmax, lmax, hlen, klen, llen, iter_flag ;
int num_hkl, max_num_r2hkl, max_num_hkl2r, max_r_block_sz ;
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
    int i, iter, if_cal_intens = 0 ;
    double t1, t2 ;

    t1 = MPI_Wtime() ;
    logP_MIN = log(P_MIN) ;
    if (argc == 4){
        if (strstr(argv[1], "-intens"))
            if_cal_intens = 1 ;
        sprintf(config_file, "%s", argv[2]) ;
        iter = atoi(argv[3]) ;
        if (myid == 0)
            printf("update intensity only\n") ;
    }
    else{
        sprintf(config_file, "%s", argv[1]) ;
        iter = atoi(argv[2])*2 ;
        if (myid == 0)
            printf("update intensity and phi alternately\n") ;
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
            fprintf(fp, "nproc = %d    num_rot = %d    num_data = %d    num_pix = %d\n", \
                nproc, num_rot, num_data, num_pix) ;
            fclose(fp) ;
        }

        if (if_cal_intens == 1)
            iter_flag = 2*i - 1 ;
        else
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
