#ifndef emc_hdr
#define emc_hdr

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include <sys/stat.h>

extern const long long int_max ;
extern const double golden_ratio ;
extern const double CP_THRES ;
extern const double P_MIN ;
extern const int mem_step_sz ;

typedef struct{
    int npix ;
    int *pix_id ;
    int *vox_ct ;
    int *vox_id ;
    double *weight ;
} E_matrix ;

extern char config_file[256], data_file[256], peak2r_file[256] ;
extern char start_intens_file[256], home_dir[256], prob_dir[256] ;
extern double **ave_bg, **log_bg, *phi, *updated_phi ;
extern double **intens1, *logP_offset, **pval_d2r ;
extern double rlatt_vec[3][3], proj[3][3], rot[3][3], logP_MIN ;
extern double qval_max, rmin2, rmax2, wl_D, gw, min_rcell_sz ;
extern int *r_block_id, *r_block_sz_mpi, **r_block_id_mpi ;
extern int *rid_map, *pix_map, **prob_orien, *num_prob_orien ;
extern int *num_r2hkl, **r2hkl_mpi, **r2hkl, *num_hkl2vox_id ;
extern int *if_visit_r, *num_hkl2r, *hkl2r_mpi, *relevant_pix_id ;
extern int *len_d2r, *count_d2r, **idx_d2r, **box_vox_id, *tomo_idx ;
extern int **my_Aidx, **my_pix_ct, **my_pval_id, *my_pval_ct ;
extern int *my_Acount, *my_Alen, *my_pval_len, *qid_map, *hkl_table ;
extern int num_row, num_col, num_data, s_num_data, iter_data_block ;
extern int VN, qlen, data_block_sz, box_half_len, box_len, box_len2 ;
extern int nproc, myid, num_rot, s_num_rot, num_pix, full_num_pix ;
extern int Nx, Ny, Nz, hmax, kmax, lmax, hlen, klen, llen, iter_flag ;
extern int num_hkl, max_num_r2hkl, max_num_hkl2r, max_r_block_sz ;
extern float (*quat)[5], (*pix)[3], *polar, *inv_polar, *log_polar ;
extern float *tomo_val, **my_Aval, **my_Bval ;
extern long long NyNz, **hkl2vox_id ;
extern short **data_frame ;
extern E_matrix *Emat ;

/* setup.c */
int setup() ;
void set_intens() ;
void free_mem() ;
void initialize() ;
void print_intens() ;

/* tomo.c */
void make_rot( int ) ;
void calprob( int ) ;
void compress( int ) ;

/* max.c */
void emc() ;
void maximize( int, int ) ;
void sym_intens() ;

/* Aval.c */
void calAval( int ) ;

/* merge-sort.c */
void merge_sort( double*, double*, int*, int*, int ) ;
void TopDownSplitMerge( double*, double*, int*, int*, int, int ) ;
void TopDownMerge( double*, double*, int*, int*, int, int, int ) ;
void merge_sort_idx( double*, double*, int*, int*, int ) ;
void TopDownSplitMergeIdx( double*, double*, int*, int*, int, int ) ;
void TopDownMergeIdx( double*, double*, int*, int*, int, int, int ) ;
void CopyArray( double*, double*, int*, int*, int, int ) ;
#endif
