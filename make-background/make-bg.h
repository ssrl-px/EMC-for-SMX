#ifndef bg_hdr
#define bg_hdr

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
//#include <openssl/md5.h>
//#include <openssl/bio.h>
//#include <openssl/evp.h>

extern const double cdf_thres ;
extern const double vmax ;
extern const double dv ;
extern const int num_iter ;
extern const int max_num_peak_on_ring ;
/* threshold for ring detection */
extern const double bg_thres ;

typedef struct{
    char name[256] ;
} filename ;

typedef struct{
    int row_size ;
    int col_size ;
    double exposure ;
    double wl ;
    double detd ;
    double c_row ;
    double c_col ;
} frame_info ;

extern char *content_md5, *conversions ;
extern char config_file[256], home_dir[256] ;
extern char pixmap_file[256], rvec_file[256] ;
extern char rec2pixmap_file[256] ;

extern int num_row, num_col, num_data, qlen, total_pix ;
extern int nproc, myid, qmax, qmin, num_pix, len_val, outlier_ct ;
extern int *pix_map, *rec2pix_map, *qid_map, *queue ;
extern int *peak_id, *peak_val, *peak_label, *peak_map, *peak_list ;
extern int *det, *radial_ct, *count_thres, (*peak_location)[2] ;
extern double detd, wl, px, hot_pix_thres, dq, cx, cy ;
extern double *radial_val, *radial_weight, *ave_bg, (*pix)[4] ;
extern filename *cbf_files, *outlierfiles, *radialfiles, *peakfiles ;
extern frame_info *data_info ;

/* mpi-make-background.c */
int setup() ;
void free_mem() ;
void make_bg() ;
void peak_finder( int ) ;
void mask_rings() ;
void make_count_thres() ;
int ad_hoc_thres( double ) ;

/* read-cbf.c */
//int read_cbf( int *, int ) ;
int parse_header( char *, char *, char *, char *, int ) ;
void byte_offset( int *, char *, int ) ;
void md5sum( char *, char *, int ) ;

#endif
