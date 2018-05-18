#ifndef wr_hdr
#define wr_hdr

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <openssl/md5.h>
#include <openssl/bio.h>
#include <openssl/evp.h>

extern const int short_max ;

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

extern frame_info *data_info ;
extern filename *cbf_files, *peakfiles ;
extern char config_file[256], home_dir[256], mpi_datafile[256] ;
extern int num_row, num_col, num_raw_data, num_data, num_pix, max_patch_sz ;
extern int *data_id, (*beam_translation)[2], *rec2pix_map, *pix_map ;
extern double hot_pix_thres ;

/* write-data.c */
int setup() ;
int check_reduce() ;
void write_data() ;
void free_mem() ;

/* read-cbf.c */
int read_cbf( int *, int ) ;
int parse_header( char *, char *, char *, char *, int ) ;
void byte_offset( int *, char *, int ) ;
void md5sum( char *, char *, int ) ;

#endif
