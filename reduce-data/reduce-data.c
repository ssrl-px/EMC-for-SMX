/*
reduce data and balance the work loads between processors 

compile:
gcc reduce-data.c -O3 -lm -o reduce-data

usage:
./reduce-data [path of the config_file] >> run.log &

needs:
config_file, quat_file

makes:
prob-orien.bin, mpi_bgfile, reduced-cbf_files.dat
reduced-radialfiles.dat, reduced-peakfiles.dat, reduced-data_id.dat

*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

typedef struct{
    char name[256] ;
} filename ;

filename *orienfiles, *cbf_files ;
filename *radialfiles, *peakfiles ;
char config_file[256], home_dir[256] ;
char quat_file[256], mpi_bgfile[256] ;
int num_raw_data, num_data, nproc, num_div, num_rot, qlen ;
int *num_orien, **prob_orien, *load_balanced_data_id ;

int setup() ;
void load_balance() ;
void write_files() ;
void free_mem() ;

/* merge sort */
void merge_sort( int*, int*, int*, int*, int ) ;
void TopDownSplitMerge( int*, int*, int*, int*, int, int ) ;
void TopDownMerge( int*, int*, int*, int*, int, int, int ) ;
void CopyArray( int*, int*, int*, int*, int, int ) ;

int main(int argc, char *argv[]){
    
    time_t t1, t2 ;
    time(&t1) ;

    if (argc < 2){
        printf("The config_file is missing!!\n") ;
        exit(1) ;
    }
    else
        sprintf(config_file, "%s", argv[1]) ;

    printf("run reduce-data.c:\n\n") ;
    if (!setup()){
        printf("setup fails!!\n") ;
        exit(1) ;
    }

    load_balance() ;
    write_files() ;

    time(&t2) ;
    printf ("elapsed time: %.f sec\n\n\n", difftime(t2, t1));

    free_mem() ;
    return 0 ;
}


void load_balance(){
    
    int i, j, d, r, rid, idx, idx_offset ;
    int s_num_data, count_d2r_d, *prob_orien_d ;
    int *count_d2r, *cpy_count_d2r, *data_id, *cpy_data_id ;
    int *frame_count, *proc_num_orien, *cpy_proc_num_orien ;
    int **proc_occupancy, **proc_data_id, *proc_id, *cpy_proc_id ;

    count_d2r = malloc(num_data * sizeof(int)) ;
    cpy_count_d2r = malloc(num_data * sizeof(int)) ;
    data_id = malloc(num_data * sizeof(int)) ;
    cpy_data_id = malloc(num_data * sizeof(int)) ;

    num_data = 0 ;
    for (d = 0 ; d < num_raw_data ; d++){
        if (num_orien[d] > 0){
            count_d2r[num_data] = num_orien[d] ;
            data_id[num_data] = d ;
            num_data += 1 ;
        }
    }

    merge_sort(count_d2r, cpy_count_d2r, data_id, cpy_data_id, num_data) ;
    if (num_data % nproc == 0)
        s_num_data = num_data / nproc ;
    else
        s_num_data = num_data / nproc + 1 ;

    frame_count = calloc(nproc, sizeof(int)) ;
    proc_num_orien = calloc(nproc, sizeof(int)) ;
    cpy_proc_num_orien = malloc(nproc * sizeof(int)) ;
    proc_occupancy = malloc(nproc * sizeof(int *)) ;
    proc_data_id = malloc(nproc * sizeof(int *)) ;
    for (i = 0 ; i < nproc ; i++){
        proc_occupancy[i] = calloc(num_rot, sizeof(int)) ;
        proc_data_id[i] = malloc(s_num_data * sizeof(int)) ;
    }

    proc_id = malloc(nproc * sizeof(int)) ;
    cpy_proc_id = malloc(nproc * sizeof(int)) ;
    for (i = 0 ; i < nproc ; i++)
        proc_id[i] = i ;
    
    for (d = 0 ; d < s_num_data-1 ; d++){
        idx_offset = d*nproc ;
        merge_sort(proc_num_orien, cpy_proc_num_orien, proc_id, cpy_proc_id, nproc) ;
        for (i = 0 ; i < nproc ; i++){
            idx = idx_offset + i ;
            count_d2r_d = count_d2r[idx] ;
            prob_orien_d = prob_orien[data_id[idx]] ;
            j = nproc - 1 - i ;
            for (r = 0 ; r < count_d2r_d ; r++){
                rid = prob_orien_d[r] ;
                if (proc_occupancy[proc_id[j]][rid] == 0){
                    proc_occupancy[proc_id[j]][rid] = 1 ;
                    proc_num_orien[j] += 1 ;
                }
            }
            proc_data_id[proc_id[j]][d] = data_id[idx] ;
        }
    }
    
    merge_sort(proc_num_orien, cpy_proc_num_orien, proc_id, cpy_proc_id, nproc) ;
    for (i = 0 ; i < nproc ; i++)
        frame_count[i] = s_num_data-1 ;

    idx_offset = (s_num_data-1)*nproc ;
    for (d = idx_offset ; d < num_data ; d++){
        j = nproc - 1 - d % nproc ;
        count_d2r_d = count_d2r[d] ;
        prob_orien_d = prob_orien[data_id[d]] ;
        for (r = 0 ; r < count_d2r_d ; r++){
            rid = prob_orien_d[r] ;
            if (proc_occupancy[proc_id[j]][rid] == 0){
                proc_occupancy[proc_id[j]][rid] = 1 ;
                proc_num_orien[j] += 1 ;
            }
        }
        proc_data_id[proc_id[j]][s_num_data-1] = data_id[d] ;
        frame_count[j] += 1 ;
    }

    load_balanced_data_id = malloc(num_data * sizeof(int)) ;
    num_data = 0 ;
    for (i = 0 ; i < nproc ; i++){
        j = nproc - 1 - i ;
        for (d = 0 ; d < frame_count[j] ; d++){
            idx = proc_data_id[proc_id[j]][d] ;
            load_balanced_data_id[num_data] = idx ;
            num_data += 1 ;
        }
        printf("myid = %d: frame_count = %d\n", i, frame_count[j]) ;
    }
    
    free(count_d2r) ;
    free(cpy_count_d2r) ;
    free(data_id) ;
    free(cpy_data_id) ;
    free(frame_count) ;
    free(proc_num_orien) ;
    free(cpy_proc_num_orien) ;
    for (i = 0 ; i < nproc ; i++){
        free(proc_occupancy[i]) ;
        free(proc_data_id[i]) ;
    }
    free(proc_occupancy) ;
    free(proc_data_id) ;
    free(proc_id) ;
    free(cpy_proc_id) ;
}


void write_files(){
    
    FILE *fp, *fn ;
    char outfile[256] ;
    int d, idx ;
    double *ave_bg ;

    sprintf(outfile, "%s/aux/prob-orien.bin", home_dir) ;
    fp = fopen(outfile, "wb") ;
    for (d = 0 ; d < num_data ; d++){
        idx = load_balanced_data_id[d] ;
        fwrite(&num_orien[idx], sizeof(int), 1, fp) ;
        fwrite(prob_orien[idx], sizeof(int), num_orien[idx], fp) ;
    }
    fclose(fp) ;
    
    sprintf(outfile, "%s/aux/reduced-cbf_files.dat", home_dir) ;
    fp = fopen(outfile, "w") ;
    for (d = 0 ; d < num_data ; d++){
        idx = load_balanced_data_id[d] ;
        fprintf(fp, "%s\n", cbf_files[idx].name) ;
    }
    fclose(fp) ;

    sprintf(outfile, "%s/aux/reduced-radialfiles.dat", home_dir) ;
    fp = fopen(outfile, "w") ;
    for (d = 0 ; d < num_data ; d++){
        idx = load_balanced_data_id[d] ;
        fprintf(fp, "%s\n", radialfiles[idx].name) ;
    }
    fclose(fp) ;
    
    sprintf(outfile, "%s/aux/reduced-peakfiles.dat", home_dir) ;
    fp = fopen(outfile, "w") ;
    for (d = 0 ; d < num_data ; d++){
        idx = load_balanced_data_id[d] ;
        fprintf(fp, "%s\n", peakfiles[idx].name) ;
    }
    fclose(fp) ;

    ave_bg = malloc(qlen * sizeof(double)) ;
    fp = fopen(mpi_bgfile, "wb") ;
    for (d = 0 ; d < num_data ; d++){
        idx = load_balanced_data_id[d] ;
        fn = fopen(radialfiles[idx].name, "rb") ;
        fread(ave_bg, sizeof(double), qlen, fn) ;
        fclose(fn) ;
        fwrite(ave_bg, sizeof(double), qlen, fp) ;
    }
    fclose(fp) ;

    fp = fopen("reduced-data_id.dat", "w") ;
    fprintf(fp, "%d\n", num_data) ;
    for (d = 0 ; d < num_data ; d++)
        fprintf(fp, "%d\n", load_balanced_data_id[d]) ;
    fclose(fp) ;

    free(ave_bg) ;
}


int setup(){
    
    FILE *fp ;
    int i, d, num ;
    char *token, line[256], infile[256] ;
    char orienfilelist[256], cbf_filelist[256] ;
    char radialfilelist[256], peakfilelist[256] ;

    fp = fopen(config_file, "r") ;
    if (!fp){
        printf("The config_file %s is not found!!\n", config_file) ;
        return 0 ;
    }
    else{
        while (fgets(line, 256, fp) != NULL){
            token = strtok(line, " =") ;
            if (token[0] == '#' || token[0] == '\n' || token[0] == '[')
                continue ;
            if (strcmp(token, "num_raw_data") == 0)
                num_raw_data = atoi(strtok(NULL, " =")) ;
            else if (strcmp(token, "nproc") == 0)
                nproc = atoi(strtok(NULL, " =")) ;
            else if (strcmp(token, "qlen") == 0)
                qlen = atoi(strtok(NULL, " =")) ;
            else if (strcmp(token, "home_dir") == 0)
                strcpy(home_dir, strtok(NULL, " =\n")) ;
            else if (strcmp(token, "quat_file") == 0)
                strcpy(quat_file, strtok(NULL, " =\n")) ;
            else if (strcmp(token, "mpi_bgfile") == 0)
                strcpy(mpi_bgfile, strtok(NULL, " =\n")) ;
        }
        fclose(fp) ;
    }

    fp = fopen(quat_file, "rb") ;
    fread(&num_rot, sizeof(int), 1, fp) ;
    fclose(fp) ;
    
    token = strtok(quat_file, "/-.") ;
    while(token != NULL){
        token = strtok(NULL, "/-.") ;
        if (token == NULL)
            break ;
        if (strstr(token, "quaternion") != NULL)
            num_div = atoi(token + strlen("quaternion")) ;
    }

    orienfiles = malloc(num_raw_data * sizeof(filename)) ;
    sprintf(orienfilelist, "%s/make-background/orienlist.dat", home_dir) ;
    fp = fopen(orienfilelist, "r") ;
    for (d = 0 ; d < num_raw_data ; d++)
        fscanf(fp, "%s", orienfiles[d].name) ;
    fclose(fp) ;

    cbf_files = malloc(num_raw_data * sizeof(filename)) ;
    sprintf(cbf_filelist, "%s/make-background/cbflist.dat", home_dir) ;
    fp = fopen(cbf_filelist, "r") ;
    for (d = 0 ; d < num_raw_data ; d++)
        fscanf(fp, "%s", cbf_files[d].name) ;
    fclose(fp) ;

    radialfiles = malloc(num_raw_data * sizeof(filename)) ;
    sprintf(radialfilelist, "%s/make-background/radiallist.dat", home_dir) ;
    fp = fopen(radialfilelist, "r") ;
    for (d = 0 ; d < num_raw_data ; d++)
        fscanf(fp, "%s", radialfiles[d].name) ;
    fclose(fp) ;
    
    peakfiles = malloc(num_raw_data * sizeof(filename)) ;
    sprintf(peakfilelist, "%s/make-background/peaklist.dat", home_dir) ;
    fp = fopen(peakfilelist, "r") ;
    for (d = 0 ; d < num_raw_data ; d++)
        fscanf(fp, "%s", peakfiles[d].name) ;
    fclose(fp) ;

    num_orien = calloc(num_raw_data, sizeof(int)) ;
    prob_orien = malloc(num_raw_data * sizeof(int *)) ;
    sprintf(infile, "%s/orient-peak/num_prob_orien.dat", home_dir) ;
    fp = fopen(infile, "r") ;
    for (d = 0 ; d < num_raw_data ; d++)
        fscanf(fp, "%d", &num_orien[d]) ;
    fclose(fp) ;

    for (d = 0 ; d < num_raw_data ; d++){
        if (num_orien[d] == 0)
            continue ;
        fp = fopen(orienfiles[d].name, "r") ;
        fscanf(fp, "%d", &num) ;
        if (num != num_orien[d])
            printf("error in num_orien[d]!!\n") ;
        prob_orien[d] = malloc(num_orien[d] * sizeof(int)) ;
        for (i = 0 ; i < num_orien[d] ; i++)
            fscanf(fp, "%d", &prob_orien[d][i]) ;
        fclose(fp) ;
    }

    num_data = 0 ;
    for (d = 0 ; d < num_raw_data ; d++){
        if (num_orien[d] > 0)
            num_data += 1 ;
    }
    
    printf("num_raw_data = %d, num_data = %d\n", num_raw_data, num_data) ;
    printf("num_rot = %d, num_div = %d\n", num_rot, num_div) ;
    
    return 1 ;
}


void free_mem(){

    int d ;

    free(orienfiles) ;
    free(cbf_files) ;
    free(radialfiles) ;
    free(peakfiles) ;
    for (d = 0 ; d < num_raw_data ; d++){
        if (num_orien[d] > 0)
            free(prob_orien[d]) ;
    }
    free(prob_orien) ;
    free(num_orien) ;
    free(load_balanced_data_id) ;
}


/* sort values in descending order */
void merge_sort( int *A, int *B, int *A_idx, int *B_idx, int length ){
    TopDownSplitMerge(A, B, A_idx, B_idx, 0, length) ;
}


void TopDownSplitMerge( int *A, int *B, int *A_idx, int *B_idx, int iBegin, int iEnd ){

    /* iBegin: inclusive, iEnd: exclusive */
    if (iEnd - iBegin == 1)
        return ;

    int iMiddle = (iEnd + iBegin) / 2 ;
    TopDownSplitMerge(A, B, A_idx, B_idx, iBegin, iMiddle) ;
    TopDownSplitMerge(A, B, A_idx, B_idx, iMiddle, iEnd) ;
    TopDownMerge(A, B, A_idx, B_idx, iBegin, iMiddle, iEnd) ;
    CopyArray(A, B, A_idx, B_idx, iBegin, iEnd) ;
}


void TopDownMerge( int *A, int *B, int *A_idx, int *B_idx, int iBegin, int iMiddle, int iEnd ){

    int idx, i0 = iBegin, j0 = iMiddle ;

    for (idx = iBegin ; idx < iEnd ; idx++){

        if ( (i0 < iMiddle) && ((j0 >= iEnd) || (A[i0] > A[j0])) ){
            B[idx] = A[i0] ;
            B_idx[idx] = A_idx[i0] ;
            i0 += 1 ;
        }
        else{
            B[idx] = A[j0] ;
            B_idx[idx] = A_idx[j0] ;
            j0 += 1 ;
        }
    }
}


void CopyArray( int *A, int *B, int *A_idx, int *B_idx, int iBegin, int iEnd ){

    int i ;
    for (i = iBegin ; i < iEnd ; i++){
        A[i] = B[i] ;
        A_idx[i] = B_idx[i] ;
    }
}
