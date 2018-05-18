/*
write data in short integers in binary format

usage:
./wr-data [path of the config_file] >> run.log &

needs:
config_file, prob-orien.bin, reduced-cbf_files.dat, reduced-data_id.dat
num_prob_orien.dat, rec-vectors.bin, rec2pix-map.dat

makes:
mpi_datafile 

*/

#include "write-data.h"

const int short_max = 32767 ;

frame_info *data_info ;
filename *cbf_files, *peakfiles ;
char config_file[256], home_dir[256], mpi_datafile[256] ;
int num_row, num_col, num_raw_data, num_data, num_pix, max_patch_sz ;
int *data_id, *rec2pix_map, *pix_map ;
double hot_pix_thres ;

int main(int argc, char *argv[]){
    
    time_t t1, t2 ;
    time(&t1) ;

    if (argc < 2){
        printf("The config_file is missing!!\n") ;
        exit(1) ;
    }
    else
        sprintf(config_file, "%s", argv[1]) ;

    printf("run write-data.c:\n\n") ;
    if (!setup()){
        printf("setup fails!!\n") ;
        exit(1) ;
    }

    if (!check_reduce()){
        printf("reduce_data fails!!\n") ;
        exit(1) ;
    }
    else
        printf("reduce_data succeeds!!\n") ;

    write_data() ;

    time(&t2) ;
    printf ("elapsed time: %.f sec\n", difftime(t2, t1));

    free_mem() ;
    return 0 ;
}


void write_data(){

    FILE *fp, *fn ;
    int d, t, pid, idx, pix_id, pix_val ;
    int total_pix, photon_count, patch_sz, *det ;
    short *data_frame ;

    total_pix = num_row*num_col ;
    det = malloc(total_pix * sizeof(int)) ;
    data_frame = malloc(num_pix * sizeof(short)) ;

    fp = fopen(mpi_datafile, "wb") ;
    for (d = 0 ; d < num_data ; d++){
        read_cbf(det, d) ;
        idx = data_id[d] ;
        for (t = 0 ; t < num_pix ; t++){
            pid = rec2pix_map[t] ;
            if (det[pid] < 0 || det[pid] > hot_pix_thres || det[pid] > short_max)
                photon_count = -1 ;
            else
                photon_count = det[pid] ;
            data_frame[t] = photon_count ;
        }
        
        fn = fopen(peakfiles[d].name, "r") ;
        while (1 == fscanf(fn, "%d", &patch_sz)){
            if (patch_sz > max_patch_sz){
                for (t = 0 ; t < patch_sz ; t++){
                    fscanf(fn, "%d %d", &pix_id, &pix_val) ;
                    pid = pix_map[pix_id] ;
                    if (pid < 0)
                        continue ;
                    else
                        data_frame[pid] = -1 ;
                }
            }
            else{
                for (t = 0 ; t < patch_sz ; t++)
                    fscanf(fn, "%d %d", &pix_id, &pix_val) ;
            }
        }
        fclose(fn) ;

        fwrite(data_frame, sizeof(short), num_pix, fp) ;
    }
    fclose(fp) ;

    free(det) ;
    free(data_frame) ;
}


int check_reduce(){
    
    FILE *fp ;
    char infile[256] ;
    int d, idx, flag, *count_d2r, *num_orien ;

    count_d2r = malloc(num_data * sizeof(int)) ;
    num_orien = malloc(num_raw_data * sizeof(int)) ;

    sprintf(infile, "%s/aux/prob-orien.bin", home_dir) ;
    fp = fopen(infile, "rb") ;
    for (d = 0 ; d < num_data ; d++){
        fread(&count_d2r[d], sizeof(int), 1, fp) ;
        fseek(fp, count_d2r[d]*sizeof(int), SEEK_CUR) ;
    }
    fclose(fp) ;

    sprintf(infile, "%s/orient-peak/num_prob_orien.dat", home_dir) ;
    fp = fopen(infile, "r") ;
    for (d = 0 ; d < num_raw_data ; d++)
        fscanf(fp, "%d", &num_orien[d]) ;
    fclose(fp) ;

    flag = 1 ;
    for (d = 0 ; d < num_data ; d++){
        idx = data_id[d] ;
        if (count_d2r[d] != num_orien[idx])
            flag = 0 ;
    }

    free(count_d2r) ;
    free(num_orien) ;
    
    return flag ;
}


int setup(){
    
    FILE *fp ;
    int d, t, qmax, qmin, total_pix ;
    char *token, line[256], cbf_filelist[256], peakfilelist[256], infile[256] ;

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
            if (strcmp(token, "num_row") == 0)
                num_row = atoi(strtok(NULL, " =")) ;
            else if (strcmp(token, "num_col") == 0)
                num_col = atoi(strtok(NULL, " =")) ;
            else if (strcmp(token, "num_raw_data") == 0)
                num_raw_data = atoi(strtok(NULL, " =")) ;
            else if (strcmp(token, "max_patch_sz") == 0)
                max_patch_sz = atoi(strtok(NULL, " =")) ;
            else if (strcmp(token, "hot_pix_thres") == 0)
                hot_pix_thres = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "home_dir") == 0)
                strcpy(home_dir, strtok(NULL, " =\n")) ;
            else if (strcmp(token, "mpi_datafile") == 0)
                strcpy(mpi_datafile, strtok(NULL, " =\n")) ;
        }
        fclose(fp) ;
    }

    fp = fopen("reduced-data_id.dat", "r") ;
    fscanf(fp, "%d", &num_data) ;
    data_id = malloc(num_data * sizeof(int)) ;
    for (d = 0 ; d < num_data ; d++)
        fscanf(fp, "%d", &data_id[d]) ;
    fclose(fp) ;

    data_info = malloc(num_data * sizeof(frame_info)) ;

    sprintf(cbf_filelist, "%s/aux/reduced-cbf_files.dat", home_dir) ;
    cbf_files = malloc(num_data * sizeof(filename)) ;
    fp = fopen(cbf_filelist, "r") ;
    for (d = 0 ; d < num_data ; d++)
        fscanf(fp, "%s", cbf_files[d].name) ;
    fclose(fp) ;

    sprintf(peakfilelist, "%s/aux/reduced-peakfiles.dat", home_dir) ;
    peakfiles = malloc(num_data * sizeof(filename)) ;
    fp = fopen(peakfilelist, "r") ;
    for (d = 0 ; d < num_data ; d++)
        fscanf(fp, "%s", peakfiles[d].name) ;
    fclose(fp) ;

    sprintf(infile, "%s/make-detector/rec-vectors.bin", home_dir) ;
    fp = fopen(infile, "rb") ;
    fread(&qmax, sizeof(int), 1, fp) ;
    fread(&qmin, sizeof(int), 1, fp) ;
    fread(&num_pix, sizeof(int), 1, fp) ;
    fclose(fp) ;

    rec2pix_map = malloc(num_pix * sizeof(int)) ;
    sprintf(infile, "%s/make-detector/rec2pix-map.dat", home_dir) ;
    fp = fopen(infile, "r") ;
    for (t = 0 ; t < num_pix ; t++)
        fscanf(fp, "%d", &rec2pix_map[t]) ;
    fclose(fp) ;

    total_pix = num_row*num_col ;
    pix_map = malloc(total_pix * sizeof(int)) ;
    sprintf(infile, "%s/make-detector/pix-map.dat", home_dir) ;
    fp = fopen(infile, "r") ;
    for (t = 0 ; t < total_pix ; t++)
        fscanf(fp, "%d", &pix_map[t]) ;
    fclose(fp) ;

    return 1 ;
}


void free_mem(){

    free(data_id) ;
    free(data_info) ;
    free(cbf_files) ;
    free(peakfiles) ;
    free(rec2pix_map) ;
    free(pix_map) ;
}
