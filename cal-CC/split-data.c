/*
split data into two subsets: odd and even

compile:
gcc split-data.c -O3 -lm -o split-data.o

usage:
./split-data.o [path of the config_file] [path to the data dir]

needs:
config_file, reduced-data_id.dat, rec-vectors.bin
mpi_datafile, mpi_bgfile, prob-orien.bin

makes:
mpi_datafile-A.bin, mpi-bg_model-A.bin, prob-orien-A.bin
mpi_datafile-B.bin, mpi-bg_model-B.bin, prob-orien-B.bin

*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

const int mem_step_sz = 128 ;
char config_file[256], data_dir[256], home_dir[256] ;
char mpi_datafile[256], mpi_bgfile[256] ;
int num_data, num_pix, qlen ;

int setup() ;
void split_data() ;

int main(int argc, char *argv[]){
    
    time_t t1, t2 ;
    time(&t1) ;

    if (argc < 3){
        printf("need two arguments: config_file & data_dir!!\n") ;
        exit(1) ;
    }
    else{
        sprintf(config_file, "%s", argv[1]) ;
        sprintf(data_dir, "%s", argv[2]) ;
    }

    if (!setup()){
        printf("setup fails!!\n") ;
        exit(1) ;
    }

    split_data() ;

    time(&t2) ;
    printf ("elapsed time: %.f sec\n", difftime(t2, t1));
    return 0 ;
}


void split_data(){

    FILE *fp, *fa, *fb ;
    char infile[256], outfile_a[256], outfile_b[256] ;
    int d, num, len, num_data_a, num_data_b, *prob_orien ;
    short *data_frame ;
    double *ave_bg ;

    num_data_b = num_data / 2 ;
    num_data_a = num_data - num_data_b ;
    data_frame = malloc(num_pix * sizeof(short)) ;
    ave_bg = malloc(qlen * sizeof(double)) ;

    sprintf(outfile_a, "%s/mpi-datafile-A.bin", data_dir) ;
    sprintf(outfile_b, "%s/mpi-datafile-B.bin", data_dir) ;

    fp = fopen(mpi_datafile, "rb") ;
    fa = fopen(outfile_a, "wb") ;
    fb = fopen(outfile_b, "wb") ;
    for (d = 0 ; d < num_data ; d++){
        fread(data_frame, sizeof(short), num_pix, fp) ;
        if (d % 2 == 0)
            fwrite(data_frame, sizeof(short), num_pix, fa) ;
        else
            fwrite(data_frame, sizeof(short), num_pix, fb) ;
    }
    fclose(fp) ;
    fclose(fa) ;
    fclose(fb) ;

    sprintf(outfile_a, "%s/mpi-bg_model-A.bin", data_dir) ;
    sprintf(outfile_b, "%s/mpi-bg_model-B.bin", data_dir) ;
    fp = fopen(mpi_bgfile, "rb") ;
    fa = fopen(outfile_a, "wb") ;
    fb = fopen(outfile_b, "wb") ;
    for (d = 0 ; d < num_data ; d++){
        fread(ave_bg, sizeof(double), qlen, fp) ;
        if (d % 2 == 0)
            fwrite(ave_bg, sizeof(double), qlen, fa) ;
        else
            fwrite(ave_bg, sizeof(double), qlen, fb) ;
    }
    fclose(fp) ;
    fclose(fa) ;
    fclose(fb) ;

    sprintf(infile, "%s/aux/prob-orien.bin", home_dir) ;
    sprintf(outfile_a, "%s/aux/prob-orien-A.bin", home_dir) ;
    sprintf(outfile_b, "%s/aux/prob-orien-B.bin", home_dir) ;

    len = mem_step_sz ;
    prob_orien = malloc(len * sizeof(int)) ;

    fp = fopen(infile, "rb") ;
    fa = fopen(outfile_a, "wb") ;
    fb = fopen(outfile_b, "wb") ;
    for (d = 0 ; d < num_data ; d++){
        fread(&num, sizeof(int), 1, fp) ;
        if (num > len){
            len = num ;
            prob_orien = realloc(prob_orien, len*sizeof(int)) ;
        }
        fread(&prob_orien[0], sizeof(int), num, fp) ;
        if (d % 2 == 0){
            fwrite(&num, sizeof(int), 1, fa) ;
            fwrite(&prob_orien[0], sizeof(int), num, fa) ;
        }
        else{
            fwrite(&num, sizeof(int), 1, fb) ;
            fwrite(&prob_orien[0], sizeof(int), num, fb) ;
        }
    }
    fclose(fp) ;
    fclose(fa) ;
    fclose(fb) ;

    printf("num_data_a = %d, num_data_b = %d\n", num_data_a, num_data_b) ;
    free(data_frame) ;
    free(ave_bg) ;
    free(prob_orien) ;
}


int setup(){
    
    FILE *fp ;
    int qmax, qmin ;
    char *token, line[256], infile[256] ;

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
            if (strcmp(token, "qlen") == 0)
                qlen = atoi(strtok(NULL, " =")) ;
            else if (strcmp(token, "home_dir") == 0)
                strcpy(home_dir, strtok(NULL, " =\n")) ;
            else if (strcmp(token, "mpi_datafile") == 0)
                strcpy(mpi_datafile, strtok(NULL, " =\n")) ;
            else if (strcmp(token, "mpi_bgfile") == 0)
                strcpy(mpi_bgfile, strtok(NULL, " =\n")) ;
        }
        fclose(fp) ;
    }

    sprintf(infile, "%s/reduce-data/reduced-data_id.dat", home_dir) ;
    fp = fopen(infile, "r") ;
    fscanf(fp, "%d", &num_data) ;
    fclose(fp) ;

    sprintf(infile, "%s/make-detector/rec-vectors.bin", home_dir) ;
    fp = fopen(infile, "rb") ;
    fread(&qmax, sizeof(int), 1, fp) ;
    fread(&qmin, sizeof(int), 1, fp) ;
    fread(&num_pix, sizeof(int), 1, fp) ;
    fclose(fp) ;

    printf("num_data = %d, num_pix = %d\n", num_data, num_pix) ;
    return 1 ;
}
