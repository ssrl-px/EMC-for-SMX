/*
generate quaternuons used in the local-update scheme

compile:
gcc reduce-quaternions.c -O3 -lm -o reduce.o

usage:
./reduce.o [path of the config_file] fine_quat_file

needs:
config_file, quat_file, fine_quat_file, reduced-data_id.dat
quat_embedding-%d-%d.dat % (num_coarse_div, num_div)

makes:
c-reduced-quaternion%d.bin % num_coarse_div
c-reduced-quaternion%d.bin % num_div
reduced-quaternion-%d-%d.dat % (num_coarse_div, num_div)
reduced probability files

*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

const double P_MIN = 0.01 ;
char config_file[256], home_dir[256], prob_dir[256] ;
char quat_file[256], fine_quat_file[256], quatmapfile[256] ;
int nproc, num_coarse_div, num_div, num_coarse_rot, num_rot, num_data ;
int *num_quat_neighbor, **quat_map, *num_orien, **prob_orien ;
int *if_coarse_visit, *if_visit ;
double (*coarse_quat)[5], (*quat)[5], **high_prob ;

void setup() ;
void free_mem() ;

int main(int argc, char *argv[]){

    FILE *fp ;
    char outfile[256] ;
    int i, d, r, rid, dmin, dmax, s_num_data, rank, rank_max ;
    int new_num_coarse_rot, new_num_rot ;
    int *new_coarse_idx, *new_idx ;

    sprintf(config_file, "%s", argv[1]) ;
    sprintf(fine_quat_file, "%s", argv[2]) ;
    setup() ;

    new_num_coarse_rot = 0 ;
    new_coarse_idx = malloc(num_coarse_rot * sizeof(int)) ;
    for (r = 0 ; r < num_coarse_rot ; r++){
        if (if_coarse_visit[r] == 0)
            new_coarse_idx[r] = -1 ;
        else{
            new_coarse_idx[r] = new_num_coarse_rot ;
            new_num_coarse_rot += 1 ;
        }
    }

    new_num_rot = 0 ;
    new_idx = malloc(num_rot * sizeof(int)) ;
    for (r = 0 ; r < num_rot ; r++){
        if (if_visit[r] == 0)
            new_idx[r] = -1 ;
        else{
            new_idx[r] = new_num_rot ;
            new_num_rot += 1 ;
        }
    }

    printf("new_num_coarse_rot = %d\n", new_num_coarse_rot) ;
    printf("new_num_rot = %d\n", new_num_rot) ;

    sprintf(outfile, "c-reduced-quaternion%d.bin", num_coarse_div) ;
    fp = fopen(outfile, "wb") ;
    fwrite(&new_num_coarse_rot, sizeof(int), 1, fp) ;
    for (r = 0 ; r < num_coarse_rot ; r++){
        if (if_coarse_visit[r] == 0)
            continue ;
        fwrite(coarse_quat[r], sizeof(double), 5, fp) ;
    }
    fclose(fp) ;

    sprintf(outfile, "c-reduced-quaternion%d.bin", num_div) ;
    fp = fopen(outfile, "wb") ;
    fwrite(&new_num_rot, sizeof(int), 1, fp) ;
    for (r = 0 ; r < num_rot ; r++){
        if (if_visit[r] == 0)
            continue ;
        fwrite(quat[r], sizeof(double), 5, fp) ;
    }
    fclose(fp) ;

    sprintf(outfile, "reduced-quaternion-%d-%d.dat", num_coarse_div, num_div) ;
    fp = fopen(outfile, "w") ;
    fprintf(fp, "%d\n\n", new_num_coarse_rot) ;
    for (r = 0 ; r < num_coarse_rot ; r++){
        if (if_coarse_visit[r] == 0)
            continue ;
        fprintf(fp, "%d\n", num_quat_neighbor[r]) ;
        for (i = 0 ; i < num_quat_neighbor[r] ; i++){
            rid = new_idx[quat_map[r][i]] ;
            fprintf(fp, "%d ", rid) ;
            if (rid < 0)
                printf("warning: something wrong with new_idx!!\n") ;
        }
        fprintf(fp, "\n") ;
    }
    fclose(fp) ;
    
    if (num_data % nproc == 0)
        s_num_data = num_data / nproc ;
    else
        s_num_data = num_data / nproc + 1 ;
    rank_max = num_data - nproc*(s_num_data - 1) ;

    for (rank = 0 ; rank < nproc ; rank++){
        if (rank < rank_max){
            dmin = rank*s_num_data ;
            dmax = dmin + s_num_data ;
        }
        else{
            dmin = rank_max*s_num_data + (rank - rank_max)*(s_num_data - 1) ;
            dmax = dmin + (s_num_data - 1) ;
        }
        
        sprintf(outfile, "%s/reduced-high_p-%d.dat", prob_dir, rank) ;
        fp = fopen(outfile, "w") ;
        fprintf(fp, "%d\n\n", dmax-dmin) ;
        for (d = dmin ; d < dmax ; d++){
            fprintf(fp, "%d\n", num_orien[d]) ;
            for (i = 0 ; i < num_orien[d] ; i++){
                rid = new_coarse_idx[prob_orien[d][i]] ;
                fprintf(fp, "%d %.6f ", rid, high_prob[d][i]) ;
                if (rid < 0)
                    printf("warning: something wrong with new_coarse_idx!!\n") ;
            }
            fprintf(fp, "\n") ;
        }
    }

    free(new_coarse_idx) ;
    free(new_idx) ;
    free_mem() ;
    return 0 ;
}


void setup(){
    
    FILE *fp ;
    char *token, line[256], infile[256], reduced_data_id_file[256] ;
    int i, j, d, r, rank, idx, rid, ct, num_frame ;
    double pval ;

    fp = fopen(config_file, "r") ;
    if (!fp){
        printf("The config_file %s is not found!!\n", config_file) ;
        return ;
    }
    else{
        while (fgets(line, 256, fp) != NULL){
            token = strtok(line, " =") ;
            if (token[0] == '#' || token[0] == '\n' || token[0] == '[')
                continue ;
            if (strcmp(token, "nproc") == 0)
                nproc = atoi(strtok(NULL, " =")) ;
            else if (strcmp(token, "home_dir") == 0)
                strcpy(home_dir, strtok(NULL, " =\n")) ;
            else if (strcmp(token, "quat_file") == 0)
                strcpy(quat_file, strtok(NULL, " =\n")) ;
            else if (strcmp(token, "start_prob_dir") == 0)
                strcpy(prob_dir, strtok(NULL, " =\n")) ;
            else if (strcmp(token, "reduced_data_id_file") == 0)
                strcpy(reduced_data_id_file, strtok(NULL, " =\n")) ;
        }
    }
    fclose(fp) ;

    sprintf(line, "%s", quat_file) ;
    token = strtok(line, "/") ;
    while (token != NULL){
        if (strstr(token, "quaternion") != NULL)
            sscanf(token, "%*[^0123456789]%d", &num_coarse_div) ;
        token = strtok(NULL, "/") ;
        if (token == NULL)
            break ;
    }
    
    sprintf(line, "%s", fine_quat_file) ;
    token = strtok(line, "/") ;
    while (token != NULL){
        if (strstr(token, "quaternion") != NULL)
            sscanf(token, "%*[^0123456789]%d", &num_div) ;
        token = strtok(NULL, "/") ;
        if (token == NULL)
            break ;
    }

    printf("num_coarse_div = %d, num_div = %d\n", num_coarse_div, num_div) ;
    
    fp = fopen(quat_file, "rb") ;
    fread(&num_coarse_rot, sizeof(int), 1, fp) ;
    coarse_quat = malloc(num_coarse_rot * sizeof(*coarse_quat)) ;
    for (r = 0 ; r < num_coarse_rot ; r++)
        fread(coarse_quat[r], sizeof(double), 5, fp) ;
    fclose(fp) ;
    
    fp = fopen(fine_quat_file, "rb") ;
    fread(&num_rot, sizeof(int), 1, fp) ;
    quat = malloc(num_rot * sizeof(*quat)) ;
    for (r = 0 ; r < num_rot ; r++)
        fread(quat[r], sizeof(double), 5, fp) ;
    fclose(fp) ;
    
    sprintf(quatmapfile, "quaternion-%d-%d.dat", num_coarse_div, num_div) ;
    fp = fopen(quatmapfile, "r") ;
    fscanf(fp, "%d", &num_coarse_rot) ;
    num_quat_neighbor = malloc(num_coarse_rot * sizeof(int)) ;
    quat_map = malloc(num_coarse_rot * sizeof(int *)) ;
    for (r = 0 ; r < num_coarse_rot ; r++){
        fscanf(fp, "%d", &num_quat_neighbor[r]) ;
        quat_map[r] = malloc(num_quat_neighbor[r] * sizeof(int)) ;
        for (i = 0 ; i < num_quat_neighbor[r] ; i++)
            fscanf(fp, "%d", &quat_map[r][i]) ;
    }
    fclose(fp) ;

    fp = fopen(reduced_data_id_file, "r") ;
    fscanf(fp, "%d", &num_data) ;
    fclose(fp) ;

    prob_orien = malloc(num_data * sizeof(int *)) ;
    high_prob = malloc(num_data * sizeof(double *)) ;
    num_orien = malloc(num_data * sizeof(int)) ;
    
    idx = 0 ;
    for (rank = 0 ; rank < nproc ; rank++){
        sprintf(infile, "%s/high_p-%d.dat", prob_dir, rank) ;
        fp = fopen(infile, "r") ;
        fscanf(fp, "%d", &num_frame) ;
        for (d = 0 ; d < num_frame ; d++){
            fscanf(fp, "%d", &num_orien[idx]) ;
            prob_orien[idx] = malloc(num_orien[idx] * sizeof(int)) ;
            high_prob[idx] = malloc(num_orien[idx] * sizeof(double)) ;
            ct = 0 ;
            for (i = 0 ; i < num_orien[idx] ; i++){
                fscanf(fp, "%d %lf", &rid, &pval) ;
                if (pval < P_MIN)
                    continue ;
                prob_orien[idx][ct] = rid ;
                high_prob[idx][ct] = pval ;
                ct += 1 ;
            }
            num_orien[idx] = ct ;
            idx += 1 ;
        }
        fclose(fp) ;
    }

    printf("num_data = %d\n", idx) ;

    if_coarse_visit = calloc(num_coarse_rot, sizeof(int)) ;
    if_visit = calloc(num_rot, sizeof(int)) ;
    for (d = 0 ; d < num_data ; d++){
        for (i = 0 ; i < num_orien[d] ; i++){
            rid = prob_orien[d][i] ;
            if_coarse_visit[rid] = 1 ;
            for (j = 0 ; j < num_quat_neighbor[rid] ; j++)
                if_visit[quat_map[rid][j]] = 1 ;
        }
    }
}


void free_mem(){

    int d, r ;
    free(coarse_quat) ;
    free(quat) ;
    for (r = 0 ; r < num_coarse_rot ; r++)
        free(quat_map[r]) ;
    free(quat_map) ;
    free(num_quat_neighbor) ;
    for (d = 0 ; d < num_data ; d++){
        free(prob_orien[d]) ;
        free(high_prob[d]) ;
    }
    free(prob_orien) ;
    free(high_prob) ;
    free(num_orien) ;
    free(if_coarse_visit) ;
    free(if_visit) ;
}
