/*
reject frames with no or multiple crystals

compile:
gcc reject-frames.c -O3 -lm -o rej

usage:
./rej [path of the config_file] [path to the directory storing probability files]

needs:
config_file, sym-op.dat

makes:
out-phi.dat

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

const double P_MIN = 0.05 ;
char config_file[256], home_dir[256], prob_dir[256], quat_file[256] ;
int nproc, num_rot, num_data, num_div, num_symop, *count_d2r, **idx_d2r ;
double delta2, (*quat)[5], (*symop)[4], **pval_d2r, *phi_val ;

int quat_similarity() ;
void setup() ;
void make_symop() ;
void free_mem() ;

int main(int argc, char *argv[]){

    FILE *fp ;
    int i, j, d, count_zero, count_multi, *if_reject ;
    double epsilon = 1.e-10 ;

    sprintf(config_file, "%s", argv[1]) ;
    sprintf(prob_dir, "%s", argv[2]) ;

    setup() ;
    make_symop() ;

    /* angular resolution for quat_similarity */
    delta2 = pow(0.944/num_div, 2) ;
    if_reject = calloc(num_data, sizeof(int)) ;

    for (d = 0 ; d < num_data ; d++){
        if (phi_val[d] < epsilon || count_d2r[d] == 0){
            if_reject[d] = 1 ;
            continue ;
        }

        for (i = 0 ; i < count_d2r[d] - 1 ; i++){
            if (pval_d2r[d][i] < P_MIN)
                continue ;
            for (j = i+1 ; j < count_d2r[d] ; j++){
                if (pval_d2r[d][j] < P_MIN)
                    continue ;
                if (if_reject[d] == 1)
                    break ;
                if (quat_similarity(idx_d2r[d][i], idx_d2r[d][j]) == 0){
                    if_reject[d] = 1 ;
                    printf("d = %d has multiple crystals!!\n", d) ;
                }
            }
        }
    }

    count_zero = 0 ;
    count_multi = 0 ;
    for (d = 0 ; d < num_data ; d++){
        if (if_reject[d] == 1){
            if (phi_val[d] < epsilon)
                count_zero += 1 ;
            else
                count_multi += 1 ;
        }
    }
    printf("\n%d out of %d frames are rejected:\n", count_zero + count_multi, num_data) ;
    printf("\t%d of them are blank\n", count_zero) ;
    printf("\t%d of them have multiple crystals\n", count_multi) ;

    fp = fopen("out-phi.dat", "w") ;
    for (d = 0 ; d < num_data ; d++){
        if (if_reject[d] == 1)
            phi_val[d] = 0. ;
        fprintf(fp, "%1.5e\n", phi_val[d]) ;
    }
    fclose(fp) ;

    free(if_reject) ;
    free_mem() ;
    return 0 ;
}


void setup(){

    FILE *fp ;
    char *token, line[256], infile[256], reduced_data_id_file[256] ;
    int i, d, r, idx, rank, num_frame ;

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
            else if (strcmp(token, "reduced_data_id_file") == 0)
                strcpy(reduced_data_id_file, strtok(NULL, " =\n")) ;
        }
    }
    fclose(fp) ;

    fp = fopen(quat_file, "rb") ;
    fread(&num_rot, sizeof(int), 1, fp) ;
    quat = malloc(num_rot * sizeof(*quat)) ;
    for (r = 0 ; r < num_rot ; r++)
        fread(quat[r], sizeof(double), 5, fp) ;
    fclose(fp) ;

    token = strtok(quat_file, "/") ;
    while (token != NULL){
        if (strstr(token, "quaternion") != NULL)
            sscanf(token, "%*[^0123456789]%d", &num_div) ;
        token = strtok(NULL, "/") ;
        if (token == NULL)
            break ;
    }

    fp = fopen(reduced_data_id_file, "r") ;
    fscanf(fp, "%d", &num_data) ;
    fclose(fp) ;
    
    phi_val = malloc(num_data * sizeof(double)) ;
    sprintf(infile, "%s/total-phi.dat", prob_dir) ;
    fp = fopen(infile, "r") ;
    for (d = 0 ; d < num_data ; d++)
        fscanf(fp, "%lf", &phi_val[d]) ;
    fclose(fp) ;

    count_d2r = malloc(num_data * sizeof(int)) ;
    idx_d2r = malloc(num_data * sizeof(int *)) ;
    pval_d2r = malloc(num_data * sizeof(double *)) ;

    idx = 0 ;
    for (rank = 0 ; rank < nproc ; rank++){
        sprintf(infile, "%s/high_p-%d.dat", prob_dir, rank) ;
        fp = fopen(infile, "r") ;
        fscanf(fp, "%d", &num_frame) ;
        for (d = 0 ; d < num_frame ; d++){
            fscanf(fp, "%d", &count_d2r[idx]) ;
            idx_d2r[idx] = malloc(count_d2r[idx] * sizeof(int)) ;
            pval_d2r[idx] = malloc(count_d2r[idx] * sizeof(double)) ;
            for (i = 0 ; i < count_d2r[idx] ; i++)
                fscanf(fp, "%d %lf", &idx_d2r[idx][i], &pval_d2r[idx][i]) ;
            idx += 1 ;
        }
        fclose(fp) ;
    }
}


void free_mem(){
    
    int d ;
    free(quat) ;
    free(phi_val) ;
    for (d = 0 ; d < num_data ; d++){
        free(idx_d2r[d]) ;
        free(pval_d2r[d]) ;
    }
    free(idx_d2r) ;
    free(pval_d2r) ;
    free(count_d2r) ;
}


void make_symop(){
    
    FILE *fp ;
    char infile[256] ;
    int i, j, k, ct ;
    double (*sym_mat)[3][3], mat[3][3], val, det, epsilon = 1.e-10 ;

    sprintf(infile, "%s/aux/sym-op.dat", home_dir) ;
    fp = fopen(infile, "r") ;
    fscanf(fp, "%d", &num_symop) ;

    /* ignore the operations that involve an inversion */
    ct = 0 ;
    sym_mat = malloc((num_symop/2) * sizeof(*sym_mat)) ;
    for (k = 0 ; k < num_symop ; k++){
        for (i = 0 ; i < 3 ; i++){
            for (j = 0 ; j < 3 ; j++)
                fscanf(fp, "%lf", &mat[i][j]) ;
        }

        det = mat[0][0]*mat[1][1]*mat[2][2] ; 
        det += mat[1][0]*mat[2][1]*mat[0][2] ;
        det += mat[2][0]*mat[0][1]*mat[1][2] ;
        det -= mat[0][2]*mat[1][1]*mat[2][0] ;
        det -= mat[1][2]*mat[2][1]*mat[0][0] ;
        det -= mat[2][2]*mat[0][1]*mat[1][0] ;
        if (det < 0)
            continue ;

        for (i = 0 ; i < 3 ; i++){
            for (j = 0 ; j < 3 ; j++)
                sym_mat[ct][i][j] = mat[i][j] ;
        }
        ct += 1 ;
    }
    num_symop /= 2 ;
    fclose(fp) ;

    symop = malloc(num_symop * sizeof(*symop)) ;
    for (k = 0 ; k < num_symop ; k++){
        val = 1. + sym_mat[k][0][0] + sym_mat[k][1][1] + sym_mat[k][2][2] ;
        if (fabs(val) > epsilon){
            symop[k][0] = 0.5*sqrt(val) ;
            symop[k][1] = -0.25*(sym_mat[k][2][1] - sym_mat[k][1][2])/symop[k][0] ;
            symop[k][2] = -0.25*(sym_mat[k][0][2] - sym_mat[k][2][0])/symop[k][0] ;
            symop[k][3] = -0.25*(sym_mat[k][1][0] - sym_mat[k][0][1])/symop[k][0] ;
            continue ;
        }
        val = 1. + sym_mat[k][0][0] - sym_mat[k][1][1] - sym_mat[k][2][2] ;
        if (fabs(val) > epsilon){
            symop[k][1] = 0.5*sqrt(val) ;
            symop[k][0] = -0.25*(sym_mat[k][2][1] - sym_mat[k][1][2])/symop[k][1] ;
            symop[k][2] = 0.25*(sym_mat[k][0][1] + sym_mat[k][1][0])/symop[k][1] ;
            symop[k][3] = 0.25*(sym_mat[k][0][2] + sym_mat[k][2][0])/symop[k][1] ;
            continue ;
        }
        val = 1. - sym_mat[k][0][0] + sym_mat[k][1][1] - sym_mat[k][2][2] ;
        if (fabs(val) > epsilon){
            symop[k][2] = 0.5*sqrt(val) ;
            symop[k][0] = -0.25*(sym_mat[k][0][2] - sym_mat[k][2][0])/symop[k][2] ;
            symop[k][1] = 0.25*(sym_mat[k][0][1] + sym_mat[k][1][0])/symop[k][2] ;
            symop[k][3] = 0.25*(sym_mat[k][1][2] + sym_mat[k][2][1])/symop[k][2] ;
            continue ;
        }
        val = 1. - sym_mat[k][0][0] - sym_mat[k][1][1] + sym_mat[k][2][2] ;
        symop[k][3] = 0.5*sqrt(val) ;
        symop[k][0] = -0.25*(sym_mat[k][1][0] - sym_mat[k][0][1])/symop[k][3] ;
        symop[k][1] = 0.25*(sym_mat[k][0][2] + sym_mat[k][2][0])/symop[k][3] ;
        symop[k][2] = 0.25*(sym_mat[k][1][2] + sym_mat[k][2][1])/symop[k][3] ;
    }

    double norm, rot[3][3] ;
    double q0, q1, q2, q3, q01, q02, q03, q11, q12, q13, q22, q23, q33 ;

    for (k = 0 ; k < num_symop ; k++){
        q0 = symop[k][0] ;
        q1 = symop[k][1] ;
        q2 = symop[k][2] ;
        q3 = symop[k][3] ;

        q01 = q0*q1 ;
        q02 = q0*q2 ;
        q03 = q0*q3 ;
        q11 = q1*q1 ;
        q12 = q1*q2 ;
        q13 = q1*q3 ;
        q22 = q2*q2 ;
        q23 = q2*q3 ;
        q33 = q3*q3 ;

        rot[0][0] = 1. - 2.*(q22 + q33) ;
        rot[0][1] = 2.*(q12 + q03) ;
        rot[0][2] = 2.*(q13 - q02) ;
        rot[1][0] = 2.*(q12 - q03) ;
        rot[1][1] = 1. - 2.*(q11 + q33) ;
        rot[1][2] = 2.*(q01 + q23) ;
        rot[2][0] = 2.*(q02 + q13) ;
        rot[2][1] = 2.*(q23 - q01) ;
        rot[2][2] = 1. - 2.*(q11 + q22) ;

        norm = 0. ;
        for (i = 0 ; i < 3 ; i++){
            for (j = 0 ; j < 3 ; j++)
                norm += pow(rot[i][j]-sym_mat[k][i][j], 2) ;
        }

        if (sqrt(norm) > epsilon)
            printf("error in make_symop!!\n") ;
    }

    free(sym_mat) ;
}


int quat_similarity( int r1, int r2 ){
    
    int i, k ;
    double q[4], p[4], r[4], dot_product, dq2 ;

    for (i = 0 ; i < 4 ; i++)
        q[i] = quat[r1][i] ;

    for (k = 0 ; k < num_symop ; k++){
        for (i = 0 ; i < 4 ; i++)
            p[i] = symop[k][i] ;

        /* transform q by the k-th symmetry operator */
        r[0] = q[0]*p[0] - q[1]*p[1] - q[2]*p[2] - q[3]*p[3] ;
        r[1] = q[0]*p[1] + p[0]*q[1] + q[2]*p[3] - q[3]*p[2] ;
        r[2] = q[0]*p[2] + p[0]*q[2] + q[3]*p[1] - q[1]*p[3] ;
        r[3] = q[0]*p[3] + p[0]*q[3] + q[1]*p[2] - q[2]*p[1] ;

        dot_product = 0. ;
        for (i = 0 ; i < 4 ; i++)
            dot_product += r[i]*quat[r2][i] ;

        if (dot_product < 0.)
            dq2 = 2 + 2*dot_product ;
        else
            dq2 = 2 - 2*dot_product ;

        if (dq2 < delta2)
            return 1 ;
    }

    return 0 ;
}
