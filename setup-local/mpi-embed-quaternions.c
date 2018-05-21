/*
embed a coarser quaternion sampling into a finer quaternion sampling

compile:
mpicc mpi-embed-quaternions.c -O3 -lm -Wall -o embed.o

need:
c-quaternion%d.bin, c-quaternion%d.bin % (n_coarse, n_fine)

usage:
mpirun -np [number of processors] ./embed.o n_coarse n_fine

make:
quat_embedding-%d-%d.dat % (n_coarse, n_fine)

Written in May 2017
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

const int mem_step_sz = 100 ;
const double epsilon = 1.e-10 ;
const double PI = 3.14159265358979323 ;
int nproc, myid, n_coarse, n_fine, num_coarse_rot, num_fine_rot ;
int *arg_map, *sorted_arg_map, count_q0, count_q1, count_q2 ;
int *idx_q0, *idx_q1, *idx_q2, *idx_q0q1, *idx_q1q2 ;
double **fine_quat, **coarse_quat, **sorted_fine_quat ;
double *val_q0, *val_q1, *val_q2 ;

void setup() ;
void check_sort() ;
void explore_vicinity() ;
void free_mem() ;
void merge_sort( double **, double **, int, int *, int * ) ;
void TopDownSplitMerge( double **, int, int, double **, int *, int * ) ;
void TopDownMerge( double **, int, int, int, double **, int *, int * ) ;
void CopyArray( double **, int, int, double **, int *, int * ) ;

int main(int argc, char* argv[]){

    MPI_Init(&argc, &argv) ;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc) ;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid) ;

    double t1, t2 ;
    t1 = MPI_Wtime() ;

    if (argc == 3){
        n_coarse = atoi(argv[1]) ;
        n_fine = atoi(argv[2]) ;
    }
    else{
        printf("need two arguments: n_coarse and n_fine!!\n") ;
        return 1 ;
    }

    setup() ;

    /* sort the duos in ascending order */
    merge_sort(fine_quat, sorted_fine_quat, num_fine_rot, arg_map, sorted_arg_map) ;
    check_sort() ;

    /* generate reference table */
    explore_vicinity() ;

    t2 = MPI_Wtime() ;
    if (myid == 0)
        printf("Computation time = %.0f sec\n", t2-t1) ;

    free_mem() ;
    MPI_Barrier(MPI_COMM_WORLD) ;
    MPI_Finalize() ;
    return 0 ;
}


void explore_vicinity(){

    FILE *fp ;
    int r, i, j, k, t, rid, s_num_rot, rmin, rmax ;
    int i0, j0, j1, k0, k1, idx, idx_min, idx_max ;
    int iBegin, iMiddle, iEnd, if_find, ct, len ;
    int **quat_table, *quat_table_r, *table_ct ;
    double dot_product, delta, delta2, gamma, gamma2, dq2 ;
    double q0_min, q1_min, q2_min, q0_max, q1_max, q2_max ;
    double *coarse_quat_r, *fine_quat_r ;
    char outfile[128] ;

    /* radius of the vicinity, which is made slightly coarser than the
     * angular resolution when n = n_coarse to avoid missing samples */
    if (n_coarse == n_fine)
        delta = (0.944/n_coarse) * 1.01 ;
    else
        delta = (0.944/n_coarse/2.) * 1.01 ;

    delta2 = delta*delta ;
    if (myid == 0)
        printf("delta = %1.5e\n", delta) ;

    if (num_coarse_rot % nproc == 0)
        s_num_rot = num_coarse_rot / nproc ;
    else
        s_num_rot = num_coarse_rot / nproc + 1 ;

    table_ct = calloc(s_num_rot, sizeof(int)) ;
    quat_table = malloc(s_num_rot * sizeof(int *)) ;

    rmin = myid*s_num_rot ;
    if (myid == nproc-1)
        rmax = num_coarse_rot ;
    else
        rmax = rmin + s_num_rot ;

    for (r = rmin ; r < rmax ; r++){
        rid = r % s_num_rot ;
        quat_table[rid] = malloc(mem_step_sz * sizeof(int)) ;
        quat_table_r = quat_table[rid] ;
        coarse_quat_r = coarse_quat[r] ;
        ct = table_ct[rid] ;
        len = mem_step_sz ;

        q0_max = coarse_quat_r[0] + delta ;
        if (q0_max > 1.)
            q0_max = 1. ;

        q0_min = coarse_quat_r[0] - delta ;
        if (q0_min < val_q0[0]){
            i0 = 0 ;
            q0_min = val_q0[i0] ;
        }
        else if (q0_min >= val_q0[count_q0 - 1]){
            i0 = count_q0 - 1 ;
            q0_min = val_q0[i0] ;
        }
        else{
            iBegin = 0 ;
            iEnd = count_q0 - 1 ;
            iMiddle = (iBegin + iEnd)/2 ;
            if_find = 0 ;

            while (if_find == 0){
                if (val_q0[iMiddle] <= q0_min && q0_min < val_q0[iMiddle+1]){
                    i0 = iMiddle ;
                    break ;
                }

                if (iEnd - iBegin <= 1){
                    printf("r = %d: fail to locate q0_min\n", r) ;
                    return ;
                }

                if (val_q0[iMiddle] <= q0_min)
                    iBegin = iMiddle ;
                else if (val_q0[iMiddle] > q0_min)
                    iEnd = iMiddle ;
                iMiddle = (iBegin + iEnd)/2 ;
            }
        }

        for (i = i0 ; i < count_q0 ; i++){
            if (val_q0[i] - q0_max > epsilon)
                break ;

            gamma = coarse_quat_r[0] - val_q0[i] ;
            gamma2 = delta2 - gamma*gamma ;
            if (gamma2 < 0.)
                continue ;
            gamma = sqrt(gamma2) ;

            q1_max = coarse_quat_r[1] + gamma ;
            if (q1_max > 1.)
                q1_max = 1. ;

            if (i == count_q0 - 1)
                j1 = count_q1 ;
            else
                j1 = idx_q0q1[i+1] ;

            q1_min = coarse_quat_r[1] - gamma ;
            if (q1_min < val_q1[idx_q0q1[i]]){
                j0 = idx_q0q1[i] ;
                q1_min = val_q1[j0] ;
            }
            else if (q1_min >= val_q1[j1 - 1]){
                j0 = j1 - 1 ;
                q1_min = val_q1[j0] ;
            }
            else{
                iBegin = idx_q0q1[i] ;
                iEnd = j1 - 1 ;
                iMiddle = (iBegin + iEnd)/2 ;
                if_find = 0 ;

                while (if_find == 0){
                    if (val_q1[iMiddle] <= q1_min && q1_min < val_q1[iMiddle+1]){
                        j0 = iMiddle ;
                        break ;
                    }

                    if (iEnd - iBegin <= 1){
                        printf("r = %d: fail to locate q1_min\n", r) ;
                        return ;
                    }

                    if (val_q1[iMiddle] <= q1_min)
                        iBegin = iMiddle ;
                    else if (val_q1[iMiddle] > q1_min)
                        iEnd = iMiddle ;
                    iMiddle = (iBegin + iEnd)/2 ;
                }
            }

            for (j = j0 ; j < j1 ; j++){
                if (val_q1[j] - q1_max > epsilon)
                    break ;

                gamma = coarse_quat_r[0] - val_q0[i] ;
                gamma2 = delta2 - gamma*gamma ;
                gamma = coarse_quat_r[1] - val_q1[j] ;
                gamma2 -= gamma*gamma ;
                if (gamma2 < 0.)
                    continue ;
                gamma = sqrt(gamma2) ;

                q2_max = coarse_quat_r[2] + gamma ;
                if (q2_max > 1.)
                    q2_max = 1. ;

                if (j == count_q1 - 1)
                    k1 = count_q2 ;
                else
                    k1 = idx_q1q2[j+1] ;

                q2_min = coarse_quat_r[2] - gamma ;
                if (q2_min < val_q2[idx_q1q2[j]]){
                    k0 = idx_q1q2[j] ;
                    q2_min = val_q2[k0] ;
                }
                else if (q2_min >= val_q2[k1 - 1]){
                    k0 = k1 - 1 ;
                    q2_min = val_q2[k0] ;
                }
                else{
                    iBegin = idx_q1q2[j] ;
                    iEnd = k1 - 1 ;
                    iMiddle = (iBegin + iEnd)/2 ;
                    if_find = 0 ;

                    while (if_find == 0){
                        if (val_q2[iMiddle] <= q2_min && q2_min < val_q2[iMiddle+1]){
                            k0 = iMiddle ;
                            break ;
                        }

                        if (iEnd - iBegin <= 1){
                            printf("r = %d: fail to locate q2_min\n", r) ;
                            return ;
                        }

                        if (val_q2[iMiddle] <= q2_min)
                            iBegin = iMiddle ;
                        else if (val_q2[iMiddle] > q2_min)
                            iEnd = iMiddle ;
                        iMiddle = (iBegin + iEnd)/2 ;
                    }
                }

                for (k = k0 ; k < k1 ; k++){
                    if (val_q2[k] - q2_max > epsilon)
                        break ;

                    idx_min = idx_q2[k] ;
                    if (k == count_q2 - 1)
                        idx_max = num_fine_rot ;
                    else
                        idx_max = idx_q2[k+1] ;

                    for (idx = idx_min ; idx < idx_max ; idx++){
                        fine_quat_r = fine_quat[idx] ;
                        dot_product = 0. ;
                        for (t = 0 ; t < 4 ; t++)
                            dot_product += fine_quat_r[t]*coarse_quat_r[t] ;
                        if (dot_product < 0.)
                            dq2 = 2 + 2*dot_product ;
                        else
                            dq2 = 2 - 2*dot_product ;

                        if (dq2 < delta2){
                            if (ct == len){
                                len *= 2 ;
                                quat_table[rid] = realloc(quat_table[rid], len*sizeof(int)) ;
                                quat_table_r = quat_table[rid] ;
                            }
                            quat_table_r[ct] = arg_map[idx] ;
                            ct += 1 ;
                        }
                    }
                }
            }
        }

        table_ct[rid] = ct ;
        quat_table[rid] = realloc(quat_table[rid], ct*sizeof(int)) ;
    }

    MPI_Barrier(MPI_COMM_WORLD) ;
    sprintf(outfile, "quaternion-%d-%d.dat", n_coarse, n_fine) ;
    if (myid == 0){
        fp = fopen(outfile, "w") ;
        fprintf(fp, "%d\n\n", num_coarse_rot) ;
        fclose(fp) ;
    }

    for (i = 0 ; i < nproc ; i++){
        if (i == myid){
            fp = fopen(outfile, "a") ;
            for (r = rmin ; r < rmax ; r++){
                rid = r % s_num_rot ;
                quat_table_r = quat_table[rid] ;
                ct = table_ct[rid] ;
                fprintf(fp, "%d\n", ct) ;
                for (j = 0 ; j < ct ; j++)
                    fprintf(fp, "%d ", quat_table_r[j]) ;
                fprintf(fp, "\n") ;
                free(quat_table[rid]) ;
            }
            fclose(fp) ;
        }
        MPI_Barrier(MPI_COMM_WORLD) ;
    }

    free(table_ct) ;
    free(quat_table) ;
}


void check_sort(){

    int r, i, id_min, id_max, len_q0, len_q1, len_q2, flag, ct ;
    double cur_q0, cur_q1, cur_q2, *fine_quat_r ;

    /* free the auxiliary arrays used in merge_sort */
    for (r = 0 ; r < num_fine_rot ; r++)
        free(sorted_fine_quat[r]) ;
    free(sorted_fine_quat) ;
    free(sorted_arg_map) ;

    len_q0 = mem_step_sz ;
    len_q1 = mem_step_sz ;
    len_q2 = mem_step_sz ;
    idx_q0 = malloc(mem_step_sz * sizeof(int)) ;
    idx_q1 = malloc(mem_step_sz * sizeof(int)) ;
    idx_q2 = malloc(mem_step_sz * sizeof(int)) ;

    cur_q0 = fine_quat[0][0] ;
    cur_q1 = fine_quat[0][1] ;
    cur_q2 = fine_quat[0][2] ;
    idx_q0[0] = 0 ;
    idx_q1[0] = 0 ;
    idx_q2[0] = 0 ;
    count_q0 = 1 ;
    count_q1 = 1 ;
    count_q2 = 1 ;

    for (r = 1 ; r < num_fine_rot ; r++){
        fine_quat_r = fine_quat[r] ;
        if (fabs(fine_quat_r[0] - cur_q0) < epsilon){
            if (fabs(fine_quat_r[1] - cur_q1) < epsilon){
                if (fabs(fine_quat_r[2] - cur_q2) < epsilon)
                    continue ;
                else{
                    cur_q2 = fine_quat_r[2] ;
                    if (count_q2 == len_q2){
                        len_q2 *= 2 ;
                        idx_q2 = realloc(idx_q2, len_q2*sizeof(int)) ;
                    }
                    idx_q2[count_q2] = r ;
                    count_q2 += 1 ;
                }
            }
            else{
                cur_q1 = fine_quat_r[1] ;
                if (count_q1 == len_q1){
                    len_q1 *= 2 ;
                    idx_q1 = realloc(idx_q1, len_q1*sizeof(int)) ;
                }
                idx_q1[count_q1] = r ;
                count_q1 += 1 ;

                cur_q2 = fine_quat_r[2] ;
                if (count_q2 == len_q2){
                    len_q2 *= 2 ;
                    idx_q2 = realloc(idx_q2, len_q2*sizeof(int)) ;
                }
                idx_q2[count_q2] = r ;
                count_q2 += 1 ;
            }
        }
        else{
            cur_q0 = fine_quat_r[0] ;
            if (count_q0 == len_q0){
                len_q0 *= 2 ;
                idx_q0 = realloc(idx_q0, len_q0*sizeof(int)) ;
            }
            idx_q0[count_q0] = r ;
            count_q0 += 1 ;

            cur_q1 = fine_quat_r[1] ;
            if (count_q1 == len_q1){
                len_q1 *= 2 ;
                idx_q1 = realloc(idx_q1, len_q1*sizeof(int)) ;
            }
            idx_q1[count_q1] = r ;
            count_q1 += 1 ;

            cur_q2 = fine_quat_r[2] ;
            if (count_q2 == len_q2){
                len_q2 *= 2 ;
                idx_q2 = realloc(idx_q2, len_q2*sizeof(int)) ;
            }
            idx_q2[count_q2] = r ;
            count_q2 += 1 ;
        }
    }

    idx_q0 = realloc(idx_q0, count_q0*sizeof(int)) ;
    idx_q1 = realloc(idx_q1, count_q1*sizeof(int)) ;
    idx_q2 = realloc(idx_q2, count_q2*sizeof(int)) ;

    for (i = 0 ; i < count_q2 ; i++){
        id_min = idx_q2[i] ;
        if (i == count_q2 - 1)
            id_max = num_fine_rot - 1 ;
        else
            id_max = idx_q2[i+1] - 1 ;

        for (r = id_min ; r < id_max ; r++){
            if (fine_quat[r][2] - fine_quat[r+1][2] > epsilon)
                printf("error in sorting q2!!\n") ;
            if (fabs(fine_quat[r][1] - fine_quat[r+1][1]) > epsilon || \
                fabs(fine_quat[r][0] - fine_quat[r+1][0]) > epsilon){
                printf("q0 and q1 should be the same!!\n") ;
            }
        }

        if (i < count_q2 - 1){
            flag = 0 ;
            if (fine_quat[id_max][0] - fine_quat[id_max+1][0] > epsilon)
                flag = 1 ;
            else if (fabs(fine_quat[id_max][0] - fine_quat[id_max+1][0]) < epsilon){
                if (fine_quat[id_max][1] - fine_quat[id_max+1][1] > epsilon)
                    flag = 1 ;
            }

            if (flag == 1)
                printf("q0 and q1 should be ascending!!\n") ;
        }
    }

    if (myid == 0)
        printf("finish checking idx_q2\n") ;

    for (i = 0 ; i < count_q1 ; i++){
        id_min = idx_q1[i] ;
        if (i == count_q1 - 1)
            id_max = num_fine_rot - 1 ;
        else
            id_max = idx_q1[i+1] - 1 ;

        for (r = id_min ; r < id_max ; r++){
            if (fine_quat[r][1] - fine_quat[r+1][1] > epsilon)
                printf("error in sorting q1!!\n") ;
            if (fabs(fine_quat[r][0] - fine_quat[r+1][0]) > epsilon)
                printf("q0 should be the same!!\n") ;
        }

        if (i < count_q1 - 1){
            if (fine_quat[id_max][0] - fine_quat[id_max+1][0] > epsilon)
                printf("q0 should be ascending!!\n") ;
        }
    }

    if (myid == 0)
        printf("finish checking idx_q1\n") ;

    for (i = 0 ; i < count_q0 ; i++){
        id_min = idx_q0[i] ;
        if (i == count_q0 - 1)
            id_max = num_fine_rot - 1 ;
        else
            id_max = idx_q0[i+1] - 1 ;

        for (r = id_min ; r < id_max ; r++){
            if (fine_quat[r][0] - fine_quat[r+1][0] > epsilon)
                printf("error in sorting q0!!\n") ;
        }
    }

    if (myid == 0)
        printf("finish checking idx_q0\n") ;

    for (r = 1 ; r < num_fine_rot ; r++){
        flag = 0 ;
        if (fabs(fine_quat[r-1][0] - fine_quat[r][0]) < epsilon){
            if (fabs(fine_quat[r-1][1] - fine_quat[r][1]) < epsilon){
                if (fabs(fine_quat[r-1][2] - fine_quat[r][2]) < epsilon){
                    if (fabs(fine_quat[r-1][3] - fine_quat[r][3] > epsilon))
                        flag = 1 ;
                }
                else if (fine_quat[r-1][2] - fine_quat[r][2] > epsilon)
                    flag = 1 ;
            }
            else if (fine_quat[r-1][1] - fine_quat[r][1] > epsilon)
                flag = 1 ;
        }
        else if (fine_quat[r-1][0] - fine_quat[r][0] > epsilon)
            flag = 1 ;

        if (flag == 1)
            printf("error in sorting fine_quat!!\n") ;
    }

    val_q0 = malloc(count_q0 * sizeof(double)) ;
    val_q1 = malloc(count_q2 * sizeof(double)) ;
    val_q2 = malloc(count_q2 * sizeof(double)) ;

    for (i = 0 ; i < count_q0 ; i++)
        val_q0[i] = fine_quat[idx_q0[i]][0] ;
    for (i = 0 ; i < count_q1 ; i++)
        val_q1[i] = fine_quat[idx_q1[i]][1] ;
    for (i = 0 ; i < count_q2 ; i++)
        val_q2[i] = fine_quat[idx_q2[i]][2] ;

    /* mapping between the start of indices */
    idx_q0q1 = malloc(count_q0 * sizeof(int)) ;
    idx_q1q2 = malloc(count_q1 * sizeof(int)) ;

    idx_q0q1[0] = 0 ;
    ct = 1 ;
    cur_q0 = fine_quat[0][0] ;
    for (i = 1 ; i < count_q1 ; i++){
        r = idx_q1[i] ;
        if (fabs(fine_quat[r][0] - cur_q0) > epsilon){
            cur_q0 = fine_quat[r][0] ;
            idx_q0q1[ct] = i ;
            ct += 1 ;
        }
    }

    if (ct != count_q0){
        if (myid == 0)
            printf("error in idx_q0q1\n") ;
    }

    idx_q1q2[0] = 0 ;
    ct = 1 ;
    cur_q1 = fine_quat[0][1] ;
    for (i = 1 ; i < count_q2 ; i++){
        r = idx_q2[i] ;
        if (fabs(fine_quat[r][1] - cur_q1) > epsilon){
            cur_q1 = fine_quat[r][1] ;
            idx_q1q2[ct] = i ;
            ct += 1 ;
        }
    }

    if (ct != count_q1){
        if (myid == 0)
            printf("error in idx_q1q2\n") ;
    }
}


void setup(){

    int r, i ;
    FILE *fp ;
    char infile[128] ;
    double weight ;

    sprintf(infile, "c-quaternion%d.bin", n_fine) ;
    fp = fopen(infile, "rb") ;
    fread(&num_fine_rot, sizeof(int), 1, fp) ;
    fine_quat = malloc(num_fine_rot * sizeof(double *)) ;
    for (r = 0 ; r < num_fine_rot ; r++){
        fine_quat[r] = malloc(4 * sizeof(double)) ;
        fread(fine_quat[r], sizeof(double), 4, fp) ;
        fread(&weight, sizeof(double), 1, fp) ;
    }
    fclose(fp) ;

    sprintf(infile, "c-quaternion%d.bin", n_coarse) ;
    fp = fopen(infile, "rb") ;
    fread(&num_coarse_rot, sizeof(int), 1, fp) ;
    coarse_quat = malloc(num_coarse_rot * sizeof(double *)) ;
    for (r = 0 ; r < num_coarse_rot ; r++){
        coarse_quat[r] = malloc(4 * sizeof(double)) ;
        fread(coarse_quat[r], sizeof(double), 4, fp) ;
        fread(&weight, sizeof(double), 1, fp) ;
    }
    fclose(fp) ;

    arg_map = malloc(num_fine_rot * sizeof(int)) ;
    sorted_arg_map = malloc(num_fine_rot * sizeof(int)) ;
    for (r = 0 ; r < num_fine_rot ; r++){
        arg_map[r] = r ;
        sorted_arg_map[r] = r ;
    }

    double *fine_quat_r, *coarse_quat_r ;
    for (r = 0 ; r < num_fine_rot ; r++){
        fine_quat_r = fine_quat[r] ;
        if (fine_quat_r[0] < -epsilon){
            for (i = 0 ; i < 4 ; i++)
                fine_quat_r[i] *= -1 ;
        }
    }

    for (r = 0 ; r < num_coarse_rot ; r++){
        coarse_quat_r = coarse_quat[r] ;
        if (coarse_quat_r[0] < -epsilon){
            for (i = 0 ; i < 4 ; i++)
                coarse_quat_r[i] *= -1 ;
        }
    }

    sorted_fine_quat = malloc(num_fine_rot * sizeof(double *)) ;
    for (r = 0 ; r < num_fine_rot ; r++)
        sorted_fine_quat[r] = malloc(4 * sizeof(double)) ;
}


void free_mem(){

    int r ;

    for (r = 0 ; r < num_fine_rot ; r++)
        free(fine_quat[r]) ;
    free(fine_quat) ;

    for (r = 0 ; r < num_coarse_rot ; r++)
        free(coarse_quat[r]) ;
    free(coarse_quat) ;

    free(arg_map) ;
    free(idx_q0) ;
    free(idx_q1) ;
    free(idx_q2) ;
    free(val_q0) ;
    free(val_q1) ;
    free(val_q2) ;
    free(idx_q0q1) ;
    free(idx_q1q2) ;
}


void merge_sort( double **A, double **B, int length, int *arg, int *sorted_arg ){
    TopDownSplitMerge(A, 0, length, B, arg, sorted_arg) ;
}


void TopDownSplitMerge( double **A, int iBegin, int iEnd, double **B, int *arg, int *sorted_arg ){

    /* iBegin: inclusive, iEnd: exclusive */
    if (iEnd - iBegin == 1)
        return ;

    int iMiddle = (iEnd + iBegin) / 2 ;
    TopDownSplitMerge(A, iBegin, iMiddle, B, arg, sorted_arg) ;
    TopDownSplitMerge(A, iMiddle, iEnd, B, arg, sorted_arg) ;
    TopDownMerge(A, iBegin, iMiddle, iEnd, B, arg, sorted_arg) ;
    CopyArray(B, iBegin, iEnd, A, arg, sorted_arg) ;
}


void TopDownMerge( double **A, int iBegin, int iMiddle, int iEnd, double **B, int *arg, int *sorted_arg ){

    int k, idx, i0 = iBegin, j0 = iMiddle, flag ;

    for (idx = iBegin ; idx < iEnd ; idx++){
        flag = 0 ;
        if (i0 < iMiddle){
            if (j0 >= iEnd)
                flag = 1 ;
            else if (A[i0][0] < A[j0][0] - epsilon)
                flag = 1 ;
            else if (fabs(A[i0][0] - A[j0][0]) < epsilon){
                if (A[i0][1] < A[j0][1] - epsilon)
                    flag = 1 ;
                else if (fabs(A[i0][1] - A[j0][1]) < epsilon){
                    if (A[i0][2] < A[j0][2] - epsilon)
                        flag = 1 ;
                    else if (fabs(A[i0][2] - A[j0][2]) < epsilon){
                        if (A[i0][3] < A[j0][3] - epsilon)
                            flag = 1 ;
                    }
                }
            }
        }

        if (flag == 1){
            for (k = 0 ; k < 4 ; k++)
                B[idx][k] = A[i0][k] ;
            sorted_arg[idx] = arg[i0] ;
            i0 += 1 ;
        }
        else{
            for (k = 0 ; k < 4 ; k++)
                B[idx][k] = A[j0][k] ;
            sorted_arg[idx] = arg[j0] ;
            j0 += 1 ;
        }
    }
}


void CopyArray( double **B, int iBegin, int iEnd, double **A, int *arg, int *sorted_arg ){

    int i, j ;

    for (i = iBegin ; i < iEnd ; i++){
        for (j = 0 ; j < 4 ; j++)
            A[i][j] = B[i][j] ;
        arg[i] = sorted_arg[i] ;
    }
}
