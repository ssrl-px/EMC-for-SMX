/*
exhaustive search for the compatible orientations of the
identified peaks using one-sided communication in MPI

compile:
mpicc mpi-sync-orient-peak.c -O3 -lm -o orient.o

usage:
mpirun -np [number of processors] ./orient.o [path of the config_file] > run.log &

needs:
config_file, quat_file, basis-vec.dat, peakfiles.dat, orienfiles.dat

makes:
orien_files, num_prob_orien.dat

*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

const int root = 0 ;
const int mem_step_size = 100 ;

typedef struct{
    char name[128] ;
} filename ;

typedef struct{
    double vec[3] ;
} r_vector ;

r_vector **peak_rvec ;
filename *peakfiles, *orienfiles ;
char config_file[256], home_dir[256], quat_file[256] ;
char basis_file[256], peakfilelist[256], orienfilelist[256] ;
int nproc, myid, qmax, qmin, num_row, num_col, num_data, num_rot ;
int total_pix, min_patch_sz, max_patch_sz, min_num_peak, max_num_peak ;
int hmax, kmax, lmax, hlen, klen, llen, num_hkl, *hkl_map ;
int *num_peak, **prob_orien, *num_prob_orien, *max_num_prob_orien ;
double hot_pix_thres, tol, min_rcell_sz, rmin2, rmax2 ;
double (*quat)[5], (*pix)[3], (*rot_mat_table)[9] ;
double beam_vec[3], basis_vec[3][3], rot[3][3] ;
double rlatt_vec[3][3], proj[3][3] ;

int setup() ;
void free_mem() ;
void make_rot( int ) ;
void rotate_rvec( int ) ;

int main(int argc, char *argv[]){

	MPI_Init(&argc, &argv) ;
	MPI_Comm_size(MPI_COMM_WORLD, &nproc) ;
	MPI_Comm_rank(MPI_COMM_WORLD, &myid) ;
    MPI_Win win ;

    int i, d, ct, num_load, load_id, unity = 1 ;
    FILE *fp ;
    double t1, t2 ;
    t1 = MPI_Wtime() ;

    if (argc < 2){
        printf("The config_file is missing!!\n") ;
        exit(1) ;
    }
    else
        sprintf(config_file, "%s", argv[1]) ;

    if (!setup()){
        printf("setup fails!!\n") ;
        exit(1) ;
    }

    MPI_Barrier(MPI_COMM_WORLD) ;
    if (myid == 0)
        printf("setup takes %.0f sec\n", MPI_Wtime()-t1) ;
    
    if (num_data % mem_step_size == 0)
        num_load = num_data / mem_step_size ;
    else
        num_load = num_data / mem_step_size + 1 ;

    if (myid == root){
        load_id = nproc-1 ;
        MPI_Win_create(&load_id, sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &win) ;
    }
    else{
        load_id = myid-1 ;
        MPI_Win_create(NULL, 0, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &win) ;
    }

    if (myid != root){
        for (i = 0 ; i < num_load ; i++){

            rotate_rvec( load_id ) ;
            MPI_Win_lock(MPI_LOCK_EXCLUSIVE, root, 0, win) ;
            MPI_Get(&load_id, 1, MPI_INT, root, 0, 1, MPI_INT, win) ;
            MPI_Accumulate(&unity, 1, MPI_INT, root, 0, 1, MPI_INT, MPI_SUM, win) ;
            MPI_Win_unlock(root, win) ;
            
            if (load_id > num_load-1)
                break ;
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, num_prob_orien, num_data, MPI_INT, MPI_MAX, MPI_COMM_WORLD) ;
    MPI_Barrier(MPI_COMM_WORLD) ;
    
    if (myid == 0){
        fp = fopen("num_prob_orien.dat", "w") ;
        ct = 0 ;
        for (d = 0 ; d < num_data ; d++){
            fprintf(fp, "%d\n", num_prob_orien[d]) ;
            if (num_prob_orien[d] > 0)
                ct += 1 ;
        }
        fclose(fp) ;
        printf("number of relevant frames = %d\n", ct) ;
    }

    t2 = MPI_Wtime() ;
    if (myid == nproc-1)
        printf ("elapsed time: %.0f sec\n", t2-t1);

    free_mem() ;
    MPI_Win_free(&win) ;
    MPI_Barrier(MPI_COMM_WORLD) ;
    MPI_Finalize() ;
    return 0 ;
}


void rotate_rvec( int load_id ){

    int i, j, k, d, r, t, hval, kval, lval, max_len ;
    int dmin, dmax, patch_sz, pix_id, pix_val ;
    int *hkl_list, *if_visit, hkl_ct, hkl_id, klen_llen ;
    double rvec[3], weighted_vec[3], coefft[3] ;
    double dx, dy, dz, qval2, val_sum, inv_val_sum ;
    FILE *fp ;

    klen_llen = klen*llen ;
    
    /* dmin is inclusive */
    dmin = load_id * mem_step_size ;
    if (dmin > num_data-1)
        return ;

    /* dmax is exclusive */
    dmax = (load_id+1) * mem_step_size ;
    if (dmax > num_data)
        dmax = num_data ;

    printf("myid = %d: load_id = %d, (dmin, dmax) = (%d, %d)\n", myid, load_id, dmin, dmax) ;
    hkl_list = malloc(num_hkl * sizeof(int)) ;
    if_visit = calloc(num_hkl, sizeof(int)) ;

    for (d = dmin ; d < dmax ; d++){   
        
        num_peak[d] = 0 ;
        max_len = mem_step_size ;
        peak_rvec[d] = malloc(max_len * sizeof(r_vector)) ;

        fp = fopen(peakfiles[d].name, "r") ;
        while (1 == fscanf(fp, "%d", &patch_sz)){
            val_sum = 0. ;
            for (i = 0 ; i < 3 ; i++)
                weighted_vec[i] = 0. ;
            for (t = 0 ; t < patch_sz ; t++){
                fscanf(fp, "%d %d", &pix_id, &pix_val) ;
                for (i = 0 ; i < 3 ; i++)
                    weighted_vec[i] += pix[pix_id][i]*pix_val ;
                val_sum += pix_val ;
            }
            if (patch_sz < min_patch_sz)
                continue ;
            if (patch_sz > max_patch_sz)
                continue ;

            qval2 = 0. ;
            inv_val_sum = 1./val_sum ;
            for (i = 0 ; i < 3 ; i++){
                weighted_vec[i] *= inv_val_sum ;
                qval2 += weighted_vec[i]*weighted_vec[i] ;
            }

            if (qval2 > rmax2 || qval2 < rmin2)
                continue ;

            if (num_peak[d] == max_len){
                max_len += mem_step_size ;
                peak_rvec[d] = realloc(peak_rvec[d], max_len*sizeof(r_vector)) ;
            }

            for (i = 0 ; i < 3 ; i++)
                peak_rvec[d][num_peak[d]].vec[i] = weighted_vec[i] ;

            num_peak[d] += 1 ;
        }
        fclose(fp) ;

        /* No enough peaks in the frame */
        if (num_peak[d] < min_num_peak)
            continue ;
        
        /* too many peaks in the frame */
        if (num_peak[d] > max_num_peak){
            num_peak[d] = 0 ;
            continue ;
        }
        
        prob_orien[d] = malloc(mem_step_size * sizeof(int)) ;
        max_num_prob_orien[d] = mem_step_size ;
        num_prob_orien[d] = 0 ;
    }

    for (r = 0 ; r < num_rot ; r++){
        
        t = 0 ;
        for (i = 0 ; i < 3 ; i++){
            for (j = 0 ; j < 3 ; j++){
                rot[i][j] = rot_mat_table[r][t] ;
                t += 1 ;
            }
        }

        for (d = dmin ; d < dmax ; d++){
            /* No enough peaks in the frame */
            if (num_peak[d] < min_num_peak)
                continue ;
            
            hkl_ct = 0 ;
            for (t = 0 ; t < num_peak[d] ; t++){
                for (j = 0 ; j < 3 ; j++)
                    weighted_vec[j] = peak_rvec[d][t].vec[j] ;
                for (j = 0 ; j < 3 ; j++){
                    rvec[j] = rot[j][0]*weighted_vec[0] ;
                    for (k = 1 ; k < 3 ; k++)
                        rvec[j] += rot[j][k]*weighted_vec[k] ;
                }

                for (j = 0 ; j < 3 ; j++){
                    coefft[j] = 0. ;
                    for (k = 0 ; k < 3 ; k++)
                        coefft[j] += proj[j][k]*rvec[k] ;
                }

                hval = ((int) round(coefft[0])) ;
                if (abs(hval) > hmax)
                    continue ;
                kval = ((int) round(coefft[1])) ;
                if (abs(kval) > kmax)
                    continue ;
                lval = ((int) round(coefft[2])) ;
                if (abs(lval) > lmax)
                    continue ;
                
                dx = rvec[0] - hval*rlatt_vec[0][0] - kval*rlatt_vec[0][1] - lval*rlatt_vec[0][2] ;
                dy = rvec[1] - hval*rlatt_vec[1][0] - kval*rlatt_vec[1][1] - lval*rlatt_vec[1][2] ;
                dz = rvec[2] - hval*rlatt_vec[2][0] - kval*rlatt_vec[2][1] - lval*rlatt_vec[2][2] ;
                /* rvec is out of the view of its closet Bragg peak */
                if (dx*dx + dy*dy + dz*dz > tol)
                    continue ;
                
                hkl_id = hkl_map[(hval+hmax)*klen_llen + (kval+kmax)*llen + (lval+lmax)] ;
                if (hkl_id < 0)
                    continue ;
                if (if_visit[hkl_id] == 1)
                    continue ;
                else{
                    if_visit[hkl_id] = 1 ;
                    hkl_list[hkl_ct] = hkl_id ;
                    hkl_ct += 1 ;
                }

                if (hkl_ct == min_num_peak)
                    break ;
            }
            
            for (i = 0 ; i < hkl_ct ; i++)
                if_visit[hkl_list[i]] = 0 ;

            /* No enough peaks match the reciprocal lattice in this orientation */
            if (hkl_ct < min_num_peak)
                continue ;
            
            if (num_prob_orien[d] == max_num_prob_orien[d]){
                max_num_prob_orien[d] *= 2 ;
                prob_orien[d] = realloc(prob_orien[d], max_num_prob_orien[d] * sizeof(int)) ;
            }

            prob_orien[d][num_prob_orien[d]] = r ;
            num_prob_orien[d] += 1 ;
        }
    }

    for (d = dmin ; d < dmax ; d++){
        free(peak_rvec[d]) ;
        if (num_peak[d] < min_num_peak)
            continue ;
        if (num_prob_orien[d] > 0){
            fp = fopen(orienfiles[d].name, "w") ;
            fprintf(fp, "%d\n", num_prob_orien[d]) ;
            for (i = 0 ; i < num_prob_orien[d] ; i++)
                fprintf(fp, "%d\n", prob_orien[d][i]) ;
            fclose(fp) ;
        }
        free(prob_orien[d]) ;
    }

    free(hkl_list) ;
    free(if_visit) ;
}


int setup(){
    
    FILE *fp ;
    char *token, line[256] ;
    int i, j, k, d, r, t, i0, j0, i1, j1, VN, idx ;
    int hcur, kcur, lcur, if_enclosed, indices[3] ;
    double norm, x, y, sx, sy, sz, D, rescale, volume ;
    double detd, wl, px, cx, cy, res_cutoff, Rstop, gw, det ;
    double a[3], b[3], c[3], outer_prod[3] ;
    double rvec[3], mat[3][3], cofactor[3][3] ;
    
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
                num_data = atoi(strtok(NULL, " =")) ;
            else if (strcmp(token, "min_patch_sz") == 0)
                min_patch_sz = atoi(strtok(NULL, " =")) ;
            else if (strcmp(token, "max_patch_sz") == 0)
                max_patch_sz = atoi(strtok(NULL, " =")) ;
            else if (strcmp(token, "min_num_peak") == 0)
                min_num_peak = atoi(strtok(NULL, " =")) ;
            else if (strcmp(token, "max_num_peak") == 0)
                max_num_peak = atoi(strtok(NULL, " =")) ;
            else if (strcmp(token, "VN") == 0)
                VN = atoi(strtok(NULL, " =")) ;
            else if (strcmp(token, "detd") == 0)
                detd = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "wl") == 0)
                wl = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "px") == 0)
                px = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "cx") == 0)
                cx = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "cy") == 0)
                cy = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "hot_pix_thres") == 0)
                hot_pix_thres = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "sx") == 0)
                beam_vec[0] = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "sy") == 0)
                beam_vec[1] = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "sz") == 0)
                beam_vec[2] = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "Rstop") == 0)
                Rstop = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "res_cutoff") == 0)
                res_cutoff = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "gw") == 0)
                gw = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "home_dir") == 0)
                strcpy(home_dir, strtok(NULL, " =\n")) ;
            else if (strcmp(token, "quat_file") == 0)
                strcpy(quat_file, strtok(NULL, " =\n")) ;
        }
        fclose(fp) ;
    }

    fp = fopen(quat_file, "rb") ;
    fread(&num_rot, sizeof(int), 1, fp) ;
    quat = malloc(num_rot * sizeof(*quat)) ;
    for (r = 0 ; r < num_rot ; r++)
        fread(quat[r], sizeof(double), 5, fp) ;
    fclose(fp) ;

    norm = sqrt(pow(beam_vec[0], 2) + pow(beam_vec[1], 2) + pow(beam_vec[2], 2)) ;
    D = (detd/px) * norm/fabs(beam_vec[2]) ;
    rescale = 1./wl/D ;
    
    beam_vec[0] *= D/norm ;
    beam_vec[1] *= D/norm ;
    beam_vec[2] *= D/norm ;

	qmax = ((int) ceil(wl*D/res_cutoff)) ;
    qmin = ((int) ceil(2*D*sin(atan(Rstop/D)/2))) ;
    rmax2 = pow(qmax*rescale, 2) ;
    rmin2 = pow(qmin*rescale, 2) ;
    total_pix = num_row*num_col ;
    pix = malloc(total_pix * sizeof(*pix)) ;

    /* pix[t] has unit 1/A */
	for (t = 0 ; t < total_pix ; t++){
        x = t / num_col - cx ;
		y = t % num_col - cy ;
		sx = beam_vec[0] + x ;
		sy = beam_vec[1] + y ;
		sz = beam_vec[2] ;
		norm = sqrt(sx*sx + sy*sy + sz*sz) ;
		sx *= D/norm ;
		sy *= D/norm ;
		sz *= D/norm ;
        pix[t][0] = (sx - beam_vec[0])*rescale ;
        pix[t][1] = (sy - beam_vec[1])*rescale ;
        pix[t][2] = (sz - beam_vec[2])*rescale ;
	}

    sprintf(peakfilelist, "%s/make-background/peaklist.dat", home_dir) ;
    peakfiles = malloc(num_data * sizeof(filename)) ;
    fp = fopen(peakfilelist, "r") ;
    for (d = 0 ; d < num_data ; d++)
        fscanf(fp, "%s", peakfiles[d].name) ;
    fclose(fp) ;
  
    sprintf(orienfilelist, "%s/make-background/orienlist.dat", home_dir) ;
    orienfiles = malloc(num_data * sizeof(filename)) ;
    fp = fopen(orienfilelist, "r") ;
    for (d = 0 ; d < num_data ; d++)
        fscanf(fp, "%s", orienfiles[d].name) ;
    fclose(fp) ;

    peak_rvec = malloc(num_data * sizeof(r_vector *)) ;
    num_peak = malloc(num_data * sizeof(int)) ;
    prob_orien = malloc(num_data * sizeof(int *)) ;
    num_prob_orien = calloc(num_data, sizeof(int)) ;
    max_num_prob_orien = malloc(num_data * sizeof(int)) ;
    
    rot_mat_table = malloc(num_rot * sizeof(*rot_mat_table)) ;
    for (r = 0 ; r < num_rot ; r++){
        make_rot(r) ;
        t = 0 ;
        for (i = 0 ; i < 3 ; i++){
            for (j = 0 ; j < 3 ; j++){
                rot_mat_table[r][t] = rot[i][j] ;
                t += 1 ;
            }
        }
    }

    sprintf(basis_file, "%s/aux/basis-vec.dat", home_dir) ;
    fp = fopen(basis_file, "r") ;
    for (i = 0 ; i < 3 ; i++){
        for (j = 0 ; j < 3 ; j++)
            fscanf(fp, "%lf", &basis_vec[i][j]) ;
    }
    fclose(fp) ;

    /* construct reciprocal lattice vectors */
    for (i = 0 ; i < 3 ; i++){
        a[i] = basis_vec[i][0] ;
        b[i] = basis_vec[i][1] ;
        c[i] = basis_vec[i][2] ;
    }

    outer_prod[0] = a[1]*b[2] - a[2]*b[1] ;
    outer_prod[1] = a[2]*b[0] - a[0]*b[2] ;
    outer_prod[2] = a[0]*b[1] - a[1]*b[0] ;
    
    volume = 0. ;
    for (i = 0 ; i < 3 ; i++)
        volume += outer_prod[i]*c[i] ;
    volume = fabs(volume) ;
    for (i = 0 ; i < 3 ; i++)
        rlatt_vec[i][2] = outer_prod[i] / volume ;

    outer_prod[0] = b[1]*c[2] - b[2]*c[1] ;
    outer_prod[1] = b[2]*c[0] - b[0]*c[2] ;
    outer_prod[2] = b[0]*c[1] - b[1]*c[0] ;
    for (i = 0 ; i < 3 ; i++)
        rlatt_vec[i][0] = outer_prod[i] / volume ;

    outer_prod[0] = c[1]*a[2] - c[2]*a[1] ;
    outer_prod[1] = c[2]*a[0] - c[0]*a[2] ;
    outer_prod[2] = c[0]*a[1] - c[1]*a[0] ;
    for (i = 0 ; i < 3 ; i++)
        rlatt_vec[i][1] = outer_prod[i] / volume ;

    if (myid == nproc-1){
        for (i = 0 ; i < 3 ; i++){
            printf("%d-th reciprocal lattice vector:\n\t", i) ;
            for (j = 0 ; j < 3 ; j++)
                printf("%1.3e ", rlatt_vec[j][i]) ;
            printf("\n\n") ;
        }
    }

    norm = 0. ;
    for (j = 0 ; j < 3 ; j++){
        norm = 0. ;
        for (i = 0 ; i < 3 ; i++)
            norm += pow(rlatt_vec[i][j], 2) ;
        if (j == 0)
            tol = norm ;
        else{
            if (tol > norm)
                tol = norm ;
        }
    }
    tol *= pow(gw/VN, 2) ;

    /* calculate cofactors of rlatt_vec */
    for (i = 0 ; i < 3 ; i++){
        for (j = 0 ; j < 3 ; j++){
            if ((i+j) % 2 == 0)
                cofactor[i][j] = 1 ;
            else
                cofactor[i][j] = -1 ;
            
            i1 = 0 ;
            for (i0 = 0 ; i0 < 3 ; i0++){
                if (i0 == i)
                    continue ;
                j1 = 0 ;
                for (j0 = 0 ; j0 < 3 ; j0++){
                    if (j0 == j)
                        continue ;
                    mat[i1][j1] = rlatt_vec[i0][j0] ;
                    j1 += 1 ;
                }
                i1 += 1 ;
            }

            det = mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0] ;
            cofactor[i][j] *= det ;
        }
    }

    /* determinant of rlatt_vec */
    det = 0. ;
    for (i = 0 ; i < 3 ; i++)
        det += rlatt_vec[0][i] * cofactor[0][i] ;

    /* inverse of rlatt_vec */
    for (i = 0 ; i < 3 ; i++){
        for (j = 0 ; j < 3 ; j++)
            proj[i][j] = cofactor[j][i] / det ;
    }

    /* determine hmax, kmax, lmax */
    norm = 0. ;
    for (i = 0 ; i < 3 ; i++)
        norm += pow(rlatt_vec[i][0], 2) ;
    hcur = ((int) floor(sqrt(rmax2/norm))) ;

    norm = 0. ;
    for (i = 0 ; i < 3 ; i++)
        norm += pow(rlatt_vec[i][1], 2) ;
    kcur = ((int) floor(sqrt(rmax2/norm))) ;

    norm = 0. ;
    for (i = 0 ; i < 3 ; i++)
        norm += pow(rlatt_vec[i][2], 2) ;
    lcur = ((int) floor(sqrt(rmax2/norm))) ;

    if_enclosed = 1 ;
    while (if_enclosed == 1){
        
        if_enclosed = 0 ;
        for (i = 0 ; i < 2 ; i++){
            indices[0] = (2*i-1)*hcur ;
            for (j = -kcur ; j < kcur+1 ; j++){
                indices[1] = j ;
                for (k = -lcur ; k < lcur+1 ; k++){
                    indices[2] = k ;
                    if (if_enclosed == 1)
                        break ;
                    norm = 0. ;
                    for (i0 = 0 ; i0 < 3 ; i0++){
                        rvec[i0] = 0. ;
                        for (j0 = 0 ; j0 < 3 ; j0++)
                            rvec[i0] += rlatt_vec[i0][j0]*indices[j0] ;
                        norm += rvec[i0]*rvec[i0] ;
                    }
                    if (norm < rmax2)
                        if_enclosed = 1 ;
                }
            }
        }

        if (if_enclosed == 1){
            hcur += 1 ;
            continue ;
        }

        for (i = -hcur ; i < hcur+1 ; i++){
            indices[0] = i ;
            for (j = 0 ; j < 2 ; j++){
                indices[1] = (2*j-1)*kcur ;
                for (k = -lcur ; k < lcur+1 ; k++){
                    indices[2] = k ;
                    if (if_enclosed == 1)
                        break ;
                    norm = 0. ;
                    for (i0 = 0 ; i0 < 3 ; i0++){
                        rvec[i0] = 0. ;
                        for (j0 = 0 ; j0 < 3 ; j0++)
                            rvec[i0] += rlatt_vec[i0][j0]*indices[j0] ;
                        norm += rvec[i0]*rvec[i0] ;
                    }
                    if (norm < rmax2)
                        if_enclosed = 1 ;
                }
            }
        }
        
        if (if_enclosed == 1){
            kcur += 1 ;
            continue ;
        }

        for (i = -hcur ; i < hcur+1 ; i++){
            indices[0] = i ;
            for (j = -kcur ; j < kcur+1 ; j++){
                indices[1] = j ;
                for (k = 0 ; k < 2 ; k++){
                    indices[2] = (2*k-1)*lcur ;
                    if (if_enclosed == 1)
                        break ;
                    norm = 0. ;
                    for (i0 = 0 ; i0 < 3 ; i0++){
                        rvec[i0] = 0. ;
                        for (j0 = 0 ; j0 < 3 ; j0++)
                            rvec[i0] += rlatt_vec[i0][j0]*indices[j0] ;
                        norm += rvec[i0]*rvec[i0] ;
                    }
                    if (norm < rmax2)
                        if_enclosed = 1 ;
                }
            }
        }
        
        if (if_enclosed == 1){
            lcur += 1 ;
            continue ;
        }
    }

    hmax = hcur ;
    kmax = kcur ;
    lmax = lcur ;
    hlen = 2*hmax + 1 ;
    klen = 2*kmax + 1 ;
    llen = 2*lmax + 1 ;
    hkl_map = malloc(hlen*klen*llen * sizeof(int)) ;

    num_hkl = 0 ;
    for (i = -hmax ; i < hmax+1 ; i++){
    for (j = -kmax ; j < kmax+1 ; j++){
    for (k = -lmax ; k < lmax+1 ; k++){
        indices[0] = i ;
        indices[1] = j ;
        indices[2] = k ;
        norm = 0. ;
        for (i0 = 0 ; i0 < 3 ; i0++){
            rvec[i0] = 0. ;
            for (j0 = 0 ; j0 < 3 ; j0++)
                rvec[i0] += rlatt_vec[i0][j0]*indices[j0] ;
            norm += rvec[i0]*rvec[i0] ;
        }
        idx = (i+hmax)*klen*llen + (j+kmax)*llen + (k+lmax) ;
        if (norm > rmax2 || norm < rmin2){
            hkl_map[idx] = -1 ;
            continue ;
        }
        hkl_map[idx] = num_hkl ;
        num_hkl += 1 ;
    }
    }
    }

    if (myid == nproc-1)
        printf("hmax = %d, kmax = %d, lmax = %d, num_hkl = %d\n", hmax, kmax, lmax, num_hkl) ;

    return 1 ;
}


void free_mem(){

    free(quat) ;
    free(pix) ;
    free(peakfiles) ;
    free(orienfiles) ;
    free(peak_rvec) ;
    free(num_peak) ;
    free(prob_orien) ;
    free(num_prob_orien) ;
    free(max_num_prob_orien) ;
    free(rot_mat_table) ;
    free(hkl_map) ;
}


void make_rot( int r ){

    double q0, q1, q2, q3, q01, q02, q03, q11, q12, q13, q22, q23, q33 ;
    
    q0 = quat[r][0] ;
    q1 = quat[r][1] ;
    q2 = quat[r][2] ;
    q3 = quat[r][3] ;
    
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
}
