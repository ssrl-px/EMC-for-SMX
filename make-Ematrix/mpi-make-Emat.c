/*
write maps between peaks and rotations

compile:
mpicc mpi-make-Emat.c -O3 -lm -o emat

usage:
mpirun -np [number of processors] ./emat (-low) [path of the config_file] >> run.log &
"-low": save calculation on rotation samples using prob_orien.bin for low-resolution reconstruction

needs:
config_file, quat_file, rec-vectors.bin, rec2pix-map.dat, pix-map.dat, basis-vec.dat 

makes:
r2peak_file, peak2r_file

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

typedef struct{
    int row_min ;
    int row_max ;
    int col_min ;
    int col_max ;
} box ;

const int root = 0 ;
const int mem_step_size = 10000 ;

char config_file[256], home_dir[256], quat_file[256] ;
char r2peak_file[256], peak2r_file[256] ;
int nproc, myid, num_row, num_col, num_rot, num_pix, if_low ;
int hmax, kmax, lmax, hlen, klen, llen, num_hkl, total_pix, num_load ;
int *rec2pix_map, *pix_map, *num_peak2r, **peak2r, *load2rank, *if_skip_r ;
int *hkl_map, **r2peak, *cur_r2peak, *num_r2peak, *if_visit ;
double (*quat)[5], (*pix)[3], beam_vec[3], rot[3][3] ;
double basis_vec[3][3], rlatt_vec[3][3], proj[3][3] ;
double tol, min_rcell_sz, rmax2, rmin2 ;
box *peak_r ;

int setup() ;
void free_mem() ;
void make_Emat( int ) ;
void write_Emat() ;
void make_rot( int ) ;

int main(int argc, char* argv[]){

    MPI_Init(&argc, &argv) ;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc) ;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid) ;
    MPI_Win win ;

    int i, load_id, unity = 1 ;
    double t1, t2 ;
    t1 = MPI_Wtime() ;

    if_low = 0 ;
    if (argc == 3){
        if (strstr(argv[1], "-low"))
            if_low = 1 ;
        sprintf(config_file, "%s", argv[2]) ;
    }
    else if (argc == 2)
        sprintf(config_file, "%s", argv[1]) ;
    else
        exit(1) ;

    if (myid == 0){
        if (if_low == 1)
            printf("make Ematrix for low-resolution reconstruction\n") ;
        else
            printf("make Ematrix for local update scheme\n") ;
    }

    if (!setup()){
        printf("setup fails!!") ;
        exit(1) ;
    }

    MPI_Barrier(MPI_COMM_WORLD) ;
    t2 = MPI_Wtime() ;
    if (myid == 0)
        printf("setup takes: %.0f sec\n", t2-t1) ;

    if (num_rot % mem_step_size == 0)
        num_load = num_rot / mem_step_size ;
    else
        num_load = num_rot / mem_step_size + 1 ;

    load2rank = calloc(num_load, sizeof(int)) ;
    if (myid == 0)
        printf("num_load = %d\n", num_load) ;

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
            
            make_Emat( load_id ) ;
            MPI_Win_lock(MPI_LOCK_EXCLUSIVE, root, 0, win) ;
            MPI_Get(&load_id, 1, MPI_INT, root, 0, 1, MPI_INT, win) ;
            MPI_Accumulate(&unity, 1, MPI_INT, root, 0, 1, MPI_INT, MPI_SUM, win) ;
            MPI_Win_unlock(root, win) ;
            
            if (load_id > num_load-1)
                break ;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD) ;
    write_Emat() ;
    free_mem() ;

    MPI_Win_free(&win) ;
    MPI_Barrier(MPI_COMM_WORLD) ;
    t2 = MPI_Wtime() ;
    if (myid == 0)
        printf("elapsed time: %.0f sec\n", t2-t1) ;

    MPI_Finalize() ;
    return 0 ;
}


void make_Emat( int load_id ){
    
    int i, j, r, t, rmin, rmax, peak_ct ;
    int hval, kval, lval, klen_llen, hkl_id ;
    int idx, row_id, col_id, length ;
    double dx, dy, dz, rot_pix[3], coefft[3] ;

    klen_llen = klen*llen ;
    rmin = load_id*mem_step_size ;
    if (rmin > num_rot - 1)
        return ;
    rmax = rmin + mem_step_size ;
    if (rmax > num_rot)
        rmax = num_rot ;

    load2rank[load_id] = myid ;
    printf("myid = %d: load_id = %d, (rmin, rmax) = (%d, %d)\n", \
        myid, load_id, rmin, rmax) ;
    
    for (r = rmin ; r < rmax ; r++){
        
        make_rot(r) ;
        peak_ct = 0 ;

        if (if_skip_r[r] == 1){
            num_r2peak[r] = peak_ct ;
            r2peak[r] = malloc((1 + 3*peak_ct) * sizeof(int)) ;
            r2peak[r][0] = peak_ct ;
            continue ;
        }

        for (t = 0 ; t < num_pix ; t++){
            for (i = 0 ; i < 3 ; i++){
                rot_pix[i] = rot[i][0]*pix[t][0] ;
                for (j = 1 ; j < 3 ; j++)
                    rot_pix[i] += rot[i][j]*pix[t][j] ;
            }

            for (i = 0 ; i < 3 ; i++){
                coefft[i] = proj[i][0]*rot_pix[0] ;
                for (j = 1 ; j < 3 ; j++)
                    coefft[i] += proj[i][j]*rot_pix[j] ;
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
            
            idx = (hval+hmax)*klen_llen + (kval+kmax)*llen + (lval+lmax) ;
            hkl_id = hkl_map[idx] ;
            if (hkl_id < 0)
                continue ;
            
            dx = rot_pix[0] - hval*rlatt_vec[0][0] - kval*rlatt_vec[0][1] - lval*rlatt_vec[0][2] ;
            dy = rot_pix[1] - hval*rlatt_vec[1][0] - kval*rlatt_vec[1][1] - lval*rlatt_vec[1][2] ;
            dz = rot_pix[2] - hval*rlatt_vec[2][0] - kval*rlatt_vec[2][1] - lval*rlatt_vec[2][2] ;
            /* rot_pix is out of the view of its closet Bragg peak */
            if (dx*dx + dy*dy + dz*dz > tol)
                continue ;

            row_id = rec2pix_map[t] / num_col ;
            col_id = rec2pix_map[t] % num_col ;
            
            if (if_visit[hkl_id] == 0){
                cur_r2peak[peak_ct] = hkl_id ;
                peak_ct += 1 ;
                peak_r[hkl_id].row_min = row_id ;
                peak_r[hkl_id].row_max = row_id ;
                peak_r[hkl_id].col_min = col_id ;
                peak_r[hkl_id].col_max = col_id ;
                if_visit[hkl_id] = 1 ;
            }
            else{
                if (peak_r[hkl_id].row_min > row_id)
                    peak_r[hkl_id].row_min = row_id ;
                if (peak_r[hkl_id].row_max < row_id)
                    peak_r[hkl_id].row_max = row_id ;
                if (peak_r[hkl_id].col_min > col_id)
                    peak_r[hkl_id].col_min = col_id ;
                if (peak_r[hkl_id].col_max < col_id)
                    peak_r[hkl_id].col_max = col_id ;
            }
        }

        for (t = 0 ; t < peak_ct ; t++){
            hkl_id = cur_r2peak[t] ;
            num_peak2r[hkl_id] += 1 ;
            if_visit[hkl_id] = 0 ;
        }

        num_r2peak[r] = peak_ct ;
        r2peak[r] = malloc((1 + 3*peak_ct) * sizeof(int)) ;
        r2peak[r][0] = peak_ct ;
        for (t = 0 ; t < peak_ct ; t++){
            idx = 3*t + 1 ;
            hkl_id = cur_r2peak[t] ;
            r2peak[r][idx] = hkl_id ;
            r2peak[r][idx+1] = peak_r[hkl_id].row_min*num_col + peak_r[hkl_id].col_min ;
            r2peak[r][idx+2] = peak_r[hkl_id].row_max*num_col + peak_r[hkl_id].col_max ;
        }
    }
}


void write_Emat(){

    int i, j, r, t, rmin, rmax, length, idx, hkl_id ;
    long long *peak2r_offset ;

    MPI_File fh ;
    MPI_Status status ;
    MPI_Offset offset ;

    MPI_Allreduce(MPI_IN_PLACE, load2rank, num_load, MPI_INT, MPI_MAX, MPI_COMM_WORLD) ;
    MPI_Allreduce(MPI_IN_PLACE, num_r2peak, num_rot, MPI_INT, MPI_MAX, MPI_COMM_WORLD) ;

    offset = 0 ;
    MPI_File_open(MPI_COMM_WORLD, r2peak_file, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh) ;
    for (i = 0 ; i < num_load ; i++){
        rmin = i*mem_step_size ;
        if (rmin > num_rot - 1)
            break ;
        rmax = rmin + mem_step_size ;
        if (rmax > num_rot)
            rmax = num_rot ;
        
        for (r = rmin ; r < rmax ; r++){
            length = 1 + 3*num_r2peak[r] ;
            if (myid == load2rank[i])
                MPI_File_write_at(fh, offset, r2peak[r], length, MPI_INT, &status) ;
            offset += length*sizeof(int) ;
        }
        MPI_Barrier(MPI_COMM_WORLD) ;
    }
    MPI_File_close(&fh) ;

    peak2r = malloc(num_hkl * sizeof(int *)) ;
    for (i = 0 ; i < num_hkl ; i++)
        peak2r[i] = malloc(3*num_peak2r[i] * sizeof(int)) ;

    MPI_Allreduce(MPI_IN_PLACE, num_peak2r, num_hkl, MPI_INT, MPI_SUM, MPI_COMM_WORLD) ;
    MPI_File_open(MPI_COMM_WORLD, peak2r_file, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh) ;
    if (myid == 0){
        offset = 0 ;
        for (i = 0 ; i < num_hkl ; i++){
            MPI_File_write_at(fh, offset, &num_peak2r[i], 1, MPI_INT, &status) ;
            offset += (1 + 3*num_peak2r[i])*sizeof(int) ;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD) ;
    
    peak2r_offset = calloc(num_hkl, sizeof(long long)) ;
    peak2r_offset[0] = sizeof(int) ;
    for (i = 1 ; i < num_hkl ; i++)
        peak2r_offset[i] = peak2r_offset[i-1] + (1 + 3*num_peak2r[i-1])*sizeof(int) ;

    for (i = 0 ; i < num_load ; i++){
        rmin = i*mem_step_size ;
        if (rmin > num_rot - 1)
            break ;
        rmax = rmin + mem_step_size ;
        if (rmax > num_rot)
            rmax = num_rot ;

        memset(num_peak2r, 0, num_hkl*sizeof(int)) ;
        if (myid == load2rank[i]){
            for (r = rmin ; r < rmax ; r++){
                length = r2peak[r][0] ;
                for (t = 0 ; t < length ; t++){
                    idx = 3*t + 1 ;
                    hkl_id = r2peak[r][idx] ;
                    j = 3*num_peak2r[hkl_id] ;
                    peak2r[hkl_id][j] = r ;
                    peak2r[hkl_id][j+1] = r2peak[r][idx+1] ;
                    peak2r[hkl_id][j+2] = r2peak[r][idx+2] ;
                    num_peak2r[hkl_id] += 1 ;
                }
                free(r2peak[r]) ;
            }

            for (t = 0 ; t < num_hkl ; t++){
                offset = peak2r_offset[t] ;
                length = 3*num_peak2r[t] ;
                MPI_File_write_at(fh, offset, peak2r[t], length, MPI_INT, &status) ;
            }
        }

        MPI_Allreduce(MPI_IN_PLACE, num_peak2r, num_hkl, MPI_INT, MPI_MAX, MPI_COMM_WORLD) ;
        for (t = 0 ; t < num_hkl ; t++){
            length = 3*num_peak2r[t] ;
            peak2r_offset[t] += length*sizeof(int) ;
        }
    }
    MPI_File_close(&fh) ;    

    for (i = 0 ; i < num_hkl ; i++)
        free(peak2r[i]) ;
    free(peak2r) ;
    free(peak2r_offset) ;
}


int setup(){

    FILE *fp ;
    char *token, line[256], infile[256] ;
    char reduced_data_id_file[256], prob_orien_file[256] ;
    int i, j, k, d, r, t, x, y, z, i0, j0, i1, j1, qmax, qmin ;
    int hcur, kcur, lcur, max_num_hkl, idx, if_enclosed, num_data ;
    int tx, ty, tz, gw_ceil, VN, num, indices[3], *is_partial ;
    double detd, wl, px, cx, cy, res_cutoff, Rstop ;
    double vx, vy, vz, dx, dy, dz, dr2, gw, gw2 ;
    double norm, D, rescale, volume, polar, det ;
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
            else if (strcmp(token, "sx") == 0)
                beam_vec[0] = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "sy") == 0)
                beam_vec[1] = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "sz") == 0)
                beam_vec[2] = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "Rstop") == 0)
                Rstop = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "gw") == 0)
                gw = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "home_dir") == 0)
                strcpy(home_dir, strtok(NULL, " =\n")) ;
            
            if (if_low == 1){
                if (strcmp(token, "res_cutoff") == 0)
                    res_cutoff = atof(strtok(NULL, " =")) ;
                else if (strcmp(token, "quat_file") == 0)
                    strcpy(quat_file, strtok(NULL, " =\n")) ;
                else if (strcmp(token, "r2peak_file") == 0)
                    strcpy(r2peak_file, strtok(NULL, " =\n")) ;
                else if (strcmp(token, "peak2r_file") == 0)
                    strcpy(peak2r_file, strtok(NULL, " =\n")) ;
            }
            else{
                if (strcmp(token, "high_res_cutoff") == 0)
                    res_cutoff = atof(strtok(NULL, " =")) ;
                else if (strcmp(token, "fine_quat_file") == 0)
                    strcpy(quat_file, strtok(NULL, " =\n")) ;
                else if (strcmp(token, "local_r2peak_file") == 0)
                    strcpy(r2peak_file, strtok(NULL, " =\n")) ;
                else if (strcmp(token, "local_peak2r_file") == 0)
                    strcpy(peak2r_file, strtok(NULL, " =\n")) ;
            }
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
    qmax = ((int) ceil(wl*D/res_cutoff)) ;
    qmin = ((int) ceil(2*D*sin(atan(Rstop/D)/2))) ;
    rmax2 = pow(qmax*rescale, 2) ;
    rmin2 = pow(qmin*rescale, 2) ;
    
    sprintf(infile, "%s/make-detector/rec-vectors.bin", home_dir) ;
    fp = fopen(infile, "rb") ;
    fread(&qmax, sizeof(int), 1, fp) ;
    fread(&qmin, sizeof(int), 1, fp) ;
    fread(&num_pix, sizeof(int), 1, fp) ;
    pix = malloc(num_pix * sizeof(*pix)) ;
    for (t = 0 ; t < num_pix ; t++){
        norm = 0. ;
        for (i = 0 ; i < 3 ; i++){
            fread(&pix[t][i], sizeof(double), 1, fp) ;
            pix[t][i] *= rescale ;
            norm += pix[t][i]*pix[t][i] ;
        }
        fread(&polar, sizeof(double), 1, fp) ;
        if (norm > rmax2)
            break ;
    }
    num_pix = t ;
    pix = realloc(pix, num_pix * sizeof(*pix)) ;
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
    for (i = 0 ; i < total_pix ; i++){
        fscanf(fp, "%d", &pix_map[i]) ;
        if (pix_map[i] > num_pix-1)
            pix_map[i] = -1 ;
    }
    fclose(fp) ;

    sprintf(infile, "%s/aux/basis-vec.dat", home_dir) ;
    fp = fopen(infile, "r") ;
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

    if (myid == 0){
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
    min_rcell_sz = sqrt(tol) ;
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
    
    max_num_hkl = hlen*klen*llen ;
    is_partial = calloc(max_num_hkl, sizeof(int)) ;
    hkl_map = malloc(max_num_hkl * sizeof(int)) ;
    for (i = 0 ; i < max_num_hkl ; i++)
        hkl_map[i] = -1 ;

    rescale = VN/min_rcell_sz ;
    gw_ceil = ((int) ceil(gw)) ;
    gw2 = gw*gw ;

    num_hkl = 0 ;
    for (i = -hmax ; i < hmax+1 ; i++){
    for (j = -kmax ; j < kmax+1 ; j++){
    for (k = -lmax ; k < lmax+1 ; k++){
        
        indices[0] = i ;
        indices[1] = j ;
        indices[2] = k ;
        idx = (i+hmax)*klen*llen + (j+kmax)*llen + (k+lmax) ;

        for (i0 = 0 ; i0 < 3 ; i0++){
            rvec[i0] = 0. ;
            for (j0 = 0 ; j0 < 3 ; j0++)
                rvec[i0] += rlatt_vec[i0][j0]*indices[j0] ;
        }

        vx = rvec[0]*rescale ;
        vy = rvec[1]*rescale ;
        vz = rvec[2]*rescale ;

        tx = ((int) round(vx)) ;
        ty = ((int) round(vy)) ;
        tz = ((int) round(vz)) ;

        for ( x = -(gw_ceil+1) ; x < gw_ceil+2 ; x++){
        for ( y = -(gw_ceil+1) ; y < gw_ceil+2 ; y++){
        for ( z = -(gw_ceil+1) ; z < gw_ceil+2 ; z++){

            if (is_partial[idx] == 1)
                continue ;

            dx = tx + x - vx ;
            dy = ty + y - vy ;
            dz = tz + z - vz ;
            dr2 = dx*dx + dy*dy + dz*dz ;

            if (dr2 <= gw2){
                norm = (tx+x)*(tx+x) + (ty+y)*(ty+y) + (tz+z)*(tz+z) ;
                norm /= rescale*rescale ;
                if (norm > rmax2 || norm < rmin2){
                    is_partial[idx] = 1 ;
                    continue ;
                }
            }
        }
        }
        }

        if (is_partial[idx] == 0){
            hkl_map[idx] = num_hkl ;
            num_hkl += 1 ;
        }
    }
    }
    }

    if (myid == nproc-1)
        printf("hmax = %d, kmax = %d, lmax = %d, num_hkl = %d\n", hmax, kmax, lmax, num_hkl) ;

    r2peak = malloc(num_rot * sizeof(int *)) ;
    cur_r2peak = malloc(num_hkl * sizeof(int)) ;
    num_r2peak = calloc(num_rot, sizeof(int)) ;
    num_peak2r = calloc(num_hkl, sizeof(int)) ;
    peak_r = malloc(num_hkl * sizeof(box)) ;
    if_visit = calloc(num_hkl, sizeof(int)) ;
    if_skip_r = calloc(num_rot, sizeof(int)) ;

    if (if_low == 1){
        for (r = 0 ; r < num_rot ; r++)
            if_skip_r[r] = 1 ;
        if (myid == 0){
            sprintf(reduced_data_id_file, "%s/reduce-data/reduced-data_id.dat", home_dir) ;
            fp = fopen(reduced_data_id_file, "r") ;
            fscanf(fp, "%d", &num_data) ;
            fclose(fp) ;

            sprintf(prob_orien_file, "%s/aux/prob-orien.bin", home_dir) ;
            fp = fopen(prob_orien_file, "rb") ;
            for (d = 0 ; d < num_data ; d++){
                fread(&num, sizeof(int), 1, fp) ;
                for (i = 0 ; i < num ; i++){
                    fread(&r, sizeof(int), 1, fp) ;
                    if_skip_r[r] = 0 ;
                }
            }
            fclose(fp) ;

            num = 0 ;
            for (r = 0 ; r < num_rot ; r++){
                if (if_skip_r[r] == 0)
                    num += 1 ;
            }
            printf("number of rotations = %d\n", num_rot) ;
            printf("number of relevant rotations = %d\n", num) ;
        }
        MPI_Bcast(if_skip_r, num_rot, MPI_INT, 0, MPI_COMM_WORLD) ;
    }
    
    free(is_partial) ;
    return 1 ;
}


void free_mem(){
    
    free(quat) ;
    free(pix) ;
    free(rec2pix_map) ;
    free(pix_map) ;
    free(hkl_map) ;
    free(r2peak) ;
    free(cur_r2peak) ;
    free(num_r2peak) ;
    free(num_peak2r) ;
    free(peak_r) ;
    free(if_visit) ;
    free(if_skip_r) ;
    free(load2rank) ;
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
