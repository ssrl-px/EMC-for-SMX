/*
embed the reconstructed low-resolution 3D intensity map into the 
high-resolution 3D intensity map to be used in the local update scheme

compile:
gcc init-intens.c -O3 -lm -o init.o

usage:
./init.o [path of the config_file] [path to the directory storing the low-resolution 3D intensity map]

needs:
config_file, basis-vec.dat, hkl-val.dat, finish_intensity.bin

makes:
start_intensity.bin

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

char config_file[256], home_dir[256], prob_dir[256] ;
int VN, hmax, kmax, lmax, hlen, klen, llen, num_hkl ;
int Nx, Ny, Nz, box_len, box_len2, box_half_len ;
int *hkl_table, *num_hkl2vox_id, **box_vox_id ;
double rmin2, rmax2, wl_D, gw, res_cutoff, min_rcell_sz ;
double rlatt_vec[3][3], proj[3][3], **intens1, *vox_val ;
long long NyNz ;

void setup() ;
void embed_intens() ;
void free_mem() ;

int main(int argc, char *argv[]){

    sprintf(config_file, "%s", argv[1]) ;
    sprintf(prob_dir, "%s", argv[2]) ;

    setup() ;
    embed_intens() ;
    free_mem() ;
    
    return 0 ;
}


void embed_intens(){

    FILE *fp ;
    char infile[256], outfile[256] ;
    int i, j, low_hkl_id, low_num_hkl, hval, kval, lval, *low_hkl_table ;
    int nvox, low_hmax, low_kmax, low_lmax, low_hval, low_kval, low_lval ;
    double *low_hkl_val, w, tol = 1.e-5 ;

    sprintf(infile, "%s/hkl-val.dat", prob_dir) ;
    fp = fopen(infile, "r") ;
    if (!fp)
        printf("cannot open %s\n", infile) ;

    fscanf(fp, "%d %d %d %d", &low_hmax, &low_kmax, &low_lmax, &low_num_hkl) ;
    low_hkl_table = malloc(3*low_num_hkl * sizeof(int)) ;
    low_hkl_val = malloc(low_num_hkl * sizeof(double)) ;
    for (i = 0 ; i < low_num_hkl ; i++){
        fscanf(fp, "%d %d %d %lf", &low_hkl_table[3*i], &low_hkl_table[3*i+1], \
            &low_hkl_table[3*i+2], &low_hkl_val[i]) ;
    }
    fclose(fp) ;

    sprintf(infile, "%s/finish_intensity.bin", prob_dir) ;
    fp = fopen(infile, "rb") ;
    if (!fp)
        printf("cannot open %s\n", infile) ;

    low_hkl_id = 0 ;
    low_hval = low_hkl_table[3*low_hkl_id] ;
    low_kval = low_hkl_table[3*low_hkl_id+1] ;
    low_lval = low_hkl_table[3*low_hkl_id+2] ;

    for (i = 0 ; i < num_hkl ; i++){
        hval = hkl_table[3*i] ;
        kval = hkl_table[3*i+1] ;
        lval = hkl_table[3*i+2] ;
        nvox = num_hkl2vox_id[i] ;
        
        if (hval == low_hval && kval == low_kval && lval == low_lval){
            fread(vox_val, sizeof(double), nvox, fp) ;
            for (j = 0 ; j < nvox ; j++)
                intens1[i][j] = vox_val[j] ;

            w = 0. ;
            for (j = 0 ; j < nvox ; j++)
                w += vox_val[j] ;
            if (fabs(w - low_hkl_val[low_hkl_id]) > tol*fabs(w))
                printf("error: reflection values do not match!!\n") ;

            low_hkl_id += 1 ;
            if (low_hkl_id == low_num_hkl)
                break ;

            low_hval = low_hkl_table[3*low_hkl_id] ;
            low_kval = low_hkl_table[3*low_hkl_id+1] ;
            low_lval = low_hkl_table[3*low_hkl_id+2] ;
        }
    }
    fclose(fp) ;

    printf("low_hmax = %d, low_kmax = %d, low_lmax = %d, low_num_hkl = %d\n", \
        low_hmax, low_kmax, low_lmax, low_hkl_id) ;

    sprintf(outfile, "%s/aux/start_intensity.bin", home_dir) ;
    fp = fopen(outfile, "wb") ;
    for (i = 0 ; i < num_hkl ; i++)
        fwrite(intens1[i], sizeof(double), num_hkl2vox_id[i], fp) ;
    fclose(fp) ;
}


void setup(){

    FILE *fp ;
    char *token, line[256], infile[256] ;
    double norm, detd, px, wl, D, rescale, Rstop, rand_val ;
    double basis_vec[3][3], cofactor[3][3], mat[2][2], rvec[3] ;
    double a[3], b[3], c[3], outer_prod[3], beam_vec[3] ;
    double volume, det, gw2, qx, qy, qz, dx, dy, dz, dr2 ;
    int i, j, k, i0, j0, i1, j1, qmin, qmax, idx ;
    int x, y, z, tx, ty, tz, vx, vy, vz, hkl_id ;
    int hcur, kcur, lcur, if_enclosed, max_num_vox ;
    int max_hkl_len, Nx_half, Ny_half, Nz_half ;
    int *is_partial, *num_vox, gw_ceil, indices[3] ;
    long long vox_id, **tmp_hkl2vox_id, **hkl2vox_id ;

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
            if (strcmp(token, "VN") == 0)
                VN = atoi(strtok(NULL, " =")) ;
            else if (strcmp(token, "wl") == 0)
                wl = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "px") == 0)
                px = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "detd") == 0)
                detd = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "gw") == 0)
                gw = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "Rstop") == 0)
                Rstop = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "sx") == 0)
                beam_vec[0] = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "sy") == 0)
                beam_vec[1] = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "sz") == 0)
                beam_vec[2] = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "high_res_cutoff") == 0)
                res_cutoff = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "home_dir") == 0)
                strcpy(home_dir, strtok(NULL, " =\n")) ;
        }
    }
    fclose(fp) ;
    
    norm = sqrt(pow(beam_vec[0], 2) + pow(beam_vec[1], 2) + pow(beam_vec[2], 2)) ;
    D = (detd/px) * norm/fabs(beam_vec[2]) ;
    wl_D = wl*D ;
    rescale = 1./wl_D ;

    qmin = ((int) ceil(2*D*sin(atan(Rstop/D)/2))) ;
    rmin2 = pow(qmin*rescale, 2) ;
    qmax = ((int) ceil(wl*D/res_cutoff)) ;
    rmax2 = pow(qmax*rescale, 2) ;

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
    
    for (i = 0 ; i < 3 ; i++){
        printf("%d-th reciprocal lattice vector:\n\t", i) ;
        for (j = 0 ; j < 3 ; j++)
            printf("%1.3e ", rlatt_vec[j][i]) ;
        printf("\n\n") ;
    }

    norm = 0. ;
    for (j = 0 ; j < 3 ; j++){
        norm = 0. ;
        for (i = 0 ; i < 3 ; i++)
            norm += pow(rlatt_vec[i][j], 2) ;
        if (j == 0)
            min_rcell_sz = norm ;
        else{
            if (min_rcell_sz > norm)
                min_rcell_sz = norm ;
        }
    }
    min_rcell_sz = sqrt(min_rcell_sz) ;

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

    gw2 = gw*gw ;
    gw_ceil = ((int) ceil(gw)) ;
    max_num_vox = (2*gw_ceil+3)*(2*gw_ceil+3)*(2*gw_ceil+3) ;
    max_hkl_len = hlen*klen*llen ;

    vox_val = malloc(max_num_vox * sizeof(double)) ;
    is_partial = calloc(max_hkl_len, sizeof(int)) ;
    num_vox = calloc(max_hkl_len, sizeof(int)) ;
    tmp_hkl2vox_id = malloc(max_hkl_len * sizeof(long long *)) ;
    for (i = 0 ; i < max_hkl_len ; i++)
        tmp_hkl2vox_id[i] = malloc(max_num_vox * sizeof(long long)) ;

    Nx = 2*((int) ceil(sqrt(rmax2)/min_rcell_sz*VN)) + 1 ;
    Ny = Nx ;
    Nz = Nx ;
    NyNz = Ny*Nz ;
    Nx_half = (Nx-1)/2 ;
    Ny_half = (Ny-1)/2 ;
    Nz_half = (Nz-1)/2 ;
    rescale = VN/min_rcell_sz ;

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

        qx = rvec[0]*rescale ;
        qy = rvec[1]*rescale ;
        qz = rvec[2]*rescale ;

        tx = ((int) round(qx)) ;
        ty = ((int) round(qy)) ;
        tz = ((int) round(qz)) ;

        for ( x = -(gw_ceil+1) ; x < gw_ceil+2 ; x++){
        for ( y = -(gw_ceil+1) ; y < gw_ceil+2 ; y++){
        for ( z = -(gw_ceil+1) ; z < gw_ceil+2 ; z++){
            
            if (is_partial[idx] == 1)
                continue ;

            dx = tx + x - qx ;
            dy = ty + y - qy ;
            dz = tz + z - qz ;
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

        if (is_partial[idx] == 1)
            continue ;
        
        for (x = -(gw_ceil+1) ; x < gw_ceil+2 ; x++){
        for (y = -(gw_ceil+1) ; y < gw_ceil+2 ; y++){
        for (z = -(gw_ceil+1) ; z < gw_ceil+2 ; z++){
            
            dx = tx + x - qx ;
            dy = ty + y - qy ;
            dz = tz + z - qz ;
            dr2 = dx*dx + dy*dy + dz*dz ;

            if (dr2 <= gw2 && abs(tx+x) <= Nx_half && abs(ty+y) <= Ny_half && abs(tz+z) <= Nz_half){
                vx = tx + x + Nx_half ;
                vy = ty + y + Ny_half ;
                vz = tz + z + Nz_half ;
                vox_id = vx*NyNz + vy*Nz + vz ;
                tmp_hkl2vox_id[idx][num_vox[idx]] = vox_id ;
                num_vox[idx] += 1 ;
            }
        }
        }
        }

        num_hkl += 1 ;
    }
    }
    }

    printf("Nx, Ny, Nz = %d, %d, %d\n", Nx, Ny, Nz) ;
    printf("hmax = %d, kmax = %d, lmax = %d, num_hkl = %d\n", hmax, kmax, lmax, num_hkl) ;

    hkl_table = malloc(3*num_hkl * sizeof(int)) ;
    num_hkl2vox_id = malloc(num_hkl * sizeof(int)) ;
    hkl2vox_id = malloc(num_hkl * sizeof(long long *)) ;

    hkl_id = 0 ;
    for (i = 0 ; i < max_hkl_len ; i++){
        if (num_vox[i] > 0){
            num_hkl2vox_id[hkl_id] = num_vox[i] ;
            hkl2vox_id[hkl_id] = malloc(num_vox[i] * sizeof(long long)) ;
            for (j = 0 ; j < num_vox[i] ; j++)
                hkl2vox_id[hkl_id][j] = tmp_hkl2vox_id[i][j] ;
            x = i / (klen*llen) ;
            y = (i % (klen*llen)) / llen ;
            z = i % llen ;
            idx = 3*hkl_id ;
            hkl_table[idx] = x - hmax ;
            hkl_table[idx+1] = y - kmax ;
            hkl_table[idx+2] = z - lmax ;
            hkl_id += 1 ;
        }
    }

    box_len = 2*gw_ceil + 3 ;
    box_len2 = box_len*box_len ;
    box_half_len = (box_len - 1)/2 ;
    box_vox_id = malloc(num_hkl * sizeof(int *)) ;
    for (i = 0 ; i < num_hkl ; i++){
        box_vox_id[i] = malloc(max_num_vox * sizeof(int)) ;
        for (j = 0 ; j < max_num_vox ; j++)
            box_vox_id[i][j] = -1 ;
            
        idx = 3*i ;
        for (j = 0 ; j < 3 ; j++)
            indices[j] = hkl_table[idx+j] ;

        for (i0 = 0 ; i0 < 3 ; i0++){
            rvec[i0] = 0. ;
            for (j0 = 0 ; j0 < 3 ; j0++)
                rvec[i0] += rlatt_vec[i0][j0]*indices[j0] ;
        }

        qx = rvec[0]*rescale ;
        qy = rvec[1]*rescale ;
        qz = rvec[2]*rescale ;

        tx = ((int) round(qx)) + Nx_half ;
        ty = ((int) round(qy)) + Ny_half ;
        tz = ((int) round(qz)) + Nz_half ;

        for (j = 0 ; j < num_hkl2vox_id[i] ; j++){
            vox_id = hkl2vox_id[i][j] ;
            vx = vox_id / NyNz ;
            vy = (vox_id % NyNz) / Nz ;
            vz = vox_id % Nz ;
            x = vx - tx + box_half_len ;
            y = vy - ty + box_half_len ;
            z = vz - tz + box_half_len ;
            if (abs(x-box_half_len) > box_half_len || abs(y-box_half_len) > box_half_len || abs(z-box_half_len) > box_half_len){
                printf("error in constructing box_vox_id[%d]!!\n", i) ;
                break ;
            }
            box_vox_id[i][x*box_len2 + y*box_len + z] = j ;
        }
    }

    intens1 = malloc(num_hkl * sizeof(double *)) ;
    for (i = 0 ; i < num_hkl ; i++)
        intens1[i] = malloc(num_hkl2vox_id[i] * sizeof(double)) ;

    for (i = 0 ; i < num_hkl ; i++){
        idx = 3*i ;
        for (j = 0 ; j < 3 ; j++)
            indices[j] = hkl_table[idx+j] ;

        for (i0 = 0 ; i0 < 3 ; i0++){
            rvec[i0] = 0. ;
            for (j0 = 0 ; j0 < 3 ; j0++)
                rvec[i0] += rlatt_vec[i0][j0]*indices[j0] ;
        }

        qx = rvec[0]*rescale ;
        qy = rvec[1]*rescale ;
        qz = rvec[2]*rescale ;

        rand_val = ((double) rand()) / RAND_MAX ;
        for (j = 0 ; j < num_hkl2vox_id[i] ; j++){
            vox_id = hkl2vox_id[i][j] ;
            vx = vox_id / NyNz - Nx_half ;
            vy = (vox_id % NyNz) / Nz - Ny_half ;
            vz = vox_id % Nz - Nz_half ;
            dx = vx - qx ;
            dy = vy - qy ;
            dz = vz - qz ;
            dr2 = dx*dx + dy*dy + dz*dz ;
            intens1[i][j] = rand_val*exp(-dr2/gw2) ;
        }
    }
}


void free_mem(){

    int i ;

    free(vox_val) ;
    for (i = 0 ; i < num_hkl ; i++)
        free(intens1[i]) ;
    free(intens1) ;
}
