#include "emc.h"

int setup(){
    
    FILE *fp ;
    char *token, line[256], infile[256], quat_file[256] ;
    char bg_file[256], r2peak_file[256], start_phi_file[256] ;
    char reduced_data_id_file[256], prob_orien_file[256] ;
    int i, j, d, r, t, num, f_id, idx, idx_start, rank_max ;
    int dmin, dmax, qmin, qmax, total_pix, length, r_block_sz ;
    double wl, px, detd, rescale, norm, D, val, dist, dq, Rstop ;
    double res_cutoff, beam_vec[3], quat_r[5], epsilon = 1.e-10 ;
    long long byte_offset ;

    MPI_File fh ;
    MPI_Status status ;
    MPI_Offset offset ;

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
            else if (strcmp(token, "qlen") == 0)
                qlen = atoi(strtok(NULL, " =")) ;
            else if (strcmp(token, "iter_data_block") == 0)
                iter_data_block = atoi(strtok(NULL, " =")) ;
            else if (strcmp(token, "wl") == 0)
                wl = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "px") == 0)
                px = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "detd") == 0)
                detd = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "gw") == 0)
                gw = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "res_cutoff") == 0)
                res_cutoff = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "Rstop") == 0)
                Rstop = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "sx") == 0)
                beam_vec[0] = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "sy") == 0)
                beam_vec[1] = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "sz") == 0)
                beam_vec[2] = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "home_dir") == 0)
                strcpy(home_dir, strtok(NULL, " =\n")) ;
            else if (strcmp(token, "prob_dir") == 0)
                strcpy(prob_dir, strtok(NULL, " =\n")) ;
            else if (strcmp(token, "quat_file") == 0)
                strcpy(quat_file, strtok(NULL, " =\n")) ;
            else if (strcmp(token, "mpi_datafile") == 0)
                strcpy(data_file, strtok(NULL, " =\n")) ;
            else if (strcmp(token, "mpi_bgfile") == 0)
                strcpy(bg_file, strtok(NULL, " =\n")) ;
            else if (strcmp(token, "r2peak_file") == 0)
                strcpy(r2peak_file, strtok(NULL, " =\n")) ;
            else if (strcmp(token, "peak2r_file") == 0)
                strcpy(peak2r_file, strtok(NULL, " =\n")) ;
            else if (strcmp(token, "start_phi_file") == 0)
                strcpy(start_phi_file, strtok(NULL, " =\n")) ;
            else if (strcmp(token, "start_intens_file") == 0)
                strcpy(start_intens_file, strtok(NULL, " =\n")) ;
            else if (strcmp(token, "reduced_data_id_file") == 0)
                strcpy(reduced_data_id_file, strtok(NULL, " =\n")) ;
            else if (strcmp(token, "prob_orien_file") == 0)
                strcpy(prob_orien_file, strtok(NULL, " =\n")) ;
        }
        fclose(fp) ;
    }

    fp = fopen(quat_file, "rb") ;
    fread(&num_rot, sizeof(int), 1, fp) ;
    quat = malloc(num_rot * sizeof(*quat)) ;
    rid_map = malloc(num_rot * sizeof(int)) ;
    for (r = 0 ; r < num_rot ; r++){
        fread(quat_r, sizeof(double), 5, fp) ;
        for (i = 0 ; i < 5 ; i++)
            quat[r][i] = quat_r[i] ;
    }
    fclose(fp) ;

    norm = sqrt(pow(beam_vec[0], 2) + pow(beam_vec[1], 2) + pow(beam_vec[2], 2)) ;
    D = (detd/px) * norm/fabs(beam_vec[2]) ;
    rescale = 1./wl/D ;
    qmax = ((int) ceil(wl*D/res_cutoff)) ;
    qmin = ((int) ceil(2*D*sin(atan(Rstop/D)/2))) ;
    rmax2 = pow(qmax*rescale, 2) ;
    rmin2 = pow(qmin*rescale, 2) ;
    wl_D = wl*D ;

    sprintf(infile, "%s/make-detector/rec-vectors.bin", home_dir) ;
    fp = fopen(infile, "rb") ;
    fread(&qmax, sizeof(int), 1, fp) ;
    fread(&qmin, sizeof(int), 1, fp) ;
    fread(&full_num_pix, sizeof(int), 1, fp) ;
    qval_max = qmax ;

    pix = malloc(full_num_pix * sizeof(*pix)) ;
    polar = malloc(full_num_pix * sizeof(float)) ;
    for (t = 0 ; t < full_num_pix ; t++){
        norm = 0. ;
        for (i = 0 ; i < 3 ; i++){
            fread(&val, sizeof(double), 1, fp) ;
            norm += pow(val*rescale, 2) ;
            pix[t][i] = val ;
        }
        fread(&val, sizeof(double), 1, fp) ;
        polar[t] = val ;
        if (norm > rmax2)
            break ;
    }
    num_pix = t ;
    pix = realloc(pix, num_pix * sizeof(*pix)) ;
    polar = realloc(polar, num_pix * sizeof(float)) ;
    fclose(fp) ;

    dq = qval_max / qlen ;
    inv_polar = malloc(num_pix * sizeof(float)) ;
    log_polar = malloc(num_pix * sizeof(float)) ;
    qid_map = malloc(num_pix * sizeof(int)) ;
    for (t = 0 ; t < num_pix ; t++){
        inv_polar[t] = 1./polar[t] ;
        log_polar[t] = log(polar[t]) ;
        dist = 0. ;
        for (i = 0 ; i < 3 ; i++)
            dist += pix[t][i]*pix[t][i] ;
        dist = sqrt(dist) ;
        idx = ((int) round(dist/dq - 0.5)) ;
        if (idx < 0 || idx > qlen-1)
            qid_map[t] = -1 ;
        else
            qid_map[t] = idx ;
    }

    total_pix = num_row*num_col ;
    pix_map = malloc(total_pix * sizeof(int)) ;
    sprintf(infile, "%s/make-detector/pix-map.dat", home_dir) ;
    fp = fopen(infile, "r") ;
    for (t = 0 ; t < total_pix ; t++){
        fscanf(fp, "%d", &pix_map[t]) ;
        if (pix_map[t] > num_pix-1)
            pix_map[t] = -1 ;
    }
    fclose(fp) ;

    fp = fopen(reduced_data_id_file, "r") ;
    fscanf(fp, "%d", &num_data) ;
    fclose(fp) ;

    if (num_data % nproc == 0)
        s_num_data = num_data / nproc ;
    else
        s_num_data = num_data / nproc + 1 ;

    rank_max = num_data - nproc*(s_num_data - 1) ;
    if (myid < rank_max){
        dmin = myid*s_num_data ;
        dmax = dmin + s_num_data ;
    }
    else{
        dmin = rank_max*s_num_data + (myid - rank_max)*(s_num_data - 1) ;
        dmax = dmin + (s_num_data - 1) ;
    }

    /* scale factor that accounts for the variation in crystal size */
    phi = malloc(s_num_data * sizeof(double)) ;
    updated_phi = malloc(s_num_data * sizeof(double)) ;
    fp = fopen(start_phi_file, "r") ;
    if (!fp){
        if (myid == 0)
            printf("cannot open %s!!\n", start_phi_file) ;
        return 0 ;
    }
    else{
        if (myid == 0)
            printf("initialize phi from: %s\n", start_phi_file) ;
        f_id = 0 ;
        for (d = 0 ; d < num_data ; d++){
            fscanf(fp, "%lf", &val) ;
            if (d < dmin || d > dmax-1)
                continue ;
            phi[f_id] = val ;
            f_id += 1 ;
        }
        fclose(fp) ;
    }

    prob_orien = malloc(s_num_data * sizeof(int *)) ;
    num_prob_orien = calloc(s_num_data, sizeof(int)) ;
    fp = fopen(prob_orien_file, "rb") ;
    f_id = 0 ;
    for (d = 0 ; d < num_data ; d++){
        fread(&num, sizeof(int), 1, fp) ;
        if (d < dmin || d > dmax-1){
            fseek(fp, num*sizeof(int), SEEK_CUR) ;
            continue ;
        }
        num_prob_orien[f_id] = num ;
        prob_orien[f_id] = malloc(num * sizeof(int)) ;
        fread(prob_orien[f_id], sizeof(int), num, fp) ;
        f_id += 1 ;
    }
    fclose(fp) ;

    MPI_File_open(MPI_COMM_WORLD, bg_file, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh) ;
    ave_bg = malloc(s_num_data * sizeof(double *)) ;
    log_bg = malloc(s_num_data * sizeof(double *)) ;
    
    f_id = 0 ;
    for (d = dmin ; d < dmax ; d++){
        offset = d ;
        offset *= qlen*sizeof(double) ;
        ave_bg[f_id] = malloc(qlen * sizeof(double)) ;
        log_bg[f_id] = malloc(qlen * sizeof(double)) ;
        MPI_File_read_at(fh, offset, ave_bg[f_id], qlen, MPI_DOUBLE, &status) ;
        for (t = 0 ; t < qlen ; t++){
            val = ave_bg[f_id][t] ;
            if (val < epsilon)
                log_bg[f_id][t] = 0. ;
            else
                log_bg[f_id][t] = log(val) ;
        }
        f_id += 1 ;
    }
    MPI_File_close(&fh) ;

    if (s_num_data % iter_data_block == 0)
        data_block_sz = s_num_data / iter_data_block ;
    else
        data_block_sz = s_num_data / iter_data_block + 1 ;

    data_frame = malloc(data_block_sz * sizeof(short *)) ;
    for (d = 0 ; d < data_block_sz ; d++)
        data_frame[d] = malloc(num_pix * sizeof(short)) ;

    if_visit_r = calloc(num_rot, sizeof(int)) ;
    r_block_id = malloc(num_rot * sizeof(int)) ;
    r_block_sz_mpi = calloc(nproc*iter_data_block, sizeof(int)) ;
    r_block_id_mpi = malloc(nproc*iter_data_block * sizeof(int *)) ;
    for (i = 0 ; i < iter_data_block ; i++){
        idx_start = i*data_block_sz ;
        for (d = idx_start ; d < idx_start + data_block_sz ; d++){
            if (d + dmin > dmax-1)
                break ;
            for (j = 0 ; j < num_prob_orien[d] ; j++)
                if_visit_r[prob_orien[d][j]] = 1 ;
        }
        
        r_block_sz = 0 ;
        for (r = 0 ; r < num_rot ; r++){
            if (if_visit_r[r] == 1){
                r_block_id[r_block_sz] = r ;
                r_block_sz += 1 ;
                if_visit_r[r] = 0 ;
            }
        }

        idx = myid*iter_data_block + i ;
        r_block_sz_mpi[idx] = r_block_sz ;
        r_block_id_mpi[idx] = malloc(r_block_sz * sizeof(int)) ;
        for (r = 0 ; r < r_block_sz ; r++)
            r_block_id_mpi[idx][r] = r_block_id[r] ;
    }
    MPI_Allreduce(MPI_IN_PLACE, r_block_sz_mpi, nproc*iter_data_block, MPI_INT, MPI_MAX, MPI_COMM_WORLD) ;
    free(r_block_id) ;

    max_r_block_sz = 0 ;
    for (i = 0 ; i < nproc ; i++){
        for (j = 0 ; j < iter_data_block ; j++){
            idx = i*iter_data_block + j ;
            if (myid != i)
                r_block_id_mpi[idx] = malloc(r_block_sz_mpi[idx] * sizeof(int)) ;
            MPI_Bcast(r_block_id_mpi[idx], r_block_sz_mpi[idx], MPI_INT, i, MPI_COMM_WORLD) ;
            if (max_r_block_sz < r_block_sz_mpi[idx])
                max_r_block_sz = r_block_sz_mpi[idx] ;
        }
    }

    if (num_rot % nproc == 0)
        s_num_rot = num_rot / nproc ;
    else
        s_num_rot = num_rot / nproc + 1 ;

    r2hkl_mpi = malloc(num_rot * sizeof(int *)) ;
    num_r2hkl = malloc(num_rot * sizeof(int)) ;
    if (myid == 0){
        fp = fopen(r2peak_file, "rb") ;
        for (r = 0 ; r < num_rot ; r++){
            fread(&num_r2hkl[r], sizeof(int), 1, fp) ;
            byte_offset = num_r2hkl[r]*3*sizeof(int) ;
            fseek(fp, byte_offset, SEEK_CUR) ;
        }
        fclose(fp) ;
    }
    MPI_Bcast(num_r2hkl, num_rot, MPI_INT, 0, MPI_COMM_WORLD) ;

    max_num_r2hkl = 0 ;
    for (r = 0 ; r < num_rot ; r++){
        if (num_r2hkl[r] > max_num_r2hkl)
            max_num_r2hkl = num_r2hkl[r] ;
    }

    byte_offset = 0 ;
    for (j = 0 ; j < myid ; j++){
        idx_start = j*s_num_rot ;
        for (r = 0 ; r < s_num_rot ; r++){
            idx = r + idx_start ;
            byte_offset += 1 + 3*num_r2hkl[idx] ;
        }
    }
    byte_offset *= sizeof(int) ;

    r2hkl = malloc(s_num_rot * sizeof(int *)) ;
    MPI_File_open(MPI_COMM_WORLD, r2peak_file, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh) ;
    offset = byte_offset ;
    idx_start = myid*s_num_rot ;
    for (r = 0 ; r < s_num_rot ; r++){
        idx = idx_start + r ;
        if (idx == num_rot)
            break ;
        MPI_File_read_at(fh, offset, &length, 1, MPI_INT, &status) ;
        offset += sizeof(int) ;
        length *= 3 ;
        r2hkl[r] = malloc(length * sizeof(int)) ;
        MPI_File_read_at(fh, offset, r2hkl[r], length, MPI_INT, &status) ;
        offset += length*sizeof(int) ;
    }
    MPI_File_close(&fh) ;

    set_intens() ;

    num_hkl2r = malloc(num_hkl * sizeof(int)) ;
    if (myid == 0){
        fp = fopen(peak2r_file, "rb") ;
        for (i = 0 ; i < num_hkl ; i++){
            fread(&num_hkl2r[i], sizeof(int), 1, fp) ;
            byte_offset = num_hkl2r[i]*3*sizeof(int) ;
            fseek(fp, byte_offset, SEEK_CUR) ;
        }
        fclose(fp) ;
    }
    MPI_Bcast(num_hkl2r, num_hkl, MPI_INT, 0, MPI_COMM_WORLD) ;

    max_num_hkl2r = 0 ;
    for (i = 0 ; i < num_hkl ; i++){
        if (num_hkl2r[i] > max_num_hkl2r)
            max_num_hkl2r = num_hkl2r[i] ;
    }
    
    hkl2r_mpi = malloc(3*max_num_hkl2r * sizeof(int)) ;
    relevant_pix_id = calloc(num_pix, sizeof(int)) ;
    tomo_val = malloc(num_pix * sizeof(float)) ;
    tomo_idx = malloc(num_pix * sizeof(int)) ;
    logP_offset = malloc(s_num_data * sizeof(double)) ;

    len_d2r = malloc(s_num_data * sizeof(int)) ;
    count_d2r = malloc(s_num_data * sizeof(int)) ;
    pval_d2r = malloc(s_num_data * sizeof(double *)) ;
    idx_d2r = malloc(s_num_data * sizeof(int *)) ;
    for (d = 0 ; d < s_num_data ; d++){
        len_d2r[d] = num_prob_orien[d] ;
        pval_d2r[d] = malloc(len_d2r[d] * sizeof(double)) ;
        idx_d2r[d] = malloc(len_d2r[d] * sizeof(int)) ;
    }

    Emat = malloc(num_rot * sizeof(E_matrix)) ;
    for (r = 0 ; r < num_rot ; r++)
        Emat[r].npix = 0 ;

    return 1 ;
}


void set_intens(){

    FILE *fp ;
    char infile[256] ;
    double basis_vec[3][3], cofactor[3][3], mat[2][2] ;
    double a[3], b[3], c[3], outer_prod[3], rvec[3] ;
    int i0, j0, i1, j1, if_enclosed, rand_start = 0 ;
    int i, j, k, t, x, y, z, tx, ty, tz, idx, vx, vy, vz ;
    int Nx_half, Ny_half, Nz_half, hcur, kcur, lcur ;
    int hkl_id, max_hkl_len, max_num_vox, gw_ceil ;
    int *is_partial, *num_vox, indices[3] ;
    double dr2, dx, dy, dz, qx, qy, qz, rand_val ;
    double rescale, norm, gw2, volume, det, peak_radi2 ;
    long long vox_id, **tmp_hkl2vox_id ;

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
            min_rcell_sz = norm ;
        else{
            if (min_rcell_sz > norm)
                min_rcell_sz = norm ;
        }
    }
    min_rcell_sz = sqrt(min_rcell_sz) ;
    peak_radi2 = pow(min_rcell_sz*gw/VN, 2) ;

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

    if (myid == 0){
        printf("Nx, Ny, Nz = %d, %d, %d\n", Nx, Ny, Nz) ;
        printf("num_hkl = %d\n", num_hkl) ;
    }

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
                if (myid == 0)
                    printf("error in constructing box_vox_id[%d]!!\n", i) ;
                break ;
            }
            box_vox_id[i][x*box_len2 + y*box_len + z] = j ;
        }
    }

    intens1 = malloc(num_hkl * sizeof(double *)) ;
    for (i = 0 ; i < num_hkl ; i++)
        intens1[i] = malloc(num_hkl2vox_id[i] * sizeof(double)) ;

    fp = fopen(start_intens_file, "rb") ;
    if (!fp){
        rand_start = 1 ;
        if (myid == 0){
            printf("initialize intens1 from random numbers\n") ;
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

            sym_intens() ;
        }
        for (i = 0 ; i < num_hkl ; i++)
            MPI_Bcast(intens1[i], num_hkl2vox_id[i], MPI_DOUBLE, 0, MPI_COMM_WORLD) ;
    }
    else{
        if (myid == 0)
            printf("initialize intens1 from %s\n", start_intens_file) ;
        for (i = 0 ; i < num_hkl ; i++)
            fread(intens1[i], sizeof(double), num_hkl2vox_id[i], fp) ;
        fclose(fp) ;
    }

    MPI_Barrier(MPI_COMM_WORLD) ;
    if (rand_start == 1){
        if (myid == 0){
            fp = fopen(start_intens_file, "wb") ;
            for (i = 0 ; i < num_hkl ; i++)
                fwrite(intens1[i], sizeof(double), num_hkl2vox_id[i], fp) ;
            fclose(fp) ;
        }
    }

    rescale = VN/min_rcell_sz/wl_D ;
    for (t = 0 ; t < num_pix ; t++){
        for (i = 0 ; i < 3 ; i++)
            pix[t][i] *= rescale ;
    }

    for (i = 0 ; i < max_hkl_len ; i++)
        free(tmp_hkl2vox_id[i]) ;
    free(tmp_hkl2vox_id) ;
    free(is_partial) ;
    free(num_vox) ;
}


void free_mem(){

    int i, d, r, f_id, dmin, dmax, rank_max, idx_start ;

    free(quat) ;
    free(rid_map) ;
    free(pix) ;
    free(polar) ;
    free(inv_polar) ;
    free(log_polar) ;
    free(qid_map) ;
    free(pix_map) ;
    free(phi) ;
    free(updated_phi) ;

    rank_max = num_data - nproc*(s_num_data - 1) ;
    if (myid < rank_max){
        dmin = myid*s_num_data ;
        dmax = dmin + s_num_data ;
    }
    else{
        dmin = rank_max*s_num_data + (myid - rank_max)*(s_num_data - 1) ;
        dmax = dmin + (s_num_data - 1) ;
    }

    f_id = 0 ;
    for (d = dmin ; d < dmax ; d++){
        free(prob_orien[f_id]) ;
        free(ave_bg[f_id]) ;
        free(log_bg[f_id]) ;
        f_id += 1 ;
    }
    free(prob_orien) ;
    free(num_prob_orien) ;
    free(ave_bg) ;
    free(log_bg) ;

    for (d = 0 ; d < data_block_sz ; d++)
        free(data_frame[d]) ;
    free(data_frame) ;

    for (i = 0 ; i < nproc*iter_data_block ; i++)
        free(r_block_id_mpi[i]) ;
    free(r_block_id_mpi) ;
    free(r_block_sz_mpi) ;
    free(if_visit_r) ;

    idx_start = myid*s_num_rot ;
    for (r = 0 ; r < s_num_rot ; r++){
        if (r + idx_start == num_rot)
            break ;
        free(r2hkl[r]) ;
    }
    free(r2hkl) ;
    free(num_r2hkl) ;
    free(r2hkl_mpi) ;

    for (i = 0 ; i < num_hkl ; i++){
        free(hkl2vox_id[i]) ;
        free(box_vox_id[i]) ;
        free(intens1[i]) ;
    }
    free(hkl2vox_id) ;
    free(box_vox_id) ;
    free(intens1) ;
    free(num_hkl2vox_id) ;
    free(hkl_table) ;

    free(num_hkl2r) ;
    free(hkl2r_mpi) ;
    free(relevant_pix_id) ;
    free(tomo_val) ;
    free(tomo_idx) ;
    free(logP_offset) ;

    for (d = 0 ; d < s_num_data ; d++){
        free(pval_d2r[d]) ;
        free(idx_d2r[d]) ;
    }
    free(pval_d2r) ;
    free(idx_d2r) ;
    free(len_d2r) ;
    free(count_d2r) ;
    free(Emat) ;
}


void initialize(){

    int r ;
    double inv_w, total_weight = 0. ;

    /* normalize and take logarithm of the weight of quat */
    for (r = 0 ; r < num_rot ; r++)
        total_weight += quat[r][4] ;

    inv_w = 1./total_weight ;
    for (r = 0 ; r < num_rot ; r++)
        quat[r][4] *= inv_w ;

    for (r = 0 ; r < num_rot ; r++)
        quat[r][4] = log(quat[r][4]) ;
}


void print_intens(){

    FILE *fp ;
    char outfile[256] ;
    int i, j ;
    double val ;

    sprintf(outfile, "%s/iter_flag-%d/finish_intensity.bin", prob_dir, iter_flag) ;
    fp = fopen(outfile, "wb") ;
    for (i = 0 ; i < num_hkl ; i++)
        fwrite(intens1[i], sizeof(double), num_hkl2vox_id[i], fp) ;
    fclose(fp) ;

    sprintf(outfile, "%s/iter_flag-%d/hkl-val.dat", prob_dir, iter_flag) ;
    fp = fopen(outfile, "w") ;
    fprintf(fp, "%d %d %d %d\n", hmax, kmax, lmax, num_hkl) ;
    for (i = 0 ; i < num_hkl ; i++){
        val = 0. ;
        for (j = 0 ; j < num_hkl2vox_id[i] ; j++)
            val += intens1[i][j] ;
        fprintf(fp, "%d %d %d %1.7e\n", hkl_table[3*i], hkl_table[3*i+1], hkl_table[3*i+2], val) ;
    }
    fclose(fp) ;
}
