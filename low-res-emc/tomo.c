#include "emc.h"

void calprob( int data_block_id ){

    FILE *fp ;
    char out_prob_dir[256], out_prob_file[256] ;
    int i, j, k, d, r, t, m, n, x, y, z, tx, ty, tz, pid, qid, rid ;
    int idx, idx_start, idx_offset, photon_count, num_r2hkl_rid, tmplen ;
    int row_min, row_max, col_min, col_max, pix_ct, vox_ct, hkl_id ;
    int x0, y0, z0, tomo_len, r_block_sz, dmin, dmax, rank_max ;
    int high_p_count, *prob_orien_d, *idx_d2r_d, *high_p_idx, indices[3] ;
    int *r2hkl_mpi_rid, **r2d_table, *r2d_table_ct, *box_vox_id_hkl ;
    double fx, fy, fz, cx, cy, cz, vx, vy, vz, wval, dx, dy, dz, dr2, pval ;
    double logP_offset_d, rescale, high_p_sum, p_sum, inv_p_sum, w, bv ;
    double max_p, phi_d, gw2 = gw*gw, rot_pix[3], rvec[3], epsilon = 1.e-10 ;
    double *intens1_hkl, *high_p_val, *pval_d2r_d, *ave_bg_d, *log_bg_d ;
    short *data_frame_d ;
    float *pix_pid ;

    rank_max = num_data - nproc*(s_num_data - 1) ;
    if (myid < rank_max){
        dmin = myid*s_num_data ;
        dmax = dmin + s_num_data ;
    }
    else{
        dmin = rank_max*s_num_data + (myid - rank_max)*(s_num_data - 1) ;
        dmax = dmin + (s_num_data - 1) ;
    }

    idx = myid*iter_data_block + data_block_id ;
    r_block_id = r_block_id_mpi[idx] ;
    r_block_sz = r_block_sz_mpi[idx] ;
    for (r = 0 ; r < r_block_sz ; r++){
        rid = r_block_id[r] ;
        rid_map[rid] = r ;
    }

    r2d_table = malloc(r_block_sz * sizeof(int *)) ;
    r2d_table_ct = calloc(r_block_sz, sizeof(int)) ;
    idx_start = data_block_id*data_block_sz ;
    for (d = idx_start ; d < idx_start + data_block_sz ; d++){
        if (d + dmin > dmax-1)
            break ;
        prob_orien_d = prob_orien[d] ;
        for (i = 0 ; i < num_prob_orien[d] ; i++){
            r = rid_map[prob_orien_d[i]] ;
            r2d_table_ct[r] += 1 ;
        }
    }

    for (r = 0 ; r < r_block_sz ; r++)
        r2d_table[r] = malloc(r2d_table_ct[r] * sizeof(int)) ;

    memset(r2d_table_ct, 0, r_block_sz*sizeof(int)) ;
    for (d = idx_start ; d < idx_start + data_block_sz ; d++){
        if (d + dmin > dmax-1)
            break ;
        prob_orien_d = prob_orien[d] ;
        for (i = 0 ; i < num_prob_orien[d] ; i++){
            r = rid_map[prob_orien_d[i]] ;
            r2d_table[r][r2d_table_ct[r]] = d ;
            r2d_table_ct[r] += 1 ;
        }
    }

    rescale = VN/min_rcell_sz ;
    for (r = 0 ; r < r_block_sz ; r++){
        rid = r_block_id[r] ;
        r2hkl_mpi_rid = r2hkl_mpi[r] ;
        num_r2hkl_rid = num_r2hkl[rid] ;
        make_rot(rid) ;
        /* compute tomogram */
        tomo_len = 0 ;
        for (n = 0 ; n < num_r2hkl_rid ; n++){
            idx = 3*n ;
            hkl_id = r2hkl_mpi_rid[idx] ;
            pid = r2hkl_mpi_rid[idx+1] ;
            row_min = pid / num_col ;
            col_min = pid % num_col ;
            pid = r2hkl_mpi_rid[idx+2] ;
            row_max = pid / num_col + 1 ;
            col_max = pid % num_col + 1 ;
            pix_ct = 0 ;
            for (i = row_min ; i < row_max ; i++){
                idx_offset = i*num_col ;
                for (j = col_min ; j < col_max ; j++){
                    pid = pix_map[idx_offset + j] ;
                    if (pid < 0)
                        continue ;
                    relevant_pix_id[pix_ct] = pid ;
                    pix_ct += 1 ;
                }
            }

            box_vox_id_hkl = box_vox_id[hkl_id] ;
            intens1_hkl = intens1[hkl_id] ;
            idx = 3*hkl_id ;
            for (i = 0 ; i < 3 ; i++)
                indices[i] = hkl_table[idx+i] ;

            for (i = 0 ; i < 3 ; i++){
                rvec[i] = 0. ;
                for (j = 0 ; j < 3 ; j++)
                    rvec[i] += rlatt_vec[i][j]*indices[j] ;
            }

            vx = rvec[0]*rescale ;
            vy = rvec[1]*rescale ;
            vz = rvec[2]*rescale ;

            tx = ((int) round(vx)) - box_half_len ;
            ty = ((int) round(vy)) - box_half_len ;
            tz = ((int) round(vz)) - box_half_len ;

            for (t = 0 ; t < pix_ct ; t++){
                pid = relevant_pix_id[t] ;
                pix_pid = pix[pid] ;
                for (i = 0 ; i < 3 ; i++){
                    rot_pix[i] = rot[i][0]*pix_pid[0] ;
                    for (j = 1 ; j < 3 ; j++)
                        rot_pix[i] += rot[i][j]*pix_pid[j] ;
                }

                dx = rot_pix[0] - vx ;
                dy = rot_pix[1] - vy ;
                dz = rot_pix[2] - vz ;
                dr2 = dx*dx + dy*dy + dz*dz ;
                if (dr2 > gw2)
                    continue ;

                x = rot_pix[0] - tx ;
                y = rot_pix[1] - ty ;
                z = rot_pix[2] - tz ;

                fx = rot_pix[0] - floor(rot_pix[0]) ;
                fy = rot_pix[1] - floor(rot_pix[1]) ;
                fz = rot_pix[2] - floor(rot_pix[2]) ;

                wval = 0. ;
                vox_ct = 0 ;
                for (i = 0 ; i < 2 ; i++){
                for (j = 0 ; j < 2 ; j++){
                for (k = 0 ; k < 2 ; k++){
                    x0 = x+i ;
                    if (abs(x0 - box_half_len) > box_half_len)
                        continue ;
                    y0 = y+j ;
                    if (abs(y0 - box_half_len) > box_half_len)
                        continue ;
                    z0 = z+k ;
                    if (abs(z0 - box_half_len) > box_half_len)
                        continue ;
                    idx = x0*box_len2 + y0*box_len + z0 ;
                    m = box_vox_id_hkl[idx] ;
                    if (m < 0)
                        continue ;
                    if (i == 0)
                        cx = 1. - fx ;
                    else
                        cx = fx ;
                    if (j == 0)
                        cy = 1. - fy ;
                    else
                        cy = fy ;
                    if (k == 0)
                        cz = 1. - fz ;
                    else
                        cz = fz ;
                    wval += cx*cy*cz*intens1_hkl[m] ;
                    vox_ct += 1 ;
                }
                }
                }

                if (vox_ct == 0)
                    continue ;

                tomo_val[tomo_len] = wval ;
                tomo_idx[tomo_len] = pid ;
                tomo_len += 1 ;
            }
        }

        /* calculate probabilities */
        for (i = 0 ; i < r2d_table_ct[r] ; i++){
            d = r2d_table[r][i] ;
            logP_offset_d = logP_offset[d] ;
            phi_d = phi[d] ;
            pval_d2r_d = pval_d2r[d] ;
            idx_d2r_d = idx_d2r[d] ;
            pval = quat[rid][4] + logP_offset_d ;

            if (epsilon < phi_d){
                ave_bg_d = ave_bg[d] ;
                log_bg_d = log_bg[d] ;
                data_frame_d = data_frame[d - idx_start] ;
                if (iter_flag % 2 == 0){
                    for (t = 0 ; t < tomo_len ; t++){
                        pid = tomo_idx[t] ;
                        if (tomo_val[t] < epsilon)
                            continue ;
                        qid = qid_map[pid] ;
                        if (qid < 0)
                            continue ;
                        bv = ave_bg_d[qid] ;
                        if (bv < epsilon)
                            continue ;
                        photon_count = data_frame_d[pid] ;
                        /* ignore the dead pixels */
                        if (photon_count < 0)
                            continue ;
                        w = tomo_val[t]*phi_d*polar[pid] ;
                        /* the contribution from bv to pval
                           has been included in logP_offset_d */
                        pval -= w ;
                        w += bv*polar[pid] ;
                        if (photon_count > 0)
                            pval += photon_count*(log(w) - log_bg_d[qid] - log_polar[pid]) ;
                    }
                }
                else{
                    for (t = 0 ; t < tomo_len ; t++){
                        pid = tomo_idx[t] ;
                        qid = qid_map[pid] ;
                        if (qid < 0)
                            continue ;
                        bv = ave_bg_d[qid] ;
                        if (bv < epsilon)
                            continue ;
                        photon_count = data_frame_d[pid] ;
                        /* ignore the dead pixels */
                        if (photon_count < 0)
                            continue ;
                        w = tomo_val[t]*phi_d*polar[pid] ;
                        /* the contribution from bv to pval
                           has been included in logP_offset_d */
                        pval -= w ;
                        w += bv*polar[pid] ;
                        if (w < epsilon){
                            w = tomo_val[t]*phi_d*polar[pid] ;
                            pval += w ;
                            continue ;
                        }
                        if (photon_count > 0)
                            pval += photon_count*(log(w) - log_bg_d[qid] - log_polar[pid]) ;
                    }
                }
            }

            pval_d2r_d[count_d2r[d]] = pval ;
            idx_d2r_d[count_d2r[d]] = rid ;
            count_d2r[d] += 1 ;
        }
    }

    /* normalize probabilities */
    tmplen = 0 ;
    for (d = idx_start ; d < idx_start + data_block_sz ; d++){
        if (d + dmin > dmax-1)
            break ;
        if (tmplen < count_d2r[d])
            tmplen = count_d2r[d] ;
    }
    high_p_val = malloc(tmplen * sizeof(double)) ;
    high_p_idx = malloc(tmplen * sizeof(int)) ;

    for (d = idx_start ; d < idx_start + data_block_sz ; d++){
        if (d + dmin > dmax-1)
            break ;
        if (count_d2r[d] == 0)
            continue ;

        pval_d2r_d = pval_d2r[d] ;
        idx_d2r_d = idx_d2r[d] ;
        max_p = -1.e30 ;
        for (i = 0 ; i < count_d2r[d] ; i++){
            if (max_p < pval_d2r_d[i])
                max_p = pval_d2r_d[i] ;
        }

        p_sum = 0. ;
        for (i = 0 ; i < count_d2r[d] ; i++){
            pval_d2r_d[i] -= max_p ;
            pval_d2r_d[i] = exp(pval_d2r_d[i]) ;
            p_sum += pval_d2r_d[i] ;
        }

        inv_p_sum = 1./p_sum ;
        for (i = 0 ; i < count_d2r[d] ; i++)
            pval_d2r_d[i] *= inv_p_sum ;

        if (count_d2r[d] == 1)
            continue ;

        /* sort w.r.t values in descending order */
        merge_sort(pval_d2r_d, high_p_val, idx_d2r_d, high_p_idx, count_d2r[d]) ;
        high_p_sum = 0. ;
        high_p_count = 0 ;
        for (i = 0 ; i < count_d2r[d] ; i++){
            if (pval_d2r_d[i] < P_MIN)
                break ;
            high_p_sum += pval_d2r_d[i] ;
            high_p_count += 1 ;
            if (high_p_sum > CP_THRES)
                break ;
        }

        if (high_p_count == 0){
            count_d2r[d] = 0 ;
            continue ;
        }

        inv_p_sum = 1./high_p_sum ;
        count_d2r[d] = high_p_count ;
        for (i = 0 ; i < count_d2r[d] ; i++){
            pval_d2r_d[i] = high_p_val[i] * inv_p_sum ;
            idx_d2r_d[i] = high_p_idx[i] ;
        }
    }

    /* output probabilities */
    sprintf(out_prob_dir, "%s/iter_flag-%d", prob_dir, iter_flag) ;
    sprintf(out_prob_file, "%s/high_p-%d.dat", out_prob_dir, myid) ;
    if (data_block_id == 0){
        fp = fopen(out_prob_file, "w") ;
        fprintf(fp, "%d\n\n", dmax-dmin) ;
    }
    else
        fp = fopen(out_prob_file, "a") ;

    for (d = idx_start ; d < idx_start + data_block_sz ; d++){
        if (d + dmin > dmax-1)
            break ;
        idx_d2r_d = idx_d2r[d] ;
        pval_d2r_d = pval_d2r[d] ;
        fprintf(fp, "%d\n", count_d2r[d]) ;
        for (i = 0 ; i < count_d2r[d] ; i++)
            fprintf(fp, "%d %.6f ", idx_d2r_d[i], pval_d2r_d[i]) ;
        fprintf(fp, "\n") ;
        /* sort w.r.t indices in ascending order */
        if (count_d2r[d] > 1)
            merge_sort_idx(pval_d2r_d, high_p_val, idx_d2r_d, high_p_idx, count_d2r[d]) ;
    }
    fclose(fp) ;

    for (r = 0 ; r < r_block_sz ; r++)
        free(r2d_table[r]) ;
    free(r2d_table) ;
    free(r2d_table_ct) ;
    free(high_p_val) ;
    free(high_p_idx) ;
}


void compress( int hkl_id ){

    int i, j, k, r, t, m, x, y, z, tx, ty, tz, rid, pid ;
    int x0, y0, z0, pix_ct, vox_ct, idx, idx_offset ;
    int row_min, col_min, row_max, col_max, total_count ;
    int cur_num_hkl2r, indices[3], *box_vox_id_hkl ;
    double dx, dy, dz, dr2, fx, fy, fz, cx, cy, cz, vx, vy, vz ;
    double t1, dt, gw2 = gw*gw, rot_pix[3], rvec[3] ;
    E_matrix Emat_r ;

    t1 = MPI_Wtime() ;
    cur_num_hkl2r = num_hkl2r[hkl_id] ;
    box_vox_id_hkl = box_vox_id[hkl_id] ;
    idx = 3*hkl_id ;
    for (i = 0 ; i < 3 ; i++)
        indices[i] = hkl_table[idx+i] ;

    for (i = 0 ; i < 3 ; i++){
        rvec[i] = 0. ;
        for (j = 0 ; j < 3 ; j++)
            rvec[i] += rlatt_vec[i][j]*indices[j] ;
    }

    vx = rvec[0]*VN/min_rcell_sz ;
    vy = rvec[1]*VN/min_rcell_sz ;
    vz = rvec[2]*VN/min_rcell_sz ;

    tx = ((int) round(vx)) - box_half_len ;
    ty = ((int) round(vy)) - box_half_len ;
    tz = ((int) round(vz)) - box_half_len ;

    total_count = 0 ;
    /* calculate Emat */
    for (r = 0 ; r < cur_num_hkl2r ; r++){
        idx = 3*r ;
        rid = hkl2r_mpi[idx] ;
        if (if_visit_r[rid] == 0)
            continue ;
        make_rot(rid) ;
        pid = hkl2r_mpi[idx+1] ;
        row_min = pid / num_col ;
        col_min = pid % num_col ;
        pid = hkl2r_mpi[idx+2] ;
        row_max = pid / num_col + 1 ;
        col_max = pid % num_col + 1 ;
        pix_ct = 0 ;
        for (i = row_min ; i < row_max ; i++){
            idx_offset = i*num_col ;
            for (j = col_min ; j < col_max ; j++){
                pid = pix_map[idx_offset + j] ;
                if (pid < 0)
                    continue ;
                relevant_pix_id[pix_ct] = pid ;
                pix_ct += 1 ;
            }
        }
        Emat[rid].npix = pix_ct ;
        Emat[rid].pix_id = malloc(pix_ct * sizeof(int)) ;
        Emat[rid].vox_ct = malloc(pix_ct * sizeof(int)) ;
        Emat[rid].vox_id = malloc(8*pix_ct * sizeof(int)) ;
        Emat[rid].weight = malloc(8*pix_ct * sizeof(double)) ;
        Emat_r = Emat[rid] ;

        pix_ct = 0 ;
        for (t = 0 ; t < Emat_r.npix ; t++){
            pid = relevant_pix_id[t] ;
            for (i = 0 ; i < 3 ; i++){
                rot_pix[i] = rot[i][0]*pix[pid][0] ;
                for (j = 1 ; j < 3 ; j++)
                    rot_pix[i] += rot[i][j]*pix[pid][j] ;
            }

            dx = rot_pix[0] - vx ;
            dy = rot_pix[1] - vy ;
            dz = rot_pix[2] - vz ;
            dr2 = dx*dx + dy*dy + dz*dz ;
            if (dr2 > gw2)
                continue ;

            x = rot_pix[0] - tx ;
            y = rot_pix[1] - ty ;
            z = rot_pix[2] - tz ;

            fx = rot_pix[0] - floor(rot_pix[0]) ;
            fy = rot_pix[1] - floor(rot_pix[1]) ;
            fz = rot_pix[2] - floor(rot_pix[2]) ;

            vox_ct = 0 ;
            idx_offset = 8*pix_ct ;
            for (i = 0 ; i < 2 ; i++){
            for (j = 0 ; j < 2 ; j++){
            for (k = 0 ; k < 2 ; k++){
                x0 = x+i ;
                if (abs(x0 - box_half_len) > box_half_len)
                    continue ;
                y0 = y+j ;
                if (abs(y0 - box_half_len) > box_half_len)
                    continue ;
                z0 = z+k ;
                if (abs(z0 - box_half_len) > box_half_len)
                    continue ;
                m = box_vox_id_hkl[x0*box_len2 + y0*box_len + z0] ;
                if (m < 0)
                    continue ;
                if (i == 0)
                    cx = 1. - fx ;
                else
                    cx = fx ;
                if (j == 0)
                    cy = 1. - fy ;
                else
                    cy = fy ;
                if (k == 0)
                    cz = 1. - fz ;
                else
                    cz = fz ;
                idx = idx_offset + vox_ct ;
                Emat_r.vox_id[idx] = m ;
                Emat_r.weight[idx] = cx*cy*cz ;
                vox_ct += 1 ;
            }
            }
            }

            if (vox_ct == 0)
                continue ;

            Emat_r.pix_id[pix_ct] = pid ;
            Emat_r.vox_ct[pix_ct] = vox_ct ;
            pix_ct += 1 ;
        }

        total_count += pix_ct ;
        Emat[rid].npix = pix_ct ;
    }

    dt = MPI_Wtime() - t1 ;
    t1 = MPI_Wtime() ;
    if (myid == 0)
        printf("hkl_id = %d, calculating Emat takes = %.3f sec, ", hkl_id, dt) ;

    maximize(cur_num_hkl2r, hkl_id) ;

    for (r = 0 ; r < cur_num_hkl2r ; r++){
        rid = hkl2r_mpi[3*r] ;
        if (if_visit_r[rid] == 0)
            continue ;
        free(Emat[rid].pix_id) ;
        free(Emat[rid].vox_ct) ;
        free(Emat[rid].vox_id) ;
        free(Emat[rid].weight) ;
    }
    dt = MPI_Wtime() - t1 ;
    if (myid == 0)
        printf("maximize takes %.3f sec, ", dt) ;
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
