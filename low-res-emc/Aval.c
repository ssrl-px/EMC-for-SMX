#include "emc.h"

void calAval( int data_block_id ){

    /* Aval is the cross correlation between Pjk and Kik */
    int i, j, k, d, m, n, r, s, t, rid, pid, qid, peak_ct, dmin, dmax, rank_max ;
    int x, y, z, x0, y0, z0, tx, ty, tz, ct, len, tomo_len, pix_ct, vox_ct ;
    int idx, idx_offset, idx_start, count_d2r_d, hkl_id, photon_count ;
    int my_Acount_d, my_pval_ct_d, my_Alen_d, my_pval_len_d, num_r2hkl_rid ;
    int r_block_sz, row_min, col_min, row_max, col_max, search_ct, max_iter = 200 ;
    int if_overshoot, *idx_d2r_d, *r2hkl_mpi_rid, *box_vox_id_hkl, indices[3] ;
    int *hkl_npix_mpi, *my_Aidx_d, *my_pval_id_d, *my_pix_ct_d ;
    double fx, fy, fz, cx, cy, cz, vx, vy, vz, dx, dy, dz, dr2, wval ;
    double cur_phi, grad, w, Bval, dphi, rescale, lambda, inv_tau ;
    double S_a, S_b, S_c, S_d, phi_a, phi_b, phi_c, phi_d, bv, pval ;
    double *my_wval_d, *pval_d2r_d, *ave_bg_d, *log_bg_d, *intens1_hkl, rvec[3] ;
    double gw2 = gw*gw, rot_pix[3], epsilon = 1.e-10, tolerance = 1.e-6 ;
    float *my_Aval_d, *my_Aval_hkl, *my_Bval_hkl, *pix_pid ;
    short *data_frame_d ;
    long long tmplen ;

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

    hkl_npix_mpi = malloc(max_num_r2hkl*2 * sizeof(int)) ;
    if (iter_flag % 2 == 0){
        my_Alen_d = mem_step_sz ;
        my_Aval_d = malloc(mem_step_sz * sizeof(float)) ;
        my_Aidx_d = malloc(mem_step_sz * sizeof(int)) ;
        my_wval_d = malloc(mem_step_sz * sizeof(double)) ;
        my_pval_len_d = mem_step_sz ;
        my_pval_id_d = malloc(mem_step_sz * sizeof(int)) ;
        my_pix_ct_d = malloc(mem_step_sz * sizeof(int)) ;
    }

    inv_tau = 1./golden_ratio ;
    rescale = VN/min_rcell_sz ;
    idx_start = data_block_id*data_block_sz ;
    for (d = idx_start ; d < idx_start + data_block_sz ; d++){
        if (d + dmin > dmax-1)
            break ;
        count_d2r_d = count_d2r[d] ;
        pval_d2r_d = pval_d2r[d] ;
        idx_d2r_d = idx_d2r[d] ;
        data_frame_d = data_frame[d - idx_start] ;
        ave_bg_d = ave_bg[d] ;
        log_bg_d = log_bg[d] ;
        phi_d = phi[d] ;

        /* frames with phi = 0. are irrelevant when updating vox_val */
        if (phi_d < epsilon && iter_flag % 2 == 1)
            continue ;

        if (iter_flag % 2 == 0){
            my_Acount_d = 0 ;
            my_pval_ct_d = 0 ;
        }

        for (s = 0 ; s < count_d2r_d ; s++){
            rid = idx_d2r_d[s] ;
            r = rid_map[rid] ;
            pval = pval_d2r_d[s] ;
            if (rid != r_block_id[r])
                printf("error: rid != r_block_id[r]!!\n") ;
            r2hkl_mpi_rid = r2hkl_mpi[r] ;
            num_r2hkl_rid = num_r2hkl[rid] ;
            make_rot(rid) ;

            /* compute tomogram */
            tomo_len = 0 ;
            peak_ct = 0 ;
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

                len = tomo_len ;
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
                len = tomo_len - len ;
                hkl_npix_mpi[2*peak_ct] = hkl_id ;
                hkl_npix_mpi[2*peak_ct+1] = len ;
                peak_ct += 1 ;
            }

            if (iter_flag % 2 == 0){
                ct = my_Acount_d ;
                len = my_Alen_d ;
                pix_ct = 0 ;
                for (t = 0 ; t < tomo_len ; t++){
                    /* pixels with fabs(wij) = 0. are irrelevant when updating phi */
                    if (fabs(tomo_val[t]) < epsilon)
                        continue ;
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
                    /* increase array size */
                    if (ct == len){
                        tmplen = len + mem_step_sz ;
                        if (tmplen > int_max){
                            printf("myid = %d, d = %d: ", myid, d) ;
                            printf("my_Alen overflows int32!!\n") ;
                        }
                        len += mem_step_sz ;
                        my_Alen_d = len ;
                        my_Aval_d = realloc(my_Aval_d, len*sizeof(float)) ;
                        my_Aidx_d = realloc(my_Aidx_d, len*sizeof(int)) ;
                        my_wval_d = realloc(my_wval_d, len*sizeof(double)) ;
                    }
                    my_Aval_d[ct] = pval*photon_count ;
                    my_Aidx_d[ct] = pid ;
                    my_wval_d[ct] = tomo_val[t] ;
                    ct += 1 ;
                    pix_ct += 1 ;
                }
                my_Acount_d = ct ;
                if (my_pval_ct_d == my_pval_len_d){
                    tmplen = my_pval_len_d + mem_step_sz ;
                    if (tmplen > int_max){
                        printf("myid = %d, d = %d: ", myid, d) ;
                        printf("my_pval_len_d overflows int32!!\n") ;
                    }
                    my_pval_len_d += mem_step_sz ;
                    my_pix_ct_d = realloc(my_pix_ct_d, my_pval_len_d*sizeof(int)) ;
                    my_pval_id_d = realloc(my_pval_id_d, my_pval_len_d*sizeof(int)) ;
                }
                my_pix_ct_d[my_pval_ct_d] = pix_ct ;
                my_pval_id_d[my_pval_ct_d] = s ;
                my_pval_ct_d += 1 ;
            }
            else{
                idx_offset = 0 ;
                for (k = 0 ; k < peak_ct ; k++){
                    idx = 2*k ;
                    hkl_id = hkl_npix_mpi[idx] ;
                    pix_ct = hkl_npix_mpi[idx+1] ;
                    ct = my_Acount[hkl_id] ;
                    len = my_Alen[hkl_id] ;
                    my_Aval_hkl = my_Aval[hkl_id] ;
                    my_Bval_hkl = my_Bval[hkl_id] ;
                    for (t = 0 ; t < pix_ct ; t++){
                        pid = tomo_idx[idx_offset + t] ;
                        qid = qid_map[pid] ;
                        photon_count = data_frame_d[pid] ;
                        /* increase array size */
                        if (ct == len){
                            tmplen = len + mem_step_sz ;
                            if (tmplen > int_max){
                                printf("myid = %d, hkl_id = %d: ", myid, hkl_id) ;
                                printf("my_Alen overflows int32!!\n") ;
                            }
                            len += mem_step_sz ;
                            my_Alen[hkl_id] = len ;
                            my_Aval[hkl_id] = realloc(my_Aval[hkl_id], len*sizeof(float)) ;
                            my_Bval[hkl_id] = realloc(my_Bval[hkl_id], len*sizeof(float)) ;
                            my_Aval_hkl = my_Aval[hkl_id] ;
                            my_Bval_hkl = my_Bval[hkl_id] ;
                        }
                        if (photon_count < 0)
                            my_Aval_hkl[ct] = -1 ;
                        else
                            my_Aval_hkl[ct] = pval*photon_count ;

                        if (qid < 0)
                            my_Bval_hkl[ct] = -1 ;
                        else{
                            bv = ave_bg_d[qid] ;
                            if (bv < epsilon)
                                my_Bval_hkl[ct] = -1 ;
                            else
                                my_Bval_hkl[ct] = bv/phi_d ;
                        }
                        ct += 1 ;
                    }
                    my_Acount[hkl_id] = ct ;
                    idx_offset += pix_ct ;
                    ct = my_pval_ct[hkl_id] ;
                    len = my_pval_len[hkl_id] ;
                    /* increase array size */
                    if (ct == len){
                        tmplen = len + mem_step_sz ;
                        if (tmplen > int_max){
                            printf("myid = %d, hkl_id = %d: ", myid, hkl_id) ;
                            printf("my_pval_len overflows int32!!\n") ;
                        }
                        len += mem_step_sz ;
                        my_pval_len[hkl_id] = len ;
                        my_pval_id[hkl_id] = realloc(my_pval_id[hkl_id], len*sizeof(int)) ;
                    }
                    my_pval_id[hkl_id][ct] = d ;
                    my_pval_id[hkl_id][ct+1] = s ;
                    my_pval_ct[hkl_id] += 2 ;
                }
            }
        }

        /* maximize phi */
        if (iter_flag % 2 == 0){
            if (count_d2r_d == 0){
                updated_phi[d] = 0. ;
                continue ;
            }

            S_a = 0. ;
            grad = 0. ;
            phi_a = 0. ;
            idx_offset = 0 ;
            for (j = 0 ; j < my_pval_ct_d ; j++){
                i = my_pval_id_d[j] ;
                pval = pval_d2r_d[i] ;
                pix_ct = my_pix_ct_d[j] ;
                for (t = 0 ; t < pix_ct ; t++){
                    idx = idx_offset + t ;
                    pid = my_Aidx_d[idx] ;
                    qid = qid_map[pid] ;
                    w = my_wval_d[idx] ;
                    bv = ave_bg_d[qid] ;
                    grad += (pval*polar[pid] - my_Aval_d[idx]/bv)*w ;
                    S_a += pval*(bv*polar[pid]) - my_Aval_d[idx]*(log_bg_d[qid] + log_polar[pid]) ;
                }
                idx_offset += pix_ct ;
            }

            if (grad >= -epsilon){
                updated_phi[d] = 0. ;
                continue ;
            }

            /* bracket the minimum */
            S_b = S_a ;
            lambda = 1. ;
            if_overshoot = 0 ;
            while (S_b <= S_a && if_overshoot == 0){
                phi_b = phi_a + lambda ;
                S_b = 0. ;
                idx_offset = 0 ;
                for (j = 0 ; j < my_pval_ct_d ; j++){
                    i = my_pval_id_d[j] ;
                    pval = pval_d2r_d[i] ;
                    pix_ct = my_pix_ct_d[j] ;
                    if (if_overshoot == 1)
                        break ;
                    for (t = 0 ; t < pix_ct ; t++){
                        idx = idx_offset + t ;
                        pid = my_Aidx_d[idx] ;
                        qid = qid_map[pid] ;
                        Bval = (ave_bg_d[qid] + phi_b*my_wval_d[idx])*polar[pid] ;
                        if (Bval < epsilon){
                            if_overshoot = 1 ;
                            break ;
                        }
                        S_b += pval*Bval - my_Aval_d[idx]*log(Bval) ;
                    }
                    idx_offset += pix_ct ;
                }
                if (S_b > S_a || if_overshoot == 1)
                    break ;
                lambda *= 2 ;
            }

            phi_b = phi_a + lambda ;
            if (if_overshoot == 1){
                idx_offset = 0 ;
                for (j = 0 ; j < my_pval_ct_d ; j++){
                    i = my_pval_id_d[j] ;
                    pval = pval_d2r_d[i] ;
                    pix_ct = my_pix_ct_d[j] ;
                    for (t = 0 ; t < pix_ct ; t++){
                        idx = idx_offset + t ;
                        pid = my_Aidx_d[idx] ;
                        qid = qid_map[pid] ;
                        w = my_wval_d[idx] ;
                        if (w < 0.){
                            cur_phi = -ave_bg_d[qid]/w ;
                            if (phi_b > cur_phi)
                                phi_b = cur_phi ;
                        }
                    }
                    idx_offset += pix_ct ;
                }
            }

            /* golden section search */
            search_ct = 0 ;
            dphi = phi_b - phi_a ;
            if (dphi < 0.)
                printf("error: myid = %d, d = %d, phi_b < phi_a!!\n", myid, d) ;

            cur_phi = (phi_b + phi_a)/2. ;
            while (dphi > cur_phi*tolerance && search_ct < max_iter){
                dphi = phi_b - phi_a ;
                phi_c = phi_b - dphi*inv_tau ;
                phi_d = phi_a + dphi*inv_tau ;
                S_c = 0. ;
                S_d = 0. ;
                idx_offset = 0 ;
                for (j = 0 ; j < my_pval_ct_d ; j++){
                    i = my_pval_id_d[j] ;
                    pval = pval_d2r_d[i] ;
                    pix_ct = my_pix_ct_d[j] ;
                    for (t = 0 ; t < pix_ct ; t++){
                        idx = idx_offset + t ;
                        pid = my_Aidx_d[idx] ;
                        qid = qid_map[pid] ;
                        w = my_wval_d[idx] ;
                        Bval = (ave_bg_d[qid] + phi_c*w)*polar[pid] ;
                        S_c += pval*Bval - my_Aval_d[idx]*log(Bval) ;
                        Bval = (ave_bg_d[qid] + phi_d*w)*polar[pid] ;
                        S_d += pval*Bval - my_Aval_d[idx]*log(Bval) ;
                    }
                    idx_offset += pix_ct ;
                }

                if (S_c < S_d)
                    phi_b = phi_d ;
                else
                    phi_a = phi_c ;

                dphi = phi_b - phi_a ;
                cur_phi = (phi_b + phi_a)/2. ;
                search_ct += 1 ;
            }

            updated_phi[d] = cur_phi ;
            if (updated_phi[d] < 0.)
                printf("error: myid = %d, d = %d, updated_phi < 0!!\n", myid, d) ;
        }
    }

    free(hkl_npix_mpi) ;
    if (iter_flag % 2 == 0){
        free(my_Aval_d) ;
        free(my_Aidx_d) ;
        free(my_wval_d) ;
        free(my_pval_id_d) ;
        free(my_pix_ct_d) ;
    }
}
