#include "emc.h"

void calAval( int data_block_id ){

    /* Aval is the cross correlation between Pjk and Kik */
    int i, j, k, d, r, s, t, m, n, idx, idx_offset, idx_start ;
    int x, y, z, tx, ty, tz, pid, qid, rid, dmin, dmax, hkl_id, rank_max ;
    int x0, y0, z0, vox_ct, row_min, col_min, row_max, col_max, len, ct ;
    int tomo_len, pix_ct, peak_ct, photon_count, count_d2r_d, d2r_table_ct_d ;
    int r_block_sz, num_r2hkl_rid, *r_block_id, **d2r_table, *d2r_table_d ;
    int *hkl_npix_mpi, *idx_d2r_d, *r2hkl_mpi_rid, *box_vox_id_hkl, indices[3] ;
    double fx, fy, fz, cx, cy, cz, vx, vy, vz, dx, dy, dz, dr2, wval ;
    double bv, rescale, phi_d, pval, gw2 = gw*gw, epsilon = 1.e-10 ;
    double *pval_d2r_d, *ave_bg_d, *intens1_hkl, rot_pix[3], rvec[3] ;
    float *my_Aval_hkl, *my_Bval_hkl, *pix_pid ;
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
    d2r_table = malloc(data_block_sz * sizeof(int *)) ;
    memset(if_visit_r, 0, num_rot*sizeof(int)) ;
    idx_start = data_block_id*data_block_sz ;
    rescale = VN/min_rcell_sz ;

    for (d = idx_start ; d < idx_start + data_block_sz ; d++){
        if (d + dmin > dmax-1)
            break ;
        count_d2r_d = count_d2r[d] ;
        pval_d2r_d = pval_d2r[d] ;
        idx_d2r_d = idx_d2r[d] ;
        idx = d % data_block_sz ;
        d2r_table[idx] = malloc(count_d2r_d * sizeof(int)) ;
        d2r_table_d = d2r_table[idx] ;
        for (i = 0 ; i < count_d2r_d ; i++){
            r = rid_map[idx_d2r_d[i]] ;
            if_visit_r[r] = 1 ;
        }

        d2r_table_ct_d = 0 ;
        for (r = 0 ; r < r_block_sz ; r++){
            if (if_visit_r[r] == 1){
                d2r_table_d[d2r_table_ct_d] = r ;
                d2r_table_ct_d += 1 ;
                if_visit_r[r] = 0 ;
            }
        }

        if (d2r_table_ct_d != count_d2r_d)
            printf("error: d2r_table_ct_d != count_d2r_d\n!!") ;

        phi_d = phi[d] ;
        ave_bg_d = ave_bg[d] ;
        data_frame_d = data_frame[d % data_block_sz] ;
        /* frames with phi = 0. are irrelevant when updating vox_val */
        if (phi_d < epsilon)
            continue ;

        for (s = 0 ; s < count_d2r_d ; s++){
            r = d2r_table_d[s] ;
            rid = idx_d2r_d[s] ;
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
                        my_Bval[hkl_id] = realloc(my_Bval[hkl_id], len*sizeof(double)) ;
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

    for (d = idx_start ; d < idx_start + data_block_sz ; d++){
        if (d + dmin > dmax-1)
            break ;
        free(d2r_table[d % data_block_sz]) ;
    }
    free(d2r_table) ;
    free(hkl_npix_mpi) ;
}
