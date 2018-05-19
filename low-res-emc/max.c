#include "emc.h"

void emc(){

    FILE *fp ;
    char out_prob_dir[256], outfile[256] ;
    int i, j, k, d, r, t, qid, rid, ct, idx, idx_start ;
    int dmin, dmax, rank_max, length, proc_id, result ;
    int flag, tag, request_id, r_block_sz, my_r_block_sz ;
    int num_byte, photon_count, *r_block_id_mpi_idx, *if_complete ;
    int **r2hkl_send, **r2hkl_recv, *r2hkl_send_ct, *r2hkl_recv_ct ;
    double w, bv, rescale, dot_product, arg_angle, intens_sum ;
    double logP_offset_d, acceptance, epsilon = 1.e-10 ;
    double t0, t1, t2, dt, t_read, t_logP_offset, t_calprob, t_calAval ;
    double *ave_bg_d, *log_bg_d, *old_peak_val, *new_peak_val ;
    short *data_frame_d ;

    MPI_File fh ;
    MPI_Status status ;
    MPI_Offset offset ;
    MPI_Status *status_arr ;
    MPI_Request *request_arr ;

    sprintf(out_prob_dir, "%s/iter_flag-%d", prob_dir, iter_flag) ;
    if (myid == 0)
        result = mkdir(out_prob_dir, 0777) ;
    MPI_Barrier(MPI_COMM_WORLD) ;

    t0 = MPI_Wtime() ;
    t_read = 0. ;
    t_logP_offset = 0. ;
    t_calprob = 0. ;
    t_calAval = 0. ;
    for (d = 0 ; d < s_num_data ; d++)
        count_d2r[d] = 0 ;

    r2hkl_send = malloc(nproc * sizeof(int *)) ;
    r2hkl_recv = malloc(nproc * sizeof(int *)) ;
    r2hkl_send_ct = malloc(nproc * sizeof(int)) ;
    r2hkl_recv_ct = malloc(nproc * sizeof(int)) ;
    status_arr = malloc(2*(nproc-1) * sizeof(MPI_Status)) ;
    request_arr = malloc(2*(nproc-1) * sizeof(MPI_Request)) ;
    if_complete = malloc(2*(nproc-1) * sizeof(int)) ;

    if (iter_flag % 2 == 1){
        my_Acount = calloc(num_hkl, sizeof(int)) ;
        my_Alen = malloc(num_hkl * sizeof(int)) ;
        my_pval_len = malloc(num_hkl * sizeof(int)) ;
        my_pval_ct = calloc(num_hkl, sizeof(int)) ;
        my_Aval = malloc(num_hkl * sizeof(float *)) ;
        my_Bval = malloc(num_hkl * sizeof(float *)) ;
        my_pval_id = malloc(num_hkl * sizeof(int *)) ;
        for (i = 0 ; i < num_hkl ; i++){
            my_Alen[i] = mem_step_sz ;
            my_pval_len[i] = mem_step_sz ;
            my_Aval[i] = malloc(mem_step_sz * sizeof(float)) ;
            my_Bval[i] = malloc(mem_step_sz * sizeof(float)) ;
            my_pval_id[i] = malloc(mem_step_sz * sizeof(int)) ;
        }
    }

    rank_max = num_data - nproc*(s_num_data - 1) ;
    if (myid < rank_max){
        dmin = myid*s_num_data ;
        dmax = dmin + s_num_data ;
    }
    else{
        dmin = rank_max*s_num_data + (myid - rank_max)*(s_num_data - 1) ;
        dmax = dmin + (s_num_data - 1) ;
    }

    for (i = 0 ; i < iter_data_block ; i++){
        t1 = MPI_Wtime() ;
        memset(r2hkl_send_ct, 0, nproc*sizeof(int)) ;
        memset(r2hkl_recv_ct, 0, nproc*sizeof(int)) ;
        /* number of items to be sent to the j-th processor */
        for (j = 0 ; j < nproc ; j++){
            idx = j*iter_data_block + i ;
            r_block_sz = r_block_sz_mpi[idx] ;
            r_block_id_mpi_idx = r_block_id_mpi[idx] ;
            for (r = 0 ; r < r_block_sz ; r++){
                rid = r_block_id_mpi_idx[r] ;
                proc_id = rid / s_num_rot ;
                if (myid == proc_id)
                    r2hkl_send_ct[j] += 3*num_r2hkl[rid] ;
            }
        }

        request_id = 0 ;
        for (j = 0 ; j < nproc ; j++){
            if (myid == j)
                continue ;
            r2hkl_send[j] = malloc(r2hkl_send_ct[j] * sizeof(int)) ;
            idx = j*iter_data_block + i ;
            r_block_sz = r_block_sz_mpi[idx] ;
            r_block_id_mpi_idx = r_block_id_mpi[idx] ;
            ct = 0 ;
            for (r = 0 ; r < r_block_sz ; r++){
                rid = r_block_id_mpi_idx[r] ;
                proc_id = rid / s_num_rot ;
                if (myid == proc_id){
                    length = 3*num_r2hkl[rid] ;
                    num_byte = length*sizeof(int) ;
                    memcpy(&r2hkl_send[j][ct], r2hkl[rid % s_num_rot], num_byte) ;
                    ct += length ;
                }
            }
            if (ct != r2hkl_send_ct[j])
                printf("error in r2hkl_send_ct!!\n") ;

            tag = myid*nproc + j ;
            MPI_Isend(r2hkl_send[j], r2hkl_send_ct[j], MPI_INT, j, tag, MPI_COMM_WORLD, &request_arr[request_id]) ;
            request_id += 1 ;
        }

        idx = myid*iter_data_block + i ;
        my_r_block_sz = r_block_sz_mpi[idx] ;
        r_block_id_mpi_idx = r_block_id_mpi[idx] ;
        /* number of items received from the j-th processor */
        for (r = 0 ; r < my_r_block_sz ; r++){
            rid = r_block_id_mpi_idx[r] ;
            j = rid / s_num_rot ;
            r2hkl_recv_ct[j] += 3*num_r2hkl[rid] ;
        }

        request_id = nproc-1 ;
        for (j = 0 ; j < nproc ; j++){
            if (myid == j)
                continue ;
            r2hkl_recv[j] = malloc(r2hkl_recv_ct[j] * sizeof(int)) ;
            tag = j*nproc + myid ;
            MPI_Irecv(r2hkl_recv[j], r2hkl_recv_ct[j], MPI_INT, j, tag, MPI_COMM_WORLD, &request_arr[request_id]) ;
            request_id += 1 ;
        }

        /* read in data frames */
        MPI_File_open(MPI_COMM_WORLD, data_file, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh) ;
        idx_start = i*data_block_sz ;
        for (d = idx_start ; d < idx_start + data_block_sz ; d++){
            offset = d + dmin ;
            if (offset > dmax-1)
                break ;
            offset *= full_num_pix*sizeof(short) ;
            idx = d - idx_start ;
            MPI_File_read_at(fh, offset, data_frame[idx], num_pix, MPI_SHORT, &status) ;
        }
        MPI_File_close(&fh) ;

        idx = myid*iter_data_block + i ;
        r_block_id = r_block_id_mpi[idx] ;
        for (r = 0 ; r < my_r_block_sz ; r++){
            rid = r_block_id[r] ;
            r2hkl_mpi[r] = malloc(3*num_r2hkl[rid] * sizeof(int)) ;
            proc_id = rid / s_num_rot ;
            if (myid == proc_id){
                num_byte = 3*num_r2hkl[rid]*sizeof(int) ;
                memcpy(r2hkl_mpi[r], r2hkl[rid % s_num_rot], num_byte) ;
            }
        }

        ct = 0 ;
        length = nproc - 1 ;
        memset(if_complete, 0, 2*length*sizeof(int)) ;
        while (ct < 2*length){
            /* check send request */
            for (j = 0 ; j < nproc ; j++){
                if (j == myid)
                    continue ;
                else if (j < myid)
                    idx = j ;
                else
                    idx = j - 1 ;
                if (if_complete[idx] == 1)
                    continue ;
                MPI_Test(&request_arr[idx], &flag, &status_arr[idx]) ;
                if (!flag)
                    continue ;
                if_complete[idx] = 1 ;
                ct += 1 ;
                free(r2hkl_send[j]) ;
            }
            /* check recv request */
            for (j = 0 ; j < nproc ; j++){
                if (j == myid)
                    continue ;
                else if (j < myid)
                    idx = j + length ;
                else
                    idx = j + length - 1 ;
                if (if_complete[idx] == 1)
                    continue ;
                MPI_Test(&request_arr[idx], &flag, &status_arr[idx]) ;
                if (!flag)
                    continue ;
                t = 0 ;
                for (r = 0 ; r < my_r_block_sz ; r++){
                    rid = r_block_id[r] ;
                    proc_id = rid / s_num_rot ;
                    if (proc_id != j)
                        continue ;
                    num_byte = 3*num_r2hkl[rid]*sizeof(int) ;
                    memcpy(r2hkl_mpi[r], &r2hkl_recv[j][t], num_byte) ;
                    t += 3*num_r2hkl[rid] ;
                }
                if_complete[idx] = 1 ;
                ct += 1 ;
                free(r2hkl_recv[j]) ;
            }
        }

        dt = MPI_Wtime() - t1 ;
        t_read += dt ;
        if (myid == 0)
            printf("data_block_id = %d:\n\treading data takes %.1f sec\n", i, dt) ;

        /* initialize logP_offset, which stores \sum (-w + k log w) */
        if (iter_flag == 1){
            t1 = MPI_Wtime() ;
            for (d = idx_start ; d < idx_start + data_block_sz ; d++){
                if (d + dmin > dmax-1)
                    break ;
                ave_bg_d = ave_bg[d] ;
                log_bg_d = log_bg[d] ;
                data_frame_d = data_frame[d - idx_start] ;
                logP_offset_d = 0. ;
                for (t = 0 ; t < num_pix ; t++){
                    qid = qid_map[t] ;
                    if (qid < 0)
                        continue ;
                    bv = ave_bg_d[qid] ;
                    if (bv < epsilon)
                        continue ;
                    photon_count = data_frame_d[t] ;
                    /* ignore the dead pixels */
                    if (photon_count < 0)
                        continue ;
                    logP_offset_d += photon_count*(log_bg_d[qid] + log_polar[t]) - bv*polar[t] ;
                }
                logP_offset[d] = logP_offset_d ;
            }
            t_logP_offset += MPI_Wtime() - t1 ;
        }

        t1 = MPI_Wtime() ;
        calprob(i) ;
        dt = MPI_Wtime() - t1 ;
        t_calprob += dt ;
        if (myid == 0)
            printf("\tcalprob takes %.2f sec\n", dt) ;

        t1 = MPI_Wtime() ;
        calAval(i) ;
        for (r = 0 ; r < my_r_block_sz ; r++)
            free(r2hkl_mpi[r]) ;
        t_calAval += MPI_Wtime() - t1 ;
    }
    free(r2hkl_send) ;
    free(r2hkl_recv) ;
    free(r2hkl_send_ct) ;
    free(r2hkl_recv_ct) ;
    free(status_arr) ;
    free(request_arr) ;
    free(if_complete) ;

    dt = MPI_Wtime() - t0 ;
    t0 = MPI_Wtime() ;
    if (myid == 0){
        printf("expand, calprob & normalize_prob take %.0f sec\n", dt) ;
        if (iter_flag == 1)
            printf("calculating logP_offset takes %.0f sec\n", t_logP_offset) ;
        printf("t_calprob = %.0f sec, t_calAval = %.0f sec\n", t_calprob, t_calAval) ;
    }

    if (iter_flag % 2 == 0){

        idx = 0 ;
        w = 0. ;
        dot_product = 0. ;
        for (d = dmin ; d < dmax ; d++){
            dot_product += phi[idx]*updated_phi[idx] ;
            w += phi[idx]*phi[idx] ;
            idx += 1 ;
        }
        MPI_Allreduce(MPI_IN_PLACE, &dot_product, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;
        MPI_Allreduce(MPI_IN_PLACE, &w, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;

        dot_product /= sqrt(w) ;
        w = 0. ;

        idx = 0 ;
        rescale = 0. ;
        acceptance = 0. ;
        for (d = dmin ; d < dmax ; d++){
            phi[idx] = updated_phi[idx] ;
            w += phi[idx]*phi[idx] ;
            rescale += updated_phi[idx] ;
            if (updated_phi[idx] > epsilon)
                acceptance += 1 ;
            idx += 1 ;
        }
        MPI_Allreduce(MPI_IN_PLACE, &acceptance, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;
        MPI_Allreduce(MPI_IN_PLACE, &rescale, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;
        MPI_Allreduce(MPI_IN_PLACE, &w, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;
        arg_angle = acos(dot_product / sqrt(w))*180/3.1416 ;

        idx = 0 ;
        rescale = num_data/rescale ;
        for (d = dmin ; d < dmax ; d++){
            phi[idx] *= rescale ;
            idx += 1 ;
        }

        rescale = 1./rescale ;
        for (i = 0 ; i < num_hkl ; i++){
            for (j = 0 ; j < num_hkl2vox_id[i] ; j++)
                intens1[i][j] *= rescale ;
        }

        sprintf(outfile, "%s/total-phi.dat", out_prob_dir) ;
        if (myid == 0){
            fp = fopen(outfile, "w") ;
            fclose(fp) ;
        }

        for (i = 0 ; i < nproc ; i++){
            if (myid == i){
                idx = 0 ;
                fp = fopen(outfile, "a") ;
                for (d = dmin ; d < dmax ; d++){
                    fprintf(fp, "%12.5e\n", phi[idx]) ;
                    idx += 1 ;
                }
                fclose(fp) ;
            }
            MPI_Barrier(MPI_COMM_WORLD) ;
        }

        if (myid == 0){
            acceptance *= 100. / num_data ;
            printf("acceptance rate = %.4f percent\n", acceptance) ;
            fp = fopen("EMC.log", "a") ;
            fprintf(fp, "iter = %d:\n", (iter_flag+1)/2) ;
            fprintf(fp, "\tupdating phi takes %.0f sec, arg_angle = %12.5e degrees\n\n", dt, arg_angle) ;
            fclose(fp) ;
        }

        return ;
    }

    if (myid == 0){
        intens_sum = 0. ;
        for (i = 0 ; i < num_hkl ; i++){
            for (j = 0 ; j < num_hkl2vox_id[i] ; j++)
                intens_sum += intens1[i][j] ;
        }
        printf("sum of intens1 = %12.5e\n", intens_sum) ;
    }

    old_peak_val = calloc(num_hkl, sizeof(double)) ;
    for (i = 0 ; i < num_hkl ; i++){
        w = 0. ;
        for (j = 0 ; j < num_hkl2vox_id[i] ; j++)
            w += intens1[i][j] ;
        old_peak_val[i] = w ;
    }

    memset(if_visit_r, 0, num_rot*sizeof(int)) ;
    for (d = 0 ; d < dmax-dmin ; d++){
        if (phi[d] < epsilon)
            continue ;
        for (i = 0 ; i < count_d2r[d] ; i++)
            if_visit_r[idx_d2r[d][i]] = 1 ;
    }
    MPI_Allreduce(MPI_IN_PLACE, if_visit_r, num_rot, MPI_INT, MPI_MAX, MPI_COMM_WORLD) ;

    sprintf(outfile, "%s/intens-mean-std.dat", out_prob_dir) ;
    if (myid == 0){
        fp = fopen(outfile, "w") ;
        fclose(fp) ;
    }

    MPI_File_open(MPI_COMM_WORLD, peak2r_file, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh) ;
    offset = sizeof(int) ;
    /* start maximize to update vox_val (called in compress) */
    for (k = 0 ; k < num_hkl ; k++){
        t1 = MPI_Wtime() ;
        length = num_hkl2r[k]*3 ;
        MPI_File_read_at(fh, offset, &hkl2r_mpi[0], length, MPI_INT, &status) ;
        offset += (length+1)*sizeof(int) ;
        t2 = MPI_Wtime() ;
        if (myid == 0)
            printf("Initializing Emat takes = %.4f sec\n", t2-t1) ;

        t1 = MPI_Wtime() ;
        compress(k) ;
        t2 = MPI_Wtime() ;
        if (myid == 0)
            printf("elapsed time = %.3f sec\n\n", t2-t1) ;
    }
    MPI_File_close(&fh) ;

    for (k = 0 ; k < num_hkl ; k++){
        free(my_Aval[k]) ;
        free(my_Bval[k]) ;
        free(my_pval_id[k]) ;
    }
    free(my_Aval) ;
    free(my_Bval) ;
    free(my_pval_id) ;
    free(my_pval_ct) ;
    free(my_pval_len) ;
    free(my_Alen) ;
    free(my_Acount) ;

    if (myid == 0)
        sym_intens() ;
    for (i = 0 ; i < num_hkl ; i++)
        MPI_Bcast(intens1[i], num_hkl2vox_id[i], MPI_DOUBLE, 0, MPI_COMM_WORLD) ;

    t2 = MPI_Wtime() ;
    if (myid == 0)
        printf("maximize & compress take %.0f sec\n", t2-t0) ;

    new_peak_val = calloc(num_hkl, sizeof(double)) ;
    for (i = 0 ; i < num_hkl ; i++){
        w = 0. ;
        for (j = 0 ; j < num_hkl2vox_id[i] ; j++)
            w += intens1[i][j] ;
        new_peak_val[i] = w ;
    }

    dot_product = 0. ;
    for (i = 0 ; i < num_hkl ; i++)
        dot_product += new_peak_val[i]*old_peak_val[i] ;

    w = 0. ;
    for (i = 0 ; i < num_hkl ; i++)
        w += new_peak_val[i]*new_peak_val[i] ;
    dot_product /= sqrt(w) ;

    w = 0. ;
    for (i = 0 ; i < num_hkl ; i++)
        w += old_peak_val[i]*old_peak_val[i] ;
    dot_product /= sqrt(w) ;
    arg_angle = acos(dot_product)*180/3.1416 ;

    dt += (t2-t0) ;
    if (myid == 0){
        print_intens() ;
        fp = fopen("EMC.log", "a") ;
        fprintf(fp, "iter = %d:\n", (iter_flag+1)/2) ;
        fprintf(fp, "\tupdating vox_val takes = %.0f sec, arg_angle = %12.5e degrees\n", dt, arg_angle) ;
        fclose(fp) ;
        printf("arg_angle = %12.5e degrees, total elapsed time = %.0f sec\n", arg_angle, dt) ;
    }
    
    free(old_peak_val) ;
    free(new_peak_val) ;
}


void maximize( int cur_num_hkl2r, int hkl_id ){

    FILE *fp ;
    char outfile[256] ;
    int i, j, k, k0, d, r, t, rid, pid, pix_ct, vox_ct, vox_id, idx, nvox ;
    int my_pval_ct_hkl, total_pval_ct_hkl, my_start_id, rcount, s_rcount, entry_count ;
    int if_overshoot, search_ct, *cur_Acount, *if_skip_tomo, max_iter = 200 ;
    int *proc_pval_ct_hkl, *mpi_rid, *proc_Acount, *my_pval_id_hkl, *entry_map ;
    double pval, w, val, phi_val, norm, dv, diff, lambda, grad ;
    double S_a, S_b, S_c, S_d, w_a, w_b, w_c, w_d, inv_tau ;
    double **cur_Cval, **cur_Dval, **cur_Pval, tolerance = 1.e-3, epsilon = 1.e-10 ;
    double *mpi_pval, *mpi_phi, *updated_tomo, *vox_val, *inter_weight ;
    float *my_Aval_hkl, *mpi_Aval, **cur_Aval ;
    float *my_Bval_hkl, *mpi_Bval, **cur_Bval ;
    long long m, i0, total_Acount ;
    E_matrix Emat_r ;

    /* load balancing between processors */
    proc_pval_ct_hkl = calloc(nproc, sizeof(int)) ;
    my_pval_ct_hkl = my_pval_ct[hkl_id]/2 ;
    proc_pval_ct_hkl[myid] = my_pval_ct_hkl ;
    MPI_Allreduce(MPI_IN_PLACE, proc_pval_ct_hkl, nproc, MPI_INT, MPI_MAX, MPI_COMM_WORLD) ;

    total_pval_ct_hkl = 0 ;
    my_start_id = 0 ;
    for (i = 0 ; i < nproc ; i++){
        total_pval_ct_hkl += proc_pval_ct_hkl[i] ;
        if (i < myid)
            my_start_id += proc_pval_ct_hkl[i] ;
    }

    if (total_pval_ct_hkl == 0){
        if (myid == 0)
            printf("warning: total_pval_ct_hkl = 0; this peak is not intersected!!\n") ;
        nvox = num_hkl2vox_id[hkl_id] ;
        for (i = 0 ; i < nvox ; i++)
            intens1[hkl_id][i] = 0. ;
        free(proc_pval_ct_hkl) ;
        return ;
    }

    mpi_pval = malloc(total_pval_ct_hkl * sizeof(double)) ;
    mpi_phi = malloc(total_pval_ct_hkl * sizeof(double)) ;
    mpi_rid = malloc(total_pval_ct_hkl * sizeof(int)) ;
    proc_Acount = calloc(nproc, sizeof(int)) ;

    my_pval_id_hkl = my_pval_id[hkl_id] ;
    for (j = 0 ; j < my_pval_ct_hkl ; j++){
        idx = 2*j ;
        d = my_pval_id_hkl[idx] ;
        i = my_pval_id_hkl[idx+1] ;
        t = my_start_id + j ;
        rid = idx_d2r[d][i] ;
        mpi_rid[t] = rid ;
        mpi_pval[t] = pval_d2r[d][i] ;
        mpi_phi[t] = phi[d] ;
        proc_Acount[myid] += Emat[rid].npix ;
    }
    MPI_Allreduce(MPI_IN_PLACE, proc_Acount, nproc, MPI_INT, MPI_MAX, MPI_COMM_WORLD) ;

    if (proc_Acount[myid] != my_Acount[hkl_id]){
        printf("error: myid = %d, proc_Acount[myid] != my_Acount[hkl_id]!!\n", myid) ;
        return ;
    }

    my_start_id = 0 ;
    for (i = 0 ; i < nproc ; i++){
        MPI_Bcast(&mpi_rid[my_start_id], proc_pval_ct_hkl[i], MPI_INT, i, MPI_COMM_WORLD) ;
        MPI_Bcast(&mpi_pval[my_start_id], proc_pval_ct_hkl[i], MPI_DOUBLE, i, MPI_COMM_WORLD) ;
        MPI_Bcast(&mpi_phi[my_start_id], proc_pval_ct_hkl[i], MPI_DOUBLE, i, MPI_COMM_WORLD) ;
        my_start_id += proc_pval_ct_hkl[i] ;
    }

    my_start_id = 0 ;
    total_Acount = 0 ;
    for (i = 0 ; i < nproc ; i++){
        total_Acount += proc_Acount[i] ;
        if (i < myid)
            my_start_id += proc_Acount[i] ;
    }

    if (total_Acount > int_max){
        if (myid == 0)
            printf("warning: total_Acount overflows int32!!\n") ;
        return ;
    }
    if (myid == 0)
        printf("total_Acount = %lld\n", total_Acount) ;

    my_Aval_hkl = my_Aval[hkl_id] ;
    my_Bval_hkl = my_Bval[hkl_id] ;
    mpi_Aval = malloc(total_Acount * sizeof(float)) ;
    mpi_Bval = malloc(total_Acount * sizeof(float)) ;
    for (j = 0 ; j < proc_Acount[myid] ; j++){
        idx = my_start_id + j ;
        mpi_Aval[idx] = my_Aval_hkl[j] ;
        mpi_Bval[idx] = my_Bval_hkl[j] ;
    }

    my_start_id = 0 ;
    for (i = 0 ; i < nproc ; i++){
        MPI_Bcast(&mpi_Aval[my_start_id], proc_Acount[i], MPI_FLOAT, i, MPI_COMM_WORLD) ;
        MPI_Bcast(&mpi_Bval[my_start_id], proc_Acount[i], MPI_FLOAT, i, MPI_COMM_WORLD) ;
        my_start_id += proc_Acount[i] ;
    }

    rcount = 0 ;
    for (r = 0 ; r < cur_num_hkl2r ; r++){
        rid = hkl2r_mpi[3*r] ;
        if (if_visit_r[rid] == 0)
            continue ;
        rid_map[rid] = rcount ;
        rcount += 1 ;
    }

    if (rcount % nproc == 0)
        s_rcount = rcount / nproc ;
    else
        s_rcount = rcount / nproc + 1 ;

    entry_count = 0 ;
    entry_map = malloc(s_rcount * sizeof(int)) ;
    for (r = 0 ; r < cur_num_hkl2r ; r++){
        rid = hkl2r_mpi[3*r] ;
        if (if_visit_r[rid] == 0)
            continue ;
        if (rid_map[rid] / s_rcount == myid){
            entry_map[rid_map[rid] % s_rcount] = entry_count ;
            entry_count += Emat[rid].npix ;
        }
    }

    cur_Acount = calloc(entry_count, sizeof(int)) ;
    cur_Aval = malloc(entry_count * sizeof(float *)) ;
    cur_Bval = malloc(entry_count * sizeof(float *)) ;
    cur_Cval = malloc(entry_count * sizeof(double *)) ;
    cur_Dval = malloc(entry_count * sizeof(double *)) ;
    cur_Pval = malloc(entry_count * sizeof(double *)) ;

    for (j = 0 ; j < total_pval_ct_hkl ; j++){
        rid = mpi_rid[j] ;
        pix_ct = Emat[rid].npix ;
        if (rid_map[rid] / s_rcount == myid){
            k0 = entry_map[rid_map[rid] % s_rcount] ;
            for (t = 0 ; t < pix_ct ; t++)
                cur_Acount[k0 + t] += 1 ;
        }
    }

    for (i = 0 ; i < entry_count ; i++){
        cur_Aval[i] = malloc(cur_Acount[i] * sizeof(float)) ;
        cur_Bval[i] = malloc(cur_Acount[i] * sizeof(float)) ;
        cur_Cval[i] = malloc(cur_Acount[i] * sizeof(double)) ;
        cur_Dval[i] = malloc(cur_Acount[i] * sizeof(double)) ;
        cur_Pval[i] = malloc(cur_Acount[i] * sizeof(double)) ;
        cur_Acount[i] = 0 ;
    }

    m = 0 ;
    for (j = 0 ; j < total_pval_ct_hkl ; j++){
        rid = mpi_rid[j] ;
        Emat_r = Emat[rid] ;
        pix_ct = Emat_r.npix ;
        if (rid_map[rid] / s_rcount != myid){
            m += pix_ct ;
            continue ;
        }
        pval = mpi_pval[j] ;
        phi_val = mpi_phi[j] ;
        val = pval*phi_val ;
        k0 = entry_map[rid_map[rid] % s_rcount] ;
        for (t = 0 ; t < pix_ct ; t++){
            idx = k0 + t ;
            i0 = m + t ;
            pid = Emat_r.pix_id[t] ;
            cur_Aval[idx][cur_Acount[idx]] = mpi_Aval[i0] ;
            cur_Bval[idx][cur_Acount[idx]] = mpi_Bval[i0] ;
            cur_Cval[idx][cur_Acount[idx]] = val*polar[pid] ;
            cur_Dval[idx][cur_Acount[idx]] = log(phi_val*polar[pid]) ;
            cur_Pval[idx][cur_Acount[idx]] = pval ;
            cur_Acount[idx] += 1 ;
        }
        m += pix_ct ;
    }

    inv_tau = 1./golden_ratio ;
    if_skip_tomo = malloc(entry_count * sizeof(int)) ;
    updated_tomo = malloc(entry_count * sizeof(double)) ;
    for (i = 0 ; i < entry_count ; i++){
        /* minimize Sij using line search */
        w_a = 0. ;
        S_a = 0. ;
        grad = 0. ;
        if_skip_tomo[i] = 1 ;
        for (k = 0 ; k < cur_Acount[i] ; k++){
            if (cur_Aval[i][k] < 0. || cur_Bval[i][k] < 0.)
                continue ;
            w = cur_Bval[i][k] + w_a ;
            S_a += cur_Cval[i][k]*w - cur_Aval[i][k]*(log(w) + cur_Dval[i][k]) ;
            grad += cur_Cval[i][k] - cur_Aval[i][k]/w ;
            if_skip_tomo[i] = 0 ;
        }

        if (if_skip_tomo[i] == 1)
            continue ;
        if (fabs(grad) < epsilon){
            updated_tomo[i] = 0. ;
            continue ;
        }

        if_skip_tomo[i] = 1 ;
        for (k = 0 ; k < cur_Acount[i] ; k++){
            if (cur_Aval[i][k] < 0. || cur_Bval[i][k] < 0.)
                continue ;
            if (cur_Aval[i][k] > epsilon)
                if_skip_tomo[i] = 0 ;
        }

        /* when all photon counts = 0 */
        if (if_skip_tomo[i] == 1){
            for (k = 0 ; k < cur_Acount[i] ; k++){
                if (cur_Aval[i][k] < 0. || cur_Bval[i][k] < 0.)
                    continue ;
                updated_tomo[i] = -cur_Bval[i][k] ;
                break ;
            }
            for (k = 0 ; k < cur_Acount[i] ; k++){
                if (cur_Aval[i][k] < 0. || cur_Bval[i][k] < 0.)
                    continue ;
                if (updated_tomo[i] < -cur_Bval[i][k])
                    updated_tomo[i] = -cur_Bval[i][k] ;
            }
            if_skip_tomo[i] = 0 ;
            continue ;
        }

        /* bracket the minimum */
        S_b = S_a ;
        lambda = 1. ;
        if_overshoot = 0 ;
        while (S_b <= S_a && if_overshoot == 0){
            w_b = w_a - lambda*grad ;
            S_b = 0. ;
            for (k = 0 ; k < cur_Acount[i] ; k++){
                if (cur_Aval[i][k] < 0. || cur_Bval[i][k] < 0.)
                    continue ;
                w = cur_Bval[i][k] + w_b ;
                if (w < epsilon){
                    if_overshoot = 1 ;
                    break ;
                }
                S_b += cur_Cval[i][k]*w - cur_Aval[i][k]*(log(w) + cur_Dval[i][k]) ;
            }
            if (if_overshoot == 1 || S_b > S_a)
                break ;
            lambda *= 2. ;
        }

        if (if_overshoot == 1){
            for (k = 0 ; k < cur_Acount[i] ; k++){
                if (cur_Aval[i][k] < 0. || cur_Bval[i][k] < 0.)
                    continue ;
                w = cur_Bval[i][k] + w_a ;
                if (grad > 0.){
                    val = w/grad ;
                    if (lambda > val)
                        lambda = val ;
                }
            }
        }

        /* golden section search */
        w_b = w_a - lambda*grad ;
        dv = w_b - w_a ;
        diff = sqrt(dv*dv) ;
        w = (w_b + w_a)/2. ;
        norm = sqrt(w*w) ;
        search_ct = 0 ;
        while (diff > tolerance*norm && search_ct < max_iter){
            dv = (w_b - w_a)*inv_tau ;
            w_c = w_b - dv ;
            w_d = w_a + dv ;

            S_c = 0. ;
            S_d = 0. ;
            for (k = 0 ; k < cur_Acount[i] ; k++){
                if (cur_Aval[i][k] < 0. || cur_Bval[i][k] < 0.)
                    continue ;
                w = cur_Bval[i][k] + w_c ;
                S_c += cur_Cval[i][k]*w - cur_Aval[i][k]*(log(w) + cur_Dval[i][k]) ;
                w = cur_Bval[i][k] + w_d ;
                S_d += cur_Cval[i][k]*w - cur_Aval[i][k]*(log(w) + cur_Dval[i][k]) ;
            }

            if (S_c < S_d)
                w_b = w_d ;
            else
                w_a = w_c ;

            dv = w_b - w_a ;
            diff = sqrt(dv*dv) ;
            w = (w_b + w_a)/2. ;
            norm = sqrt(w*w) ;
            search_ct += 1 ;
        }

        if (search_ct == 0)
            updated_tomo[i] = (w_a + w_b)/2. ;
        else
            updated_tomo[i] = (w_c + w_d)/2. ;
    }

    /* calculate the variance of wij */
    int if_all_zero ;
    double inv_w, Cyy, inv_Cyy, Jval, *var_w, *Cxy, *cov_x ;

    var_w = calloc(entry_count, sizeof(double)) ;
    for (i = 0 ; i < entry_count ; i++){
        if (if_skip_tomo[i] == 1)
            continue ;
        w = updated_tomo[i] ;

        Cyy = 0. ;
        if_all_zero = 0 ;
        cov_x = calloc(cur_Acount[i], sizeof(double)) ;
        Cxy = calloc(cur_Acount[i], sizeof(double)) ;
        for (k = 0 ; k < cur_Acount[i] ; k++){
            if (cur_Aval[i][k] < 0. || cur_Bval[i][k] < 0.){
                cov_x[k] = -1 ;
                continue ;
            }

            /* all photons are zero */
            if (cur_Bval[i][k] + w < epsilon){
                if_all_zero = 1 ;
                break ;
            }
            inv_w = 1. / (cur_Bval[i][k] + w) ;
            Cyy += cur_Aval[i][k]*inv_w*inv_w ;
            Cxy[k] = -cur_Pval[i][k]*inv_w ;
            cov_x[k] = cur_Aval[i][k]/cur_Pval[i][k] ;
        }

        if (if_all_zero == 1){
            free(cov_x) ;
            free(Cxy) ;
            continue ;
        }

        inv_Cyy = 1./Cyy ;
        for (k = 0 ; k < cur_Acount[i] ; k++){
            if (cov_x[k] < 0.)
                continue ;
            Jval = -Cxy[k]*inv_Cyy ;
            var_w[i] += cov_x[k]*Jval*Jval ;
        }

        free(cov_x) ;
        free(Cxy) ;
    }

    nvox = num_hkl2vox_id[hkl_id] ;
    vox_val = calloc(nvox, sizeof(double)) ;
    inter_weight = calloc(nvox, sizeof(double)) ;

    for (j = 0 ; j < total_pval_ct_hkl ; j++){
        rid = mpi_rid[j] ;
        Emat_r = Emat[rid] ;
        pix_ct = Emat_r.npix ;
        if (rid_map[rid] / s_rcount != myid)
            continue ;
        pval = mpi_pval[j] ;
        phi_val = mpi_phi[j] ;
        val = pval*phi_val ;
        idx = entry_map[rid_map[rid] % s_rcount] ;
        for (t = 0 ; t < pix_ct ; t++){
            if (if_skip_tomo[idx] == 1){
                idx += 1 ;
                continue ;
            }
            k0 = 8*t ;
            vox_ct = Emat_r.vox_ct[t] ;
            for (k = 0 ; k < vox_ct ; k++){
                vox_id = Emat_r.vox_id[k0] ;
                vox_val[vox_id] += val*Emat_r.weight[k0]*updated_tomo[idx] ;
                inter_weight[vox_id] += val*Emat_r.weight[k0] ;
                k0 += 1 ;
            }
            idx += 1 ;
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, vox_val, nvox, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;
    MPI_Allreduce(MPI_IN_PLACE, inter_weight, nvox, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;
    
    double var_hkl, I_hkl, *alpha, *inv_inter_weight ;
    inv_inter_weight = calloc(nvox, sizeof(double)) ;
    
    for (i = 0 ; i < nvox ; i++){
        if (inter_weight[i] > epsilon){
            inv_inter_weight[i] = 1. / inter_weight[i] ;
            vox_val[i] *= inv_inter_weight[i] ;
        }
    }

    alpha = calloc(entry_count, sizeof(double)) ;
    for (j = 0 ; j < total_pval_ct_hkl ; j++){
        rid = mpi_rid[j] ;
        Emat_r = Emat[rid] ;
        pix_ct = Emat_r.npix ;
        if (rid_map[rid] / s_rcount != myid)
            continue ;
        pval = mpi_pval[j] ;
        phi_val = mpi_phi[j] ;
        val = pval*phi_val ;
        idx = entry_map[rid_map[rid] % s_rcount] ;
        for (t = 0 ; t < pix_ct ; t++){
            if (if_skip_tomo[idx] == 1){
                idx += 1 ;
                continue ;
            }
            k0 = 8*t ;
            vox_ct = Emat_r.vox_ct[t] ;
            for (k = 0 ; k < vox_ct ; k++){
                vox_id = Emat_r.vox_id[k0] ;
                alpha[idx] += val*Emat_r.weight[k0]*inv_inter_weight[vox_id] ;
                k0 += 1 ;
            }
            idx += 1 ;
        }
    }

    var_hkl = 0. ;
    I_hkl = 0. ;
    for (i = 0 ; i < entry_count ; i++){
        if (if_skip_tomo[i] == 1)
            continue ;
        var_hkl += alpha[i]*alpha[i]*var_w[i] ;
        I_hkl += alpha[i]*updated_tomo[i] ;
    }

    MPI_Allreduce(MPI_IN_PLACE, &var_hkl, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;
    MPI_Allreduce(MPI_IN_PLACE, &I_hkl, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;

    diff = 0. ;
    w = 0. ;
    for (i = 0 ; i < nvox ; i++){
       diff += vox_val[i] ;
       w += intens1[hkl_id][i] ;
       intens1[hkl_id][i] = vox_val[i] ;
    }

    diff -= w ;
    if (fabs(w) > epsilon){
        diff = fabs(diff)/fabs(w) ;
        if (myid == 0)
            printf("change of vox_val = %12.5e, ", diff) ;
    }

    if (myid == 0){
        printf("I_hkl = %12.5e, std_hkl = %12.5e\n", I_hkl, sqrt(var_hkl)) ;
        sprintf(outfile, "%s/iter_flag-%d/intens-mean-std.dat", prob_dir, iter_flag) ;
        fp = fopen(outfile, "a") ;
        idx = 3*hkl_id ;
        fprintf(fp, "%d %d %d %d %12.5e %12.5e\n", hkl_id, hkl_table[idx], \
            hkl_table[idx+1], hkl_table[idx+2], I_hkl, sqrt(var_hkl)) ;
        fclose(fp) ;
    }

    free(proc_pval_ct_hkl) ;
    free(mpi_pval) ;
    free(mpi_phi) ;
    free(mpi_rid) ;
    free(proc_Acount) ;
    free(mpi_Aval) ;
    free(mpi_Bval) ;
    free(entry_map) ;
    for (i = 0 ; i < entry_count ; i++){
        free(cur_Aval[i]) ;
        free(cur_Bval[i]) ;
        free(cur_Cval[i]) ;
        free(cur_Dval[i]) ;
        free(cur_Pval[i]) ;
    }
    free(cur_Aval) ;
    free(cur_Bval) ;
    free(cur_Cval) ;
    free(cur_Dval) ;
    free(cur_Pval) ;
    free(cur_Acount) ;
    free(if_skip_tomo) ;
    free(updated_tomo) ;
    free(vox_val) ;
    free(inter_weight) ;
    free(var_w) ;
    free(inv_inter_weight) ;
    free(alpha) ;
}


void sym_intens(){

    FILE *fp ;
    char sym_file[256] ;
    int i, j, k, s, i0, j0, k0, hval, kval, lval ;
    int x, y, z, x0, y0, z0, tx, ty, tz, hkl_id, idx ;
    int Nx_half, Ny_half, Nz_half, gw_ceil, vec[3] ;
    int num_sym_op, max_hkl_len, max_num_vox, indices[3] ;
    double (*sym_op)[9], **cpy_intens, **inter_weight ;
    double cx, cy, cz, fx, fy, fz, qx, qy, qz, w ;
    double rescale, rot_vec[3], rvec[3], coefft[3] ;
    long long vox_id ;

    sprintf(sym_file, "%s/aux/sym-op.dat", home_dir) ;
    fp = fopen(sym_file, "r") ;
    fscanf(fp, "%d", &num_sym_op) ;
    sym_op = malloc(num_sym_op * sizeof(*sym_op)) ;
    for (s = 0 ; s < num_sym_op ; s++){
        for (i = 0 ; i < 9 ; i++)
            fscanf(fp, "%lf", &sym_op[s][i]) ;
    }
    fclose(fp) ;

    Nx_half = (Nx-1)/2 ;
    Ny_half = (Ny-1)/2 ;
    Nz_half = (Nz-1)/2 ;
    gw_ceil = ((int) ceil(gw)) ;
    max_hkl_len = hlen*klen*llen ;
    max_num_vox = (2*gw_ceil+3)*(2*gw_ceil+3)*(2*gw_ceil+3) ;
    cpy_intens = malloc(max_hkl_len * sizeof(double *)) ;
    inter_weight = malloc(max_hkl_len * sizeof(double *)) ;
    for (i = 0 ; i < max_hkl_len ; i++){
        cpy_intens[i] = calloc(max_num_vox, sizeof(double)) ;
        inter_weight[i] = calloc(max_num_vox, sizeof(double)) ;
    }

    rescale = VN/min_rcell_sz ;
    for (s = 0 ; s < num_sym_op ; s++){
        for (k = 0 ; k < num_hkl ; k++){
            idx = 3*k ;
            for (i0 = 0 ; i0 < 3 ; i0++)
                indices[i0] = hkl_table[idx+i0] ;

            for (i0 = 0 ; i0 < 3 ; i0++){
                rvec[i0] = 0. ;
                for (j0 = 0 ; j0 < 3 ; j0++)
                    rvec[i0] += rlatt_vec[i0][j0]*indices[j0] ;
            }

            for (i0 = 0 ; i0 < 3 ; i0++){
                idx = 3*i0 ;
                rot_vec[i0] = sym_op[s][idx]*rvec[0] ;
                for (j0 = 1 ; j0 < 3 ; j0++)
                    rot_vec[i0] += sym_op[s][idx+j0]*rvec[j0] ;
            }

            for (i0 = 0 ; i0 < 3 ; i0++){
                coefft[i0] = 0. ;
                for (j0 = 0 ; j0 < 3 ; j0++)
                    coefft[i0] += proj[i0][j0]*rot_vec[j0] ;
            }

            for (i0 = 0 ; i0 < 3 ; i0++)
                indices[i0] = ((int) round(coefft[i0])) ;

            hval = indices[0] ;
            kval = indices[1] ;
            lval = indices[2] ;
            hkl_id = (hval+hmax)*klen*llen + (kval+kmax)*llen + (lval+lmax) ;

            for (i0 = 0 ; i0 < 3 ; i0++){
                rvec[i0] = 0. ;
                for (j0 = 0 ; j0 < 3 ; j0++)
                    rvec[i0] += rlatt_vec[i0][j0]*indices[j0] ;
            }

            qx = rvec[0]*rescale ;
            qy = rvec[1]*rescale ;
            qz = rvec[2]*rescale ;

            tx = ((int) round(qx)) - box_half_len ;
            ty = ((int) round(qy)) - box_half_len ;
            tz = ((int) round(qz)) - box_half_len ;

            for (j = 0 ; j < num_hkl2vox_id[k] ; j++){
                vox_id = hkl2vox_id[k][j] ;
                vec[0] = vox_id / NyNz - Nx_half ;
                vec[1] = (vox_id % NyNz) / Nz - Ny_half ;
                vec[2] = vox_id % Nz - Nz_half ;
                for (i0 = 0 ; i0 < 3 ; i0++){
                    idx = 3*i0 ;
                    rot_vec[i0] = sym_op[s][idx]*vec[0] ;
                    for (j0 = 1 ; j0 < 3 ; j0++)
                        rot_vec[i0] += sym_op[s][idx+j0]*vec[j0] ;
                }

                x = rot_vec[0] - tx ;
                y = rot_vec[1] - ty ;
                z = rot_vec[2] - tz ;

                fx = rot_vec[0] - floor(rot_vec[0]) ;
                fy = rot_vec[1] - floor(rot_vec[1]) ;
                fz = rot_vec[2] - floor(rot_vec[2]) ;

                for (i0 = 0 ; i0 < 2 ; i0++){
                for (j0 = 0 ; j0 < 2 ; j0++){
                for (k0 = 0 ; k0 < 2 ; k0++){
                    x0 = x + i0 ;
                    if (abs(x0 - box_half_len) > box_half_len)
                        continue ;
                    y0 = y + j0 ;
                    if (abs(y0 - box_half_len) > box_half_len)
                        continue ;
                    z0 = z + k0 ;
                    if (abs(z0 - box_half_len) > box_half_len)
                        continue ;
                    idx = x0*box_len2 + y0*box_len + z0 ;
                    if (i0 == 0)
                        cx = 1. - fx ;
                    else
                        cx = fx ;
                    if (j0 == 0)
                        cy = 1. - fy ;
                    else
                        cy = fy ;
                    if (k0 == 0)
                        cz = 1. - fz ;
                    else
                        cz = fz ;

                    w = cx*cy*cz ;
                    cpy_intens[hkl_id][idx] += w*intens1[k][j] ;
                    inter_weight[hkl_id][idx] += w ;
                }
                }
                }
            }
        }
    }

    for (k = 0 ; k < num_hkl ; k++){
        idx = 3*k ;
        hval = hkl_table[idx] ;
        kval = hkl_table[idx+1] ;
        lval = hkl_table[idx+2] ;
        hkl_id = (hval+hmax)*klen*llen + (kval+kmax)*llen + (lval+lmax) ;
        for (j = 0 ; j < max_num_vox ; j++){
            j0 = box_vox_id[k][j] ;
            if (j0 < 0 || inter_weight[hkl_id][j] == 0.)
                continue ;
            intens1[k][j0] = cpy_intens[hkl_id][j] / inter_weight[hkl_id][j] ;
        }
    }

    for (i = 0 ; i < max_hkl_len ; i++){
        free(cpy_intens[i]) ;
        free(inter_weight[i]) ;
    }
    free(cpy_intens) ;
    free(inter_weight) ;
    free(sym_op) ;
}
