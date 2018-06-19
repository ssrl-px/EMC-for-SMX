/*
make the azimuthally averaged background
& label the outlier pixels and Bragg peak candidates
data frames are read in the cbf format

usage:
mpirun -np [number of processors] ./ave_bg [path of the config_file] > run.log &

needs:
config_file, cbflist.dat, outlierlist.dat, radiallist.dat
peaklist.dat, pix-map.dat, rec-vectors.dat, rec2pix-map.dat

makes:
ave_bg_files, outlier_files, peak_files, data-info.dat

*/

#include "make-bg.h"

const double cdf_thres = 0.99999 ;
const double vmax = 200. ;
const double dv = 0.01 ;
const int num_iter = 5 ;
const int max_num_peak_on_ring = 5 ;
/* threshold for ring detection */
const double bg_thres = 1.05 ;

char config_file[256], home_dir[256] ;
char pixmap_file[256], rvec_file[256] ;
char rec2pixmap_file[256] ;

int num_row, num_col, num_data, qlen, total_pix ;
int nproc, myid, qmax, qmin, num_pix, len_val, outlier_ct ;
int *pix_map, *rec2pix_map, *qid_map, *queue ;
int *peak_id, *peak_val, *peak_label, *peak_map, *peak_list ;
int *det, *radial_ct, *count_thres, (*peak_location)[2] ;
double detd, wl, px, hot_pix_thres, dq, cx, cy ;
double *radial_val, *radial_weight, *ave_bg, (*pix)[4] ;
filename *cbf_files, *outlierfiles, *radialfiles, *peakfiles ;
frame_info *data_info ;

int main(int argc, char *argv[]){

	MPI_Init(&argc, &argv) ;
	MPI_Comm_size(MPI_COMM_WORLD, &nproc) ;
	MPI_Comm_rank(MPI_COMM_WORLD, &myid) ;

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

    if (myid == 0)
        printf("setup takes %.0f sec\n", MPI_Wtime()-t1) ;

    MPI_Barrier(MPI_COMM_WORLD) ;
    make_bg() ;
    free_mem() ;

    MPI_Barrier(MPI_COMM_WORLD) ;
    t2 = MPI_Wtime() ;
    if (myid == 0)
        printf("elapsed time = %.0f sec\n", t2-t1) ;

    MPI_Finalize() ;
    return 0 ;
}


void make_bg(){

    FILE *fp ;
    int i, d, t, dmin, dmax, pid, qid, idx ;
    int s_num_data, photon_count, thres, *data_frame ;
    double avg_val, inv_dv = 1./dv ;

    data_frame = malloc(num_pix * sizeof(int)) ;
    if (num_data % nproc == 0)
        s_num_data = num_data / nproc ;
    else
        s_num_data = num_data / nproc + 1 ;

    dmin = myid*s_num_data ;
    if (dmin > num_data-1)
        return ;

    dmax = dmin + s_num_data ;
    if (dmax > num_data)
        dmax = num_data ;

    for (d = dmin ; d < dmax ; d++){

        if (myid == 0){
            if ((d - dmin) % 1000 == 0)
                printf("d = %d\n", d - dmin) ;
        }

        read_cbf(det, d) ;

        for (t = 0 ; t < num_pix ; t++){
            pid = rec2pix_map[t] ;
            if (det[pid] < 0 || det[pid] > hot_pix_thres)
                photon_count = -1 ;
            else
                photon_count = det[pid] ;
            data_frame[t] = photon_count ;
        }

        memset(radial_val, 0, sizeof(double)*qlen) ;
        memset(radial_weight, 0, sizeof(double)*qlen) ;
        memset(radial_ct, 0, sizeof(int)*qlen) ;

        for (t = 0 ; t < num_pix ; t++){
            qid = qid_map[t] ;
            if (qid < 0)
                continue ;
            if (data_frame[t] < 0)
                continue ;
            radial_val[qid] += data_frame[t] ;
            radial_weight[qid] += pix[t][3] ;
            radial_ct[qid] += 1 ;
        }

        for (t = 0 ; t < qlen ; t++){
            if (radial_ct[t] > 0)
                ave_bg[t] = radial_val[t] / radial_weight[t] ;
            else
                ave_bg[t] = -1 ;
        }

        for (i = 0 ; i < num_iter ; i++){
            memset(radial_val, 0, sizeof(double)*qlen) ;
            memset(radial_weight, 0, sizeof(double)*qlen) ;
            memset(radial_ct, 0, sizeof(int)*qlen) ;
            for (t = 0 ; t < num_pix ; t++){
                qid = qid_map[t] ;
                if (qid < 0)
                    continue ;
                photon_count = data_frame[t] ;
                if (photon_count < 0 || ave_bg[qid] < 0)
                    continue ;
                avg_val = ave_bg[qid]*pix[t][3] ;
                if (avg_val < dv)
                    avg_val = dv ;
                idx = ((int) round(avg_val*inv_dv - 1.)) ;
                if (idx >= len_val)
                    thres = ad_hoc_thres(avg_val) ;
                else
                    thres = count_thres[idx] ;
                if (photon_count > thres)
                    continue ;
                radial_val[qid] += photon_count ;
                radial_weight[qid] += pix[t][3] ;
                radial_ct[qid] += 1 ;
            }

            for (t = 0 ; t < qlen ; t++){
                if (radial_ct[t] > 0)
                    ave_bg[t] = radial_val[t] / radial_weight[t] ;
                else
                    ave_bg[t] = -1 ;
            }
        }

        outlier_ct = 0 ;
        for (t = 0 ; t < num_pix ; t++){
            qid = qid_map[t] ;
            if (qid < 0)
                continue ;
            photon_count = data_frame[t] ;
            if (photon_count < 0 || ave_bg[qid] < 0)
                continue ;
            avg_val = ave_bg[qid]*pix[t][3] ;
            if (avg_val < dv)
                avg_val = dv ;
            idx = ((int) round(avg_val*inv_dv - 1.)) ;
            if (idx >= len_val)
                thres = ad_hoc_thres(avg_val) ;
            else
                thres = count_thres[idx] ;
            if (photon_count > thres){
                peak_id[outlier_ct] = rec2pix_map[t] ;
                peak_val[outlier_ct] = photon_count ;
                peak_map[peak_id[outlier_ct]] = outlier_ct ;
                peak_label[outlier_ct] = -1 ;
                outlier_ct += 1 ;
            }
        }

        mask_rings() ;
        peak_finder(d) ;

        fp = fopen(radialfiles[d].name, "wb") ;
        fwrite(ave_bg, sizeof(double), qlen, fp) ;
        fclose(fp) ;

        fp = fopen(outlierfiles[d].name, "w") ;
        for (t = 0 ; t < outlier_ct ; t++){
            pid = pix_map[peak_id[t]] ;
            qid = qid_map[pid] ;
            if (qid < 0 || ave_bg[qid] < 0)
                continue ;
            fprintf(fp, "%d %d\n", peak_id[t], peak_val[t]) ;
        }
        fclose(fp) ;
        
        /* reset peak_map */
        for (t = 0 ; t < outlier_ct ; t++)
            peak_map[peak_id[t]] = -1 ;
    }

    if (myid == 0){
        fp = fopen("data-info.dat", "w") ;
        fclose(fp) ;
    }
    MPI_Barrier(MPI_COMM_WORLD) ;

    for (i = 0 ; i < nproc ; i++){
        if (myid == i){
            fp = fopen("data-info.dat", "a") ;
            for (d = dmin ; d < dmax ; d++){
                fprintf(fp, "%d ", data_info[d].row_size) ;
                fprintf(fp, "%d ", data_info[d].col_size) ;
                fprintf(fp, "%.1f %.1f ", data_info[d].c_row, data_info[d].c_col) ;
                fprintf(fp, "%1.5e ", data_info[d].exposure) ;
                fprintf(fp, "%1.5e ", data_info[d].wl) ;
                fprintf(fp, "%1.5e\n", data_info[d].detd) ;
            }
            fclose(fp) ;
        }
        MPI_Barrier(MPI_COMM_WORLD) ;
    }

    free(data_frame) ;
}


void peak_finder( int d ){

    int i, j, t, patch_sz, num_peak, cur_idx ;
    int nn_pid[4], num_nn = 4, ct, num_patch, pid ;
    int qid, qid_min, qid_max, if_ring, ring_ct ;
    double val, val_sum, norm, inv_dq, weighted_rvec[3] ;
    FILE *fp ;

    inv_dq = 1./dq ;
    memset(radial_ct, 0, sizeof(int)*qlen) ;

    /* filter peaks based on number of connected pixels */
    num_patch = 0 ;
    num_peak = 0 ;
    ct = 0 ;
    for (i = 0 ; i < outlier_ct ; i++){
        if (peak_label[i] > -1)
            continue ;
        peak_label[i] = num_patch ;
        patch_sz = 0 ;
        queue[patch_sz] = i ;
        patch_sz += 1 ;
        cur_idx = 0 ;
        while (cur_idx < patch_sz){
            /* (x-1, y) */
            nn_pid[0] = peak_id[queue[cur_idx]] - num_col ;
            /* (x, y-1) */
            nn_pid[1] = peak_id[queue[cur_idx]] - 1 ;
            /* (x, y+1) */
            nn_pid[2] = peak_id[queue[cur_idx]] + 1 ;
            /* (x+1, y) */
            nn_pid[3] = peak_id[queue[cur_idx]] + num_col ;
            for (t = 0 ; t < num_nn ; t++){
                if (nn_pid[t] < 0 || nn_pid[t] > total_pix-1)
                    continue ;
                if (peak_map[nn_pid[t]] > -1){
                    if (peak_label[peak_map[nn_pid[t]]] > -1)
                        continue ;
                    queue[patch_sz] = peak_map[nn_pid[t]] ;
                    peak_label[queue[patch_sz]] = num_patch ;
                    patch_sz += 1 ;
                }
            }
            cur_idx += 1 ;
        }

        /* at least 2 connected pixels */
        if (patch_sz > 1){
            peak_list[ct] = patch_sz ;
            ct += 1 ;
            for (t = 0 ; t < patch_sz ; t++){
                peak_list[ct+2*t] = peak_id[queue[t]] ;
                peak_list[ct+2*t+1] = peak_val[queue[t]] ;
            }
            ct += 2*patch_sz ;
            num_peak += 1 ;
        }
        num_patch += 1 ;
    }

    ct = 0 ;
    for (i = 0 ; i < num_peak ; i++){
        patch_sz = peak_list[ct] ;
        ct += 1 ;
        val_sum = 0. ;
        for (j = 0 ; j < 3 ; j++)
            weighted_rvec[j] = 0. ;
        for (t = 0 ; t < patch_sz ; t++){
            pid = pix_map[peak_list[ct+2*t]] ;
            val = peak_list[ct+2*t+1] ;
            for (j = 0 ; j < 3 ; j++)
                weighted_rvec[j] += pix[pid][j]*val ;
            val_sum += val ;
        }
        
        norm = 0. ;
        for (j = 0 ; j < 3 ; j++){
            weighted_rvec[j] /= val_sum ;
            norm += weighted_rvec[j]*weighted_rvec[j] ;
        }
        qid = ((int) round(sqrt(norm)*inv_dq - 0.5)) ;
        if (qid >= 0 && qid < qlen)
            radial_ct[qid] += 1 ;

        ct += 2*patch_sz ;
    }

    ring_ct = 0 ;
    for (i = 0 ; i < qlen ; i++){
        if (radial_ct[i] > max_num_peak_on_ring){
            queue[ring_ct] = i ;
            ring_ct += 1 ;
        }
    }

    ct = 0 ;
    fp = fopen(peakfiles[d].name, "w") ;
    for (i = 0 ; i < num_peak ; i++){
        patch_sz = peak_list[ct] ;
        ct += 1 ;
        val_sum = 0. ;
        for (j = 0 ; j < 3 ; j++)
            weighted_rvec[j] = 0. ;
        for (t = 0 ; t < patch_sz ; t++){
            pid = pix_map[peak_list[ct+2*t]] ;
            val = peak_list[ct+2*t+1] ;
            for (j = 0 ; j < 3 ; j++)
                weighted_rvec[j] += pix[pid][j]*val ;
            val_sum += val ;
        }
        
        norm = 0. ;
        for (j = 0 ; j < 3 ; j++){
            weighted_rvec[j] /= val_sum ;
            norm += weighted_rvec[j]*weighted_rvec[j] ;
        }
        qid = ((int) round(sqrt(norm)*inv_dq - 0.5)) ;

        if_ring = 0 ;
        for (t = 0 ; t < ring_ct ; t++){
            if (qid == queue[t])
                if_ring = 1 ;
        }
        
        if (qid > -1 && qid < qlen){
            if (ave_bg[qid] < 0)
                if_ring = 1 ;
        }

        if (if_ring == 0){
            fprintf(fp, "%d\n", patch_sz) ;
            for (t = 0 ; t < patch_sz ; t++)
                fprintf(fp, "%d %d ", peak_list[ct+2*t], peak_list[ct+2*t+1]) ;
            fprintf(fp, "\n") ;
        }
        else{
            qid_min = qid ;
            qid_max = qid ;
            for (t = 0 ; t < patch_sz ; t++){
                pid = pix_map[peak_list[ct+2*t]] ;
                norm = 0. ;
                for (j = 0 ; j < 3 ; j++)
                    norm += pix[pid][j]*pix[pid][j] ;
                qid = ((int) round(sqrt(norm)*inv_dq - 0.5)) ;
                if (qid_min > qid)
                    qid_min = qid ;
                if (qid_max < qid)
                    qid_max = qid ;
            }
            qid_max += 1 ;
            if (qid_min < 0)
                qid_min = 0 ;
            if (qid_max > qlen-1)
                qid_max = qlen ;

            for (t = qid_min ; t < qid_max ; t++)
                ave_bg[t] = -1 ;
        }
        ct += 2*patch_sz ;
    }
    fclose(fp) ;
}


void mask_rings(){

    int i, j, jmin, jmax, idx, idx_start, idx_end ;
    int flag, ct, half_width = 3 ;
    double cur_val ;
    
    idx_start = ((int) ceil(qmin/dq)) ;
    idx_end = qlen ;

    /* find local maxima */
    ct = 0 ;
    for (i = idx_start ; i < idx_end ; i++){
        if (ave_bg[i] < 0)
            continue ;

        flag = 0 ;
        cur_val = ave_bg[i] ;
        for (j = 1 ; j < half_width ; j++){
            idx = i + j ;
            if (idx > idx_end-1 || ave_bg[idx] < 0)
                continue ;
            if (ave_bg[idx] < cur_val)
                cur_val = ave_bg[idx] ;
            else
                flag = 1 ;
        }

        if (flag == 1)
            continue ;
        
        flag = 0 ;
        cur_val = ave_bg[i] ;
        for (j = 1 ; j < half_width ; j++){
            idx = i - j ;
            if (idx < idx_start || ave_bg[idx] < 0)
                continue ;
            if (ave_bg[idx] < cur_val)
                cur_val = ave_bg[idx] ;
            else
                flag = 1 ;
        }

        if (flag == 1)
            continue ;
        
        jmin = i ;
        while (1){
            if (jmin == idx_start)
                break ;
            else{
                if (ave_bg[jmin-1] > ave_bg[jmin])
                    break ;
            }
            jmin -= 1 ;
        }

        jmax = i ;
        while (1){
            if (jmax == idx_end-1)
                break ;
            else{
                if (ave_bg[jmax+1] > ave_bg[jmax])
                    break ;
            }
            jmax += 1 ;
        }

        cur_val = (i-jmin)*ave_bg[jmax] + (jmax-i)*ave_bg[jmin] ;
        cur_val /= (jmax-jmin) ;

        if (ave_bg[i] > bg_thres*cur_val){
            peak_location[ct][0] = jmin ;
            peak_location[ct][1] = jmax ;
            ct += 1 ;
        }
    }

    for (i = 0 ; i < ct ; i++){
        jmin = peak_location[i][0] ;
        jmax = peak_location[i][1] ;
        for (j = jmin ; j < jmax + 1 ; j++)
            ave_bg[j] = -1 ;
    }
}


int setup(){

    FILE *fp ;
    char *token, line[256] ;
    int i, d, t, idx ;
    double qval, qval_max ;

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
            else if (strcmp(token, "qlen") == 0)
                qlen = atoi(strtok(NULL, " =")) ;
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
            else if (strcmp(token, "home_dir") == 0)
                strcpy(home_dir, strtok(NULL, " =\n")) ;
        }
        fclose(fp) ;
    }

    cbf_files = malloc(num_data * sizeof(filename)) ;
    fp = fopen("cbflist.dat", "r") ;
    for (d = 0 ; d < num_data ; d++)
        fscanf(fp, "%s", cbf_files[d].name) ;
    fclose(fp) ;

    outlierfiles = malloc(num_data * sizeof(filename)) ;
    fp = fopen("outlierlist.dat", "r") ;
    for (d = 0 ; d < num_data ; d++)
        fscanf(fp, "%s", outlierfiles[d].name) ;
    fclose(fp) ;

    radialfiles = malloc(num_data * sizeof(filename)) ;
    fp = fopen("radiallist.dat", "r") ;
    for (d = 0 ; d < num_data ; d++)
        fscanf(fp, "%s", radialfiles[d].name) ;
    fclose(fp) ;

    peakfiles = malloc(num_data * sizeof(filename)) ;
    fp = fopen("peaklist.dat", "r") ;
    for (d = 0 ; d < num_data ; d++)
        fscanf(fp, "%s", peakfiles[d].name) ;
    fclose(fp) ;

    data_info = malloc(num_data * sizeof(frame_info)) ;

    total_pix = num_row*num_col ;
    pix_map = malloc(total_pix * sizeof(int)) ;
    sprintf(pixmap_file, "%s/make-detector/pix-map.dat", home_dir) ;
    fp = fopen(pixmap_file, "r") ;
    for (t = 0 ; t < total_pix ; t++)
        fscanf(fp, "%d", &pix_map[t]) ;
    fclose(fp) ;

    peak_id = calloc(total_pix, sizeof(int)) ;
    peak_val = calloc(total_pix, sizeof(int)) ;
    peak_label = calloc(total_pix, sizeof(int)) ;
    peak_map = malloc(total_pix * sizeof(int)) ;
    queue = malloc(total_pix * sizeof(int)) ;
    peak_list = malloc(total_pix * sizeof(int)) ;
    for (i = 0 ; i < total_pix ; i++)
        peak_map[i] = -1 ;

    sprintf(rvec_file, "%s/make-detector/rec-vectors.bin", home_dir) ;
    fp = fopen(rvec_file, "rb") ;
    fread(&qmax, sizeof(int), 1, fp) ;
    fread(&qmin, sizeof(int), 1, fp) ;
    fread(&num_pix, sizeof(int), 1, fp) ;
    qval_max = qmax ;
    dq = qval_max/qlen ;

    pix = malloc(num_pix * sizeof(*pix)) ;
    qid_map = malloc(num_pix * sizeof(int)) ;
    for (t = 0 ; t < num_pix ; t++){
        qval = 0. ;
        for (i = 0 ; i < 3 ; i++){
            fread(&pix[t][i], sizeof(double), 1, fp) ;
            qval += pix[t][i]*pix[t][i] ;
        }
        fread(&pix[t][3], sizeof(double), 1, fp) ;
        qval = sqrt(qval) ;
        idx = ((int) round(qval/dq - 0.5)) ;
        if (idx < 0 || idx > qlen-1)
            qid_map[t] = -1 ;
        else
            qid_map[t] = idx ;
    }
    fclose(fp) ;

    rec2pix_map = malloc(num_pix * sizeof(int)) ;
    sprintf(rec2pixmap_file, "%s/make-detector/rec2pix-map.dat", home_dir) ;
    fp = fopen(rec2pixmap_file, "r") ;
    for (t = 0 ; t < num_pix ; t++)
        fscanf(fp, "%d", &rec2pix_map[t]) ;
    fclose(fp) ;

    det = malloc(total_pix * sizeof(int)) ;
    radial_ct = malloc(qlen * sizeof(int)) ;
    peak_location = malloc(qlen * sizeof(*peak_location)) ;
    radial_val = malloc(qlen * sizeof(double)) ;
    radial_weight = malloc(qlen * sizeof(double)) ;
    ave_bg = malloc(qlen * sizeof(double)) ;
    
    make_count_thres() ;
    return 1 ;
}


void free_mem(){

    free(cbf_files) ;
    free(outlierfiles) ;
    free(radialfiles) ;
    free(peakfiles) ;
    free(data_info) ;
    free(pix_map) ;
    free(pix) ;
    free(peak_id) ;
    free(peak_val) ;
    free(peak_label) ;
    free(peak_map) ;
    free(peak_list) ;
    free(queue) ;
    free(qid_map) ;
    free(rec2pix_map) ;
    free(det) ;
    free(radial_ct) ;
    free(peak_location) ;
    free(radial_val) ;
    free(radial_weight) ;
    free(ave_bg) ;
    free(count_thres) ;
}


void make_count_thres(){

    int i, kmax ;
    double *val, *log_val, *cdf, *log_fact, logP ;

    len_val = ((int) round(vmax/dv)) ;
    val = malloc(len_val * sizeof(double)) ;
    log_val = malloc(len_val * sizeof(double)) ;
    cdf = malloc(len_val * sizeof(double)) ;
    count_thres = malloc(len_val * sizeof(int)) ;

    kmax = ((int) ceil(vmax + 10*sqrt(vmax))) ;
    log_fact = malloc((kmax+1) * sizeof(double)) ;
    log_fact[0] = 0. ;
    for (i = 1 ; i < kmax+1 ; i++)
        log_fact[i] = log_fact[i-1] + log(i) ;
    
    for (i = 0 ; i < len_val ; i++){
        count_thres[i] = 0 ;
        val[i] = (i+1)*dv ;
        log_val[i] = log(val[i]) ;
        cdf[i] = exp(-val[i]) ;
    }

    /* count_thres is exclusive:
     * the probability to obtain count > count_thres[i] given 
     * an averge val[i] is smaller than 1-cdf_thres */
    for (i = 0 ; i < len_val ; i++){
        while (cdf[i] < cdf_thres){
            count_thres[i] += 1 ;
            logP = -val[i] + count_thres[i]*log_val[i] - log_fact[count_thres[i]] ;
            cdf[i] += exp(logP) ;
        }
    }

    free(val) ;
    free(log_val) ;
    free(cdf) ;
    free(log_fact) ;
}


int ad_hoc_thres( double ave_val ){

    int count_thres ;
    double cdf, log_fact, log_val, logP ;

    count_thres = 0 ;
    log_fact = 0. ;
    log_val = log(ave_val) ;
    cdf = exp(-ave_val) ;
    while (cdf < cdf_thres){
        count_thres += 1 ;
        log_fact += log(count_thres) ;
        logP = -ave_val + count_thres*log_val - log_fact ;
        cdf += exp(logP) ;
    }

    return count_thres ;
}
