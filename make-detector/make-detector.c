/*
map detector pixels to the reciprocal space
the output q-vectors are sorted in ascending order w.r.t. their maginitudes

compile:
gcc make-detector.c -O3 -lm -o det.o

usage:
./det.o [path of the config_file] >> run.log &

needs:
config_file, mask.dat

makes:
pix-map.dat, rec2pix-map.dat, rec-vectors.bin
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* merge sort */
void merge_sort( double*, double*, int*, int*, int ) ;
void TopDownSplitMerge( double*, double*, int*, int*, int, int ) ;
void TopDownMerge( double*, double*, int*, int*, int, int, int ) ;
void CopyArray( double*, double*, int*, int*, int, int ) ;

int main(int argc, char *argv[]){

    FILE *fp ;
    char *token, line[256], config_file[256], mask_file[256], home_dir[256] ;
    int i, t, qmax, qmin, Npix, num_row, num_col, total_pixel ;
    int *mask, *pix_map, *rec2pix_map, *rec2pix_map_cpy ;
    double detd, wl, px, beam_vec[3], cx, cy, res_max, Rstop ;
    double x, y, sx, sy, sz, qx, qy, qz, norm, theta, D ;
    double polarization, solid_angle, scale_factor, max_scale_factor ;
    double *qval_arr, *qval_arr_cpy ;

    if (argc < 2){
        printf("The config_file is missing!!\n") ;
        exit(1) ;
    }
    else{
        sprintf(config_file, "%s", argv[1]) ;
        fp = fopen(config_file, "r") ;
        if (!fp){
            printf("The config_file %s is not found!!\n", config_file) ;
            exit(1) ;
        }
        
        while (fgets(line, 256, fp) != NULL){
            token = strtok(line, " =") ;
            if (token[0] == '#' || token[0] == '\n' || token[0] == '[')
                continue ;
            if (strcmp(token, "num_row") == 0)
                num_row = atoi(strtok(NULL, " =")) ;
            else if (strcmp(token, "num_col") == 0)
                num_col = atoi(strtok(NULL, " =")) ;
            else if (strcmp(token, "detd") == 0)
                detd = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "wl") == 0)
                wl = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "px") == 0)
                px = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "sx") == 0)
                beam_vec[0] = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "sy") == 0)
                beam_vec[1] = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "sz") == 0)
                beam_vec[2] = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "cx") == 0)
                cx = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "cy") == 0)
                cy = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "res_max") == 0)
                res_max = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "Rstop") == 0)
                Rstop = atof(strtok(NULL, " =")) ;
            else if (strcmp(token, "home_dir") == 0)
                strcpy(home_dir, strtok(NULL, " =\n")) ;
        }
        fclose(fp) ;
    }

    printf("run make-detector.c:\n\n") ;
    total_pixel = num_row*num_col ;
    mask = malloc(total_pixel * sizeof(int)) ;
    pix_map = malloc(total_pixel * sizeof(int)) ;
    rec2pix_map = malloc(total_pixel * sizeof(int)) ;
    rec2pix_map_cpy = malloc(total_pixel * sizeof(int)) ;
    qval_arr = malloc(total_pixel * sizeof(double)) ;
    qval_arr_cpy = malloc(total_pixel * sizeof(double)) ;
    for(int i=0; i< total_pixel; i++)
	    mask[i] = 0;

    /* masked pixels: 1, unmasked pixels: 0 */
    fp = fopen("MASK.bin", "rb") ;
    fread(mask,sizeof(mask[0]) ,total_pixel, fp );
    //for (int i=0; i < 20; i++)
    //    printf("MASK %d\n", mask[i]);
    fclose(fp);

    norm = 0. ;
    for (i = 0 ; i < 3 ; i++)
        norm += beam_vec[i]*beam_vec[i] ;
    norm = sqrt(norm) ;

    printf("det dist=%f meters, pxsize=%f meters, wavelength=%f Angstrom\n", detd, px, wl);
    D = (detd/px) * norm/fabs(beam_vec[2]) ;
    for (i = 0 ; i < 3 ; i++)
        beam_vec[i] *= D/norm ;
    
    qmax = ((int) ceil(wl*D/res_max)) ;
    theta = atan(Rstop/D) ;
    qmin = ((int) ceil(2*D*sin(theta/2))) ;

    /* count Npix */
    Npix = 0 ;
    max_scale_factor = 0. ;
    for (i = 0 ; i < total_pixel ; i++){
        x = i / num_col - cx ;
        y = i % num_col - cy ;
        sx = beam_vec[0] + x ;
        sy = beam_vec[1] + y ;
        sz = beam_vec[2] ;
        norm = sqrt(sx*sx + sy*sy + sz*sz) ;
        /* polarization is assumed to be in the y direction for synchrotrons */
        polarization = 1 - pow(sy/norm, 2) ;
        solid_angle = fabs(sz)/pow(norm, 3) ;
        scale_factor = polarization*solid_angle ;
        sx *= D/norm ;
        sy *= D/norm ;
        sz *= D/norm ;
        qx = sx - beam_vec[0] ;
        qy = sy - beam_vec[1] ;
        qz = sz - beam_vec[2] ;
        norm = sqrt(qx*qx + qy*qy + qz*qz) ;
        qval_arr[i] = norm ;
        if (mask[i] == 0){
            if (qmin <= norm && norm <= qmax){
                Npix += 1 ;
                if (scale_factor > max_scale_factor)
                    max_scale_factor = scale_factor ;
                rec2pix_map[i] = i ;
            }
            else
                rec2pix_map[i] = -1 ;
        }
        else
            rec2pix_map[i] = -1 ;
    }

    /* sort qval_arr in ascending order */
    merge_sort(qval_arr, qval_arr_cpy, rec2pix_map, rec2pix_map_cpy, total_pixel) ;

    int count = 0 ;
    FILE *outFp, *outFp2, *outFp3 ;

    outFp = fopen("rec-vectors.bin", "wb") ;
    outFp2 = fopen("pix-map.dat", "w") ;
    outFp3 = fopen("rec2pix-map.dat", "w") ;
    fwrite(&qmax, sizeof(int), 1, outFp) ;
    fwrite(&qmin, sizeof(int), 1, outFp) ;
    fwrite(&Npix, sizeof(int), 1, outFp) ;

    for (i = 0 ; i < total_pixel ; i++)
        pix_map[i] = -1 ;

    for (t = 0 ; t < total_pixel ; t++){
        if (rec2pix_map[t] < 0)
            continue ;
        x = rec2pix_map[t] / num_col - cx ;
        y = rec2pix_map[t] % num_col - cy ;
        sx = beam_vec[0] + x ;
        sy = beam_vec[1] + y ;
        sz = beam_vec[2] ;
        norm = sqrt(sx*sx + sy*sy + sz*sz) ;
        polarization = 1 - pow(sy/norm, 2) ;
        solid_angle = fabs(sz)/pow(norm, 3) ;
        scale_factor = polarization*solid_angle/max_scale_factor ;
        sx *= D/norm ;
        sy *= D/norm ;
        sz *= D/norm ;
        qx = sx - beam_vec[0] ;
        qy = sy - beam_vec[1] ;
        qz = sz - beam_vec[2] ;
        fwrite(&qx, sizeof(double), 1, outFp) ;
        fwrite(&qy, sizeof(double), 1, outFp) ;
        fwrite(&qz, sizeof(double), 1, outFp) ;
        fwrite(&scale_factor, sizeof(double), 1, outFp) ;
        fprintf(outFp3, "%d\n", rec2pix_map[t]) ;
        pix_map[rec2pix_map[t]] = count ;
        count += 1 ;
    }

    for (t = 0 ; t < total_pixel ; t++)
        fprintf(outFp2, "%d\n", pix_map[t]) ;

    // Make binary files, the ASCII files segfaulted when read back in during make-background portion
    FILE *outFp2_bin, *outFp3_bin;

    outFp2_bin = fopen("pix-map.bin", "w") ;
    fwrite(pix_map, sizeof(int), (size_t)total_pixel, outFp2_bin);
    fclose(outFp2_bin);

    outFp3_bin = fopen("rec2pix-map.bin", "w") ;
    fwrite(rec2pix_map, sizeof(int), (size_t)total_pixel, outFp3_bin);
    fclose(outFp3_bin);

    printf("cx = %.1f, cy = %.1f, Rstop = %.1f\n", cx, cy, Rstop) ;
    printf("qmax = %d, qmin = %d, num_pix = %d\n", qmax, qmin, Npix) ;

	free(mask) ;
	free(rec2pix_map) ;
	free(rec2pix_map_cpy) ;
    free(qval_arr) ;
    free(qval_arr_cpy) ;
	fclose(outFp) ;
    fclose(outFp2) ;
    fclose(outFp3) ;
	return 0 ;
}


/* sort values in ascending order */
void merge_sort( double *A, double *B, int *A_idx, int *B_idx, int length ){
    TopDownSplitMerge(A, B, A_idx, B_idx, 0, length) ;
}


void TopDownSplitMerge( double *A, double *B, int *A_idx, int *B_idx, int iBegin, int iEnd ){

    // iBegin: inclusive, iEnd: exclusive
    if (iEnd - iBegin == 1)
        return ;

    int iMiddle = (iEnd + iBegin) / 2 ;
    TopDownSplitMerge(A, B, A_idx, B_idx, iBegin, iMiddle) ;
    TopDownSplitMerge(A, B, A_idx, B_idx, iMiddle, iEnd) ;
    TopDownMerge(A, B, A_idx, B_idx, iBegin, iMiddle, iEnd) ;
    CopyArray(A, B, A_idx, B_idx, iBegin, iEnd) ;
}


void TopDownMerge( double *A, double *B, int *A_idx, int *B_idx, int iBegin, int iMiddle, int iEnd ){

    int idx, i0 = iBegin, j0 = iMiddle ;

    for (idx = iBegin ; idx < iEnd ; idx++){

        if ( (i0 < iMiddle) && ((j0 >= iEnd) || (A[i0] < A[j0])) ){
            B[idx] = A[i0] ;
            B_idx[idx] = A_idx[i0] ;
            i0 += 1 ;
        }
        else{
            B[idx] = A[j0] ;
            B_idx[idx] = A_idx[j0] ;
            j0 += 1 ;
        }
    }
}


void CopyArray( double *A, double *B, int *A_idx, int *B_idx, int iBegin, int iEnd ){

    int i ;
    for (i = iBegin ; i < iEnd ; i++){
        A[i] = B[i] ;
        A_idx[i] = B_idx[i] ;
    }
}
