/*
make pseudo-powder patterns from the Bragg peak candidates

compile:
gcc make-powder.c -O3 -lm -o powder.o

usage:
./powder.o [path of the config_file] > run.log &

needs:
config_file, peakfiles.dat

makes:
1d-pseudo-powder.dat, 2d-pseudo-powder.dat
frame-peak-count.dat, peak-sz-count.dat

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

typedef struct{
   char name[256] ;
} filename ;

const int nbin = 3000 ;
const double sampling_sz = 2. ;
char config_file[256], home_dir[256] ;
int num_row, num_col, num_data, total_pix ;
int min_patch_sz, max_patch_sz, min_num_peak, max_num_peak ;
double wl_D, (*pix)[3] ;
filename *peakfiles ;

void setup() ;
void free_mem() ;
void make_powder() ;

int main(int argc, char *argv[]){
    
    time_t t1, t2 ;
    time(&t1) ;

    if (argc < 2){
        printf("The config_file is missing!!\n") ;
        exit(1) ;
    }
    else
        sprintf(config_file, "%s", argv[1]) ;

    setup() ;
    make_powder() ;

    time(&t2) ;
    printf ("elapsed time: %.f sec\n", difftime(t2, t1));
    
    free_mem() ;
    return 0 ;
}


void make_powder(){

    int i, d, r, t, num_pix, pix_id, pix_val, num_peak, idx ;
    int *frame_peak_count, count_relevant, *peak_sz_hist, *powder2d ;
    double (*radial_hist)[2], *reflection, (*peak_rvec)[3] ;
    double weighted_rvec[3], diff_rvec[3], val_sum, dist, rescale ;
    FILE *fp ;

    peak_rvec = malloc(total_pix * sizeof(*peak_rvec)) ;
    reflection = malloc(total_pix * sizeof(double)) ;
    frame_peak_count = calloc(num_data, sizeof(int)) ;
    peak_sz_hist = calloc(nbin, sizeof(int)) ;
    powder2d = calloc(total_pix, sizeof(int)) ;
    radial_hist = malloc(nbin * sizeof(*radial_hist)) ;
    for (i = 0 ; i < nbin ; i++){
        radial_hist[i][0] = 0. ;
        radial_hist[i][1] = 0. ;
    }

    count_relevant = 0 ;
    for (d = 0 ; d < num_data ; d++){
        fp = fopen(peakfiles[d].name, "r") ;
        num_peak = 0 ;
        while (1 == fscanf(fp, "%d", &num_pix)){
            peak_sz_hist[num_pix] += 1 ;
            if (num_pix < min_patch_sz || num_pix > max_patch_sz){
                for (t = 0 ; t < num_pix ; t++)
                    fscanf(fp, "%d %d", &pix_id, &pix_val) ;
            }
            else{
                val_sum = 0. ;
                reflection[num_peak] = 0. ;
                for (i = 0 ; i < 3 ; i++)
                    weighted_rvec[i] = 0. ;
                for (t = 0 ; t < num_pix ; t++){
                    fscanf(fp, "%d %d", &pix_id, &pix_val) ;
                    if (powder2d[pix_id] < pix_val)
                        powder2d[pix_id] = pix_val ;
                    reflection[num_peak] += pix_val ;
                    for (i = 0 ; i < 3 ; i++)
                        weighted_rvec[i] += pix[pix_id][i]*pix_val ;
                    val_sum += pix_val ;
                }
                for (i = 0 ; i < 3 ; i++)
                    peak_rvec[num_peak][i] = weighted_rvec[i]/val_sum ;
                num_peak += 1 ;
            }
        }
        fclose(fp) ;

        if (num_peak < min_num_peak || num_peak > max_num_peak)
            continue ;
        
        frame_peak_count[d] = num_peak ;
        count_relevant += 1 ;
        for (r = 1 ; r < num_peak ; r++){
            for (t = 0 ; t < r ; t++){
                dist = 0. ;
                for (i = 0 ; i < 3 ; i++){
                    diff_rvec[i] = peak_rvec[r][i] - peak_rvec[t][i] ;
                    dist += diff_rvec[i]*diff_rvec[i] ;
                }
                dist = sqrt(dist)*sampling_sz ;
                idx = ((int) round(dist)) ;
                if (idx < nbin)
                    radial_hist[idx][0] += 1 ;
            }
        }
        
        for (r = 0 ; r < num_peak ; r++){
            dist = 0. ;
            for (i = 0 ; i < 3 ; i++){
                diff_rvec[i] = peak_rvec[r][i] ;
                dist += diff_rvec[i]*diff_rvec[i] ;
            }
            dist = sqrt(dist)*sampling_sz ;
            idx = ((int) round(dist)) ;
            if (idx < nbin)
                radial_hist[idx][1] += 1 ;
        }
    }

    printf("number of relevant frames = %d\n", count_relevant) ;
    rescale = 1./wl_D/sampling_sz ;
    fp = fopen("1d-pseudo-powder.dat", "w") ;
    for (i = 0 ; i < nbin ; i++)
        fprintf(fp, "%1.5e %1.5e %1.5e\n", i*rescale, radial_hist[i][0], radial_hist[i][1]) ;
    fclose(fp) ;

    fp = fopen("frame-peak-count.dat", "w") ;
    for (d = 0 ; d < num_data ; d++)
        fprintf(fp, "%d\n", frame_peak_count[d]) ;
    fclose(fp) ;

    fp = fopen("peak-sz-count.dat", "w") ;
    for (i = min_patch_sz ; i < nbin ; i++)
        fprintf(fp, "%d %d\n", i, peak_sz_hist[i]) ;
    fclose(fp) ;

    fp = fopen("2d-pseudo-powder.dat", "w") ;
    for (i = 0 ; i < total_pix ; i++)
        fprintf(fp, "%d\n", powder2d[i]) ;
    fclose(fp) ;

    free(radial_hist) ;
    free(peak_rvec) ;
    free(reflection) ;
    free(frame_peak_count) ;
    free(peak_sz_hist) ;
    free(powder2d) ;
}


void setup(){

    int d, t ;
    char *token, line[256], peakfilelist[256] ;
    double wl, px, detd, cx, cy, beam_vec[3] ;
    double norm, x, y, sx, sy, sz, D ;
    FILE *fp ;

    fp = fopen(config_file, "r") ;
    if (!fp){
        printf("The config_file %s is not found!!\n", config_file) ;
        exit(1) ;
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
            else if (strcmp(token, "home_dir") == 0)
                strcpy(home_dir, strtok(NULL, " =\n")) ;
        }
        fclose(fp) ;
    }

    norm = sqrt(beam_vec[0]*beam_vec[0] + beam_vec[1]*beam_vec[1] + beam_vec[2]*beam_vec[2]) ;
    D = (detd/px) * norm/fabs(beam_vec[2]) ;
    wl_D = wl*D ;

	beam_vec[0] *= D/norm ;
	beam_vec[1] *= D/norm ;
	beam_vec[2] *= D/norm ;

    total_pix = num_row*num_col ;
    pix = malloc(total_pix * sizeof(*pix)) ;
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
        pix[t][0] = sx - beam_vec[0] ;
        pix[t][1] = sy - beam_vec[1] ;
        pix[t][2] = sz - beam_vec[2] ;
	}

    sprintf(peakfilelist, "%s/make-background/peaklist.dat", home_dir) ;
    peakfiles = malloc(num_data * sizeof(filename)) ;
    fp = fopen(peakfilelist, "r") ;
    for (d = 0 ; d < num_data ; d++)
        fscanf(fp, "%s", peakfiles[d].name) ;
    fclose(fp) ;
}


void free_mem(){
    free(pix) ;
    free(peakfiles) ;
}
