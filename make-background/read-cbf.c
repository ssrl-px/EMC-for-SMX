#include "make-bg.h"

int read_cbf( int *det, int d ){
    
    char *buffer, *last, *pch, *content_md5, *conversions ;
    char *ascii_header, *binary_header, *binary_data ;
    char BINARAY_SECTION[] = "--CIF-BINARY-FORMAT-SECTION--" ;
    /* the length of STARTER has to be 4 to exclude the null terminator "\0" */
    char STARTER[4] = "\x0c\x1a\x04\xd5" ;
    FILE *fin ;
    int i, j, num_iter, binary_size, memstride = 256 ;
    int len_ascii_header, len_binary_header, len_binary_data ;

    /* header size <= 8 kB */
    num_iter = 8192 / memstride ;
    buffer = malloc(sizeof(char) * memstride) ;
    last = malloc(sizeof(char) * (sizeof(BINARAY_SECTION) + memstride)) ;
    ascii_header = malloc(sizeof(char) * memstride*num_iter) ;

    fin = fopen(cbf_files[d].name, "rb") ;
    for (j = 0 ; j < num_iter ; j++){
        
        fread(buffer, sizeof(char), memstride, fin) ;
        for (i = 0 ; i < memstride ; i++){
            last[i+sizeof(BINARAY_SECTION)] = buffer[i] ;
            ascii_header[i+j*memstride] = buffer[i] ;
        }

        if (j == 0)
            pch = strstr(buffer, BINARAY_SECTION) ;
        else
            pch = strstr(last, BINARAY_SECTION) ;
        
        if (pch != NULL){
            
            if (j == 0)
                len_ascii_header = pch - buffer ;
            else
                len_ascii_header = j*memstride - sizeof(BINARAY_SECTION) + pch - last ;
            
            len_binary_header = (j+1)*memstride-len_ascii_header ;
            binary_header = malloc(sizeof(char) * len_binary_header) ;
            for (i = 0 ; i < len_binary_header ; i++)
                binary_header[i] = ascii_header[i+len_ascii_header] ;
            
            break ;
        }

        for (i = 0 ; i < sizeof(BINARAY_SECTION) ; i++)
            last[i] = buffer[i+memstride-sizeof(BINARAY_SECTION)] ;
    }

    ascii_header = realloc(ascii_header, sizeof(char)*len_ascii_header) ;

    /* look for the start of data */
    pch = strstr(binary_header, STARTER) ;
    while (pch == NULL){
        binary_header = realloc(binary_header, sizeof(char)*(len_binary_header+memstride)) ;
        fread(binary_header+len_binary_header, sizeof(char), memstride, fin) ;
        len_binary_header += memstride ;
        pch = strstr(binary_header, STARTER) ;
    }
    
    /* binary_data contains the special sequence STARTER */
    len_binary_data = len_binary_header - (pch - binary_header) ;
    binary_data = malloc(sizeof(char) * len_binary_data) ;
    len_binary_header -= len_binary_data ;
    for (i = 0 ; i < len_binary_data ; i++)
        binary_data[i] = binary_header[i+len_binary_header] ;

    binary_size = parse_header(ascii_header, binary_header, conversions, content_md5, d) ;

    /* readin the rest of binary data, with size "binary_size" bytes */
    buffer = realloc(buffer, sizeof(char) * (len_binary_data - sizeof(STARTER))) ;
    for (i = 0 ; i < len_binary_data - sizeof(STARTER) ; i++)
        buffer[i] = binary_data[i+sizeof(STARTER)] ;

    binary_data = realloc(binary_data, sizeof(char) * binary_size) ;
    for (i = 0 ; i < len_binary_data - sizeof(STARTER) ; i++)
        binary_data[i] = buffer[i] ;

    len_binary_data -= sizeof(STARTER) ;
    fread(binary_data + len_binary_data, sizeof(char), binary_size - len_binary_data, fin) ;

    md5sum(content_md5, binary_data, binary_size) ;
    byte_offset(det, binary_data, binary_size) ;

    fclose(fin) ;
    free(last) ;
    free(buffer) ;
    free(ascii_header) ;
    free(binary_header) ;
    free(binary_data) ;

    return 1 ;
}


int parse_header( char *ascii_header, char *binary_header, char *conversions, char *content_md5, int d ){

    char *token ;
    int row_size, col_size, binary_size ;
    double exposure, wl, detd, c_row, c_col ;

    token = strtok(binary_header, " :=\r\n\"") ;
    while (token != NULL){

        if (strcmp(token, "X-Binary-Size-Fastest-Dimension") == 0)
            col_size = atoi(strtok(NULL, " :\r\n\"")) ;
        if (strcmp(token, "X-Binary-Size-Second-Dimension") == 0)
            row_size = atoi(strtok(NULL, " :\r\n\"")) ;
        if (strcmp(token, "X-Binary-Size") == 0)
            binary_size = atoi(strtok(NULL, " :\r\n\"")) ;
        if (strcmp(token, "conversions") == 0)
            conversions = strtok(NULL, " :\r\n\"") ;
        if (strcmp(token, "Content-MD5") == 0)
            content_md5 = strtok(NULL, " :\r\n\"") ;
        
        if (strcmp(token, "X-Binary-Element-Type") == 0){
            token = strtok(NULL, "\"") ;
            token = strtok(NULL, "\"") ;
            if (strcmp(token, "signed 32-bit integer") != 0){
                printf("data type is not int32: %s!!\n", token) ;
                exit(1) ;
            }
        }   

        token = strtok(NULL, " :=\r\n\"") ;
    }

    token = strtok(ascii_header, " :=(),\r\n");
    while (token != NULL){

        if (strcmp(token, "Exposure_time") == 0)
            exposure = atof(strtok(NULL, " #:=(),\r\n")) ;
        if (strcmp(token, "Wavelength") == 0)
            wl = atof(strtok(NULL, " #:=(),\r\n")) ;
        if (strcmp(token, "Detector_distance") == 0)
            detd = atof(strtok(NULL, " #:=(),\r\n")) ;
        if (strcmp(token, "Beam_xy") == 0){
            c_col = atof(strtok(NULL, " #:=(),\r\n")) ;
            c_row = atof(strtok(NULL, " #:=(),\r\n")) ;
        }

        token = strtok(NULL, " #:=(),\r\n");
    }

    data_info[d].row_size = row_size ;
    data_info[d].col_size = col_size ;
    data_info[d].exposure = exposure ;
    data_info[d].wl = wl ;
    data_info[d].detd = detd ;
    data_info[d].c_row = c_row ;
    data_info[d].c_col = c_col ;

    return binary_size ;
}


void byte_offset( int *data_frame, char *binary_data, int binary_size ){

    int base = 0, mem_ct = 0, idx = 0, delta ;
    short delta16 ;

    while (mem_ct < binary_size){
        
        delta = binary_data[mem_ct] ;
        if (delta > -128){
            base += delta ;
            data_frame[idx] = base ;
            mem_ct += 1 ;
            idx += 1 ;
            continue ;
        }
        
        mem_ct += 1 ;
        delta16 = (binary_data[mem_ct] & 0x00FF) | ((binary_data[mem_ct+1] << 8) & 0xFF00) ;
        if (delta16 > -32768){
            base += delta16 ;
            data_frame[idx] = base ;
            mem_ct += 2 ;
            idx += 1 ;
            continue ;
        }
        
        mem_ct += 2 ;
        delta = (binary_data[mem_ct] & 0x000000FF) | ((binary_data[mem_ct+1] << 8) & 0x0000FF00) | \
                ((binary_data[mem_ct+2] << 16) & 0x00FF0000) | ((binary_data[mem_ct+3] << 24) & 0xFF000000) ;

        mem_ct += 4 ;
        base += delta ;
        data_frame[idx] = base ;
        idx += 1 ;
    }
}


void md5sum( char *content_md5, char *binary_data, int binary_size ){
    
    unsigned char digest[MD5_DIGEST_LENGTH] ;
    char *b64code ;

    /* calculate the MD5 hash of binary data for check */
    MD5_CTX mdContext ;
    MD5_Init(&mdContext) ;
    MD5_Update(&mdContext, binary_data, binary_size) ;
    MD5_Final(digest, &mdContext) ;
}
