#include "emc.h"

/* sort w.r.t values in descending order */
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

        if ( (i0 < iMiddle) && ((j0 >= iEnd) || (A[i0] > A[j0])) ){
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


/* sort w.r.t indices in ascending order */
void merge_sort_idx( double *A, double *B, int *A_idx, int *B_idx, int length ){
    TopDownSplitMergeIdx(A, B, A_idx, B_idx, 0, length) ;
}


void TopDownSplitMergeIdx( double *A, double *B, int *A_idx, int *B_idx, int iBegin, int iEnd ){

    // iBegin: inclusive, iEnd: exclusive
    if (iEnd - iBegin == 1)
        return ;

    int iMiddle = (iEnd + iBegin) / 2 ;
    TopDownSplitMergeIdx(A, B, A_idx, B_idx, iBegin, iMiddle) ;
    TopDownSplitMergeIdx(A, B, A_idx, B_idx, iMiddle, iEnd) ;
    TopDownMergeIdx(A, B, A_idx, B_idx, iBegin, iMiddle, iEnd) ;
    CopyArray(A, B, A_idx, B_idx, iBegin, iEnd) ;
}


void TopDownMergeIdx( double *A, double *B, int *A_idx, int *B_idx, int iBegin, int iMiddle, int iEnd ){

    int idx, i0 = iBegin, j0 = iMiddle ;

    for (idx = iBegin ; idx < iEnd ; idx++){

        if ( (i0 < iMiddle) && ((j0 >= iEnd) || (A_idx[i0] < A_idx[j0])) ){
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
