#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <complex.h>
#define block_size 3

/* Markers used to bound trace regions of interest */
volatile char MARKER_START, MARKER_END;

int main(int argc, char* argv[]) {

    /* Record marker addresses */
    FILE* marker_fp = fopen(".marker","w");
    assert(marker_fp);
    fprintf(marker_fp, "%llx\n%llx",
        (unsigned long long int) &MARKER_START,
        (unsigned long long int) &MARKER_END);
    fclose(marker_fp);

    if (argc != 3) {
        printf("Usage: ./matMul <matrix_a_file> <matrix_b_file>\n");
        exit(EXIT_FAILURE);
    }

    FILE* matrix_a_fp = fopen(argv[1], "r");
    if (!matrix_a_fp) {
        perror("fopen failed");
        return EXIT_FAILURE;
    }
    
    size_t n;
    int ret = fscanf(matrix_a_fp, "%ld\n", &n);
    assert (ret==1);
    complex* a = calloc( n*n, sizeof(complex) );
    for (size_t i = 0; i < n; i++) {
        for (size_t k = 0; k < n; k++) {
            double real, imag;
            ret = fscanf(matrix_a_fp, "(%lf%lfj) ", &real, &imag);
            assert (ret==2);
            a[i*n+k] = CMPLX(real, imag); // Directly store matrix A
        }
        ret = fscanf(matrix_a_fp, "\n");
    }
    fclose(matrix_a_fp);

    FILE* matrix_b_fp = fopen(argv[2], "r");
    if (!matrix_b_fp) {
        perror("fopen failed");
        return EXIT_FAILURE;
    }
    size_t m;
    ret = fscanf(matrix_b_fp, "%ld\n", &m);
    assert (ret==1);
    assert( n==m );
    complex* b = calloc( n*n, sizeof(complex) );
    for (size_t k = 0; k < n; k++) {
        for (size_t j = 0; j < n; j++) {
            double real, imag;
            ret = fscanf(matrix_b_fp, "(%lf%lfj) ", &real, &imag);
            assert (ret==2);
            b[k*n+j] = CMPLX(real, imag); // Directly store matrix B
        }
        ret = fscanf(matrix_b_fp, "\n");
    }
    fclose(matrix_b_fp);

    complex* c = calloc( n*n, sizeof(complex) );
    MARKER_START = 211;

    // //kij from matmul
    // for ( size_t k=0; k<n; k++ ) {
    //     for ( size_t i=0; i<n; i++ ) {
    //         complex r = a[i*n+k];
    //         for ( size_t j=0; j<n; j++ ) {
    //             c[i*n+j] += r * b[k*n+j];
    //         }
    //     }
    // }

    MARKER_END = 211;
    //Code is from textbook slides/lecture (version kij)
    for(int i = 0; i < n; i+=block_size){
        for(int j = 0; j < n; j+=block_size){
            for(int k = 0; k < n; k+= block_size){
                //B x B mini multiplication
                for(int k0 = k; (k0 < (k+block_size) && k0 < n); k0++){
                    for(int i0 = i; ((i0 < i+block_size) && i0 < n); i0++){
                        for(int j0 = j; ((j0 < j+block_size) && j0 < n); j0++){
                            c[i0*n+j0] += a[i0*n+k0] * b[k0*n+j0];
                        }
                    }
                }
            }
        }
    }

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            if (cimag(c[i*n+j]) < 0) {
                printf("(%.12lf%.12lfj) ", creal(c[i*n+j]), cimag(c[i*n+j]));
            } else {
                printf("(%.12lf+%.12lfj) ", creal(c[i*n+j]), cimag(c[i*n+j]));
            }
        }
        printf("\n");
    }

    free(c);
    free(b);
    free(a);
    exit(EXIT_SUCCESS);
}
