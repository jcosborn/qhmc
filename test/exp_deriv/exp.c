#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>

typedef complex double cdouble;

typedef struct {
  int m;
  int n;
  cdouble **elem;
} matrix;

void print_matrix(matrix *M){
  for ( int i = 0; i < M->m; i++ ) {
	for (int j = 0; j < M->n; j++ ) {
           if (cimag(M->elem[i][j])<0) 
		printf("%f-%fi ", creal(M->elem[i][j]),fabs(cimag(M->elem[i][j])));
	   else
                printf("%f+%fi ", creal(M->elem[i][j]),cimag(M->elem[i][j]));

	}
        printf("\n");
  }

  printf("\n");
}

void zero_matrix(matrix *M){
   for ( int i = 0; i < M->m; i++ ) {
	for ( int j = 0; j < M->n; j++) {
		M->elem[i][j] = 0.0 + 0.0*I;
	}
   }
}

void identity_matrix(matrix *M){
   for ( int i = 0; i < M->m; i++ ) {
        for ( int j = 0; j < M->n; j++) {
                if ( i == j ) M->elem[i][j] = 1.0 + 0.0*I;
		else M->elem[i][j] = 0.0+0.0*I;
        }
   }
}

void random_matrix(matrix *M) {
   srand(415);
   for ( int i = 0; i < M->m; i++ ) {
	for ( int j = 0 ; j < M->n; j++) {
		double re = ((double)rand())/RAND_MAX;
		double im = ((double)rand())/RAND_MAX;
                M->elem[i][j] = re+im*I;
	}
   }
}

//M^+
void adj(matrix *Ma, matrix *M) {
    if (Ma->m != M->n || Ma->n != M->m) {
	printf("Matrix dimensions incompatible\n");
	exit(-1);
    }

    for (int i = 0; i < M->n; ++i) {
	for ( int j = 0; j < M->m; ++j) {
		Ma->elem[i][j] = creal(M->elem[j][i]) - cimag(M->elem[j][i])*I;
	}
    } 

}

int main(int argc, char **argv) {
    printf("sizeof(cdouble)=%lu\n",sizeof(cdouble));
    matrix A;
    A.n = 3;
    A.m = 3;
    A.elem = malloc(sizeof(cdouble *) * A.n);

    for ( int i=0; i < A.m; i++) {
	A.elem[i] = malloc(sizeof(cdouble));
    }

    printf("%d x %d matrix created\n",A.n, A.m);
 
    zero_matrix(&A); 
    print_matrix(&A);

    identity_matrix(&A);
    print_matrix(&A);

    random_matrix(&A);
    print_matrix(&A);

    matrix Aa;
    Aa.n = 3;
    Aa.m = 3;
    Aa.elem = malloc(sizeof(cdouble *) * Aa.n);

    for ( int i=0; i < Aa.m; i++) {
        Aa.elem[i] = malloc(sizeof(cdouble));
    }

    zero_matrix(&Aa);

    adj(&Aa, &A);
    print_matrix(&Aa);

    for ( int i=0; i < A.m; ++i) {
 	//free(A.elem[i]);
    }
    printf("elem[] freed\n");

    //free(A.elem);

}
