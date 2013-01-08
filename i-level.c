#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <omp.h>

void printtable(float **A, int N);
void compute(float **A, int N);

int main (int argc, char *argv[] ){
	double **A;
	int N,count;
	int k,i,j;
	double l;
	struct timeval t_start,t_finish;
	double total_time;

	N = atoi( argv[1] );

	A=(double**)malloc(N*sizeof(double*));

	for (count=0; count<N; count++)
		A[count]=(double*)malloc(N*sizeof(double));

	srand( time(NULL) );

	for(count=0; count<N;count++){
		for(i=0;i<N;i++)
			A[count][i] = rand()%40 + 1;
	}
	//Start measure
	gettimeofday(&t_start,NULL);

	for (k = 0; k < N - 1; k++){
	#pragma omp parallel for private(i,l)
		for (i = k+1; i < N; i++){
			l = A[i][k] / A[k][k];
			#pragma omp parallel for private(j,l)
			for (j = k + 1; j < N; j++) {
				A[i][j] = A[i][j] -l*A[k][j];
			}
		}
	}
	//Time measure
	gettimeofday(&t_finish,NULL);
	total_time=(t_finish.tv_sec-t_start.tv_sec)+(t_finish.tv_usec-t_start.tv_usec)*0.000001;
	printf("Time: %lf sec\n",total_time);

	return 0;
}


void printtable(float **A, int N){
	int i,j;
	for (i=0; i<N; i++){
		for(j=0;j<N;j++){
			printf("%.2f ",A[i][j]);
		}
		printf("\n");
	}
}