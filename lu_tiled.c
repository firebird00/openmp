/***** Tiled LU decomposition*****/

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>

void mm(double ** a1,int x_1,int y_1, double ** a2, int x_2, int y_2 , double ** res, int x_3, int y_3, double ** tmp_res, int M, int N, int P, int sign);
void mm_lower(double ** a1,int x_1,int y_1, double ** a2, int x_2, int y_2 , double ** res, int x_3, int y_3, int M, int N, int P);
void mm_upper(double ** a1,int x_1,int y_1, double ** a2, int x_2, int y_2 , double ** res, int x_3, int y_3, int M, int N, int P);
void rectrtri_lower(double ** a, int x_s,int y_s, int N, int M, double ** tmp_res);
void rectrtri_upper(double ** a, int x_s, int y_s, int N, int M, double ** tmp_res);
double ** allocate(int X, int Y);
int min(int x, int y) {return x<y?x:y; }
double ** get_inv_l(double ** a, int xs, int ys, int X, int Y);
double ** get_inv_u(double **a, int xs, int ys, int X, int Y);
void input(double ** a, int X, int Y);
void lu(double ** a, int range, int B);
void lu_kernel(double ** a, int xs, int ys, int X, int Y);
void print(double ** a, int N);
void mm_update(double ** a1, int x_1, int y_1, double ** a2, int x_2, int y_2, double ** res, int x_3, int y_3, int M, int N, int P);
double ** up_res, ** low_res;

int main(int argc, char *argv[])
{
	double ** A;
	int N,B,range;
	struct timeval ts, tf;
	double time;
	
	if (argc<5) {
		printf("Usage: ./[executable] [grid_size] [block_size]   [threads_num]   [print_flag]\n");
		exit(1);
	}
	int print_flag = atoi(argv[4]);
	int threads_num = atoi(argv[3]);
	N=atoi(argv[1]);
	B=atoi(argv[2]);

	if (N%B!=0 || B==1) {
		printf("Grid must be multiple of block and greater than 1\n");
		exit(1);
	}
	
	A=allocate(N,N);
	input(A,N,N);
	range=N/B;

	low_res=allocate(B,(range-1)*B);
	up_res=allocate((range-1)*B,B);



	gettimeofday(&ts,NULL);

	omp_set_num_threads(threads_num);
/*    int procs = omp_get_num_procs();
    int nthreads = omp_get_num_threads();
    int maxt = omp_get_max_threads();
    int inpar = omp_in_parallel();
    int dynamic = omp_get_dynamic();
    int nested = omp_get_nested();
*/    /* Print environment information */
/*    printf("Number of processors = %d\n", procs);
    printf("Number of threads = %d\n", nthreads);
    printf("Max threads = %d\n", maxt);
    printf("In parallel? = %d\n", inpar);
    printf("Dynamic threads enabled? = %d\n", dynamic);
    printf("Nested parallelism supported? = %d\n", nested);     */
#pragma omp parallel
{
   #pragma omp single
	lu(A,range,B);
}
	gettimeofday(&tf,NULL);
	time=(tf.tv_sec-ts.tv_sec)+(tf.tv_usec-ts.tv_usec)*0.000001;
	printf("Tiled\t%d\t%lf\t%d\t\n", N,time,B);
	if(print_flag)
		print(A,N);
	return 0;
}

/**** Tiled LU decomposition *****/
void lu(double **a, int range, int B)
{
	int i,j,k;
	double ** l_inv, ** u_inv;
	for (k=0;k<range-1;k++) {

		/****Compute LU decomposition on upper left tile*****/
		lu_kernel(a,k*B,k*B,B,B);

		/****Compute inverted L and U matrices of upper left tile*****/
   #pragma omp task shared(l_inv)
   {//printf("l_inv thread: %d\n",omp_get_thread_num());
		l_inv=get_inv_l(a,k*B,k*B,B,B);
   }
   #pragma omp task shared(u_inv)
   {//printf("r_inv thread: %d\n",omp_get_thread_num());
		u_inv=get_inv_u(a,k*B,k*B,B,B);
   }
   #pragma omp taskwait
//printf("lu code: %d\n",omp_get_thread_num());            ola ta kanei ena thread afou sthn main exw grapsei single

		/*****Compute LU decomposition on upper horizontal frame and left vertical frame*****/
		for (i=k+1;i<range;i++) {
   #pragma omp task
			mm_lower(l_inv,0,0,a,k*B,i*B,a,k*B,i*B,B,B,B);
   #pragma omp task
			mm_upper(a,i*B,k*B,u_inv,0,0,a,i*B,k*B,B,B,B);
		}

   #pragma omp taskwait
		/*****Update trailing blocks*****/
		for (i=k+1;i<range;i++) {
			for (j=k+1;j<range;j++)
   #pragma omp task
				mm_update(a,i*B,k*B,a,k*B,j*B,a,i*B,j*B,B,B,B);
		} 

   #pragma omp taskwait
	}

	/***** Compute LU on final diagonal block *****/
	lu_kernel(a,(range-1)*B,(range-1)*B,B,B);

}

/***** Baseline LU Kernel *****/
void lu_kernel(double ** a, int xs,int ys, int X, int Y)
{
	int i,j,k;
	double l;
	for (k=0;k<min(X,Y);k++)
		for (i=k+1;i<X;i++)
		{
			l=a[i+xs][k+ys]=a[i+xs][k+ys]/a[k+xs][k+ys];
			for (j=k+1;j<Y;j++)
				a[i+xs][j+ys]-=l*a[k+xs][j+ys];
		}
			
}

/***** Computes the inverted L^(-1) matrix of upper diagonal block using rectrtri_lower() *****/
double ** get_inv_l(double ** a, int xs, int ys, int X, int Y)
{
	double ** l_inv=allocate(X,Y);
	double ** tmp_res=allocate(X,Y);
	int i,j;
	for (i=0;i<X;i++)
		for (j=0;j<Y;j++) {
			if (i>j)
				l_inv[i][j]=a[i+xs][j+ys];
			else if (i==j)
				l_inv[i][j]=1.0;
			else l_inv[i][j]=0;
		}
	rectrtri_lower(l_inv,0,0,X,Y,tmp_res);
	free(tmp_res);
	return l_inv;

}

/***** Computer the inverted U^(-1) matrix of upper diagonal block using rectrtri_upper() *****/
double ** get_inv_u(double **a, int xs, int ys, int X, int Y)
{
	double ** u_inv=allocate(X,Y);
	double ** tmp_res=allocate(X,Y);
	int i,j;
	for (i=0;i<X;i++)
		for (j=0;j<Y;j++) {
			if (i<=j)
				u_inv[i][j]=a[i+xs][j+ys];
			else u_inv[i][j]=0;
		}
	rectrtri_upper(u_inv,0,0,X,Y,tmp_res);
	free(tmp_res);
	return u_inv;
	
}


double ** allocate(int X,int Y)
{
	double ** array;
	array=(double**)calloc(X,sizeof(double*));
	int i;
	for (i=0;i<X;i++)
		array[i]=(double*)calloc(Y,sizeof(double));
	return array;
}

void input(double ** a, int X, int Y)
{
	int i,j;
	for (i=0;i<X;i++)
		for (j=0;j<Y;j++)
			a[i][j]=((double)rand()/1000000.0);
}

void print(double ** a, int N)
{
	int i,j;
	FILE *output_file = fopen("./out_tiled","w");
	for (i=0;i<N;i++) {
		for (j=0;j<N;j++)
			fprintf(output_file,"%.2lf ",a[i][j]);
		fprintf(output_file,"\n");	
	}
}

/***** Matrix multiplication of an upper triangular matrix with a full matrix *****/
void mm_upper(double ** a1, int x_1, int y_1, double ** a2, int x_2, int y_2, double ** res, int x_3, int y_3, int M, int N, int P)
{
	int i,j,k;
	double sum=0;
	for (i=0;i<M;i++)
		for (j=0;j<P;j++) {
			for (k=0;k<=j;k++)
				sum=sum+a1[i+x_1][k+y_1]*a2[k+x_2][j+y_2];
			up_res[i+x_3-M][j]=sum;
			sum=0;
		}
	for (i=0;i<M;i++)
		for (j=0;j<P;j++)
			res[i+x_3][j+y_3]=up_res[i+x_3-M][j];
}

/***** Matrix multiplication of a lower triangular matrix with a full matrix *****/
void mm_lower(double ** a1, int x_1, int y_1, double ** a2, int x_2, int y_2, double ** res, int x_3, int y_3, int M, int N, int P)
{
	int i,j,k;
	double sum=0;

	for (i=0;i<M;i++)
		for (j=0;j<P;j++) {
			for (k=0;k<=i;k++)
				sum=sum+a1[i+x_1][k+y_1]*a2[k+x_2][j+y_2];
			low_res[i][j+y_3-P]=sum;
			sum=0;
		}
	for (i=0;i<M;i++)
		for (j=0;j<P;j++)
			res[i+x_3][j+y_3]=low_res[i][j+y_3-P];
		
}

void mm_update(double ** a1, int x_1, int y_1, double ** a2, int x_2, int y_2, double ** res, int x_3, int y_3, int M, int N, int P)
{
	int i,j,k;
	double sum=0;
	
	for (i=0;i<M;i++)
		for (j=0;j<P;j++) {
			for (k=0;k<N;k++)
				sum+=a1[i+x_1][k+y_1]*a2[k+x_2][j+y_2];
			res[i+x_3][j+y_3]-=sum;
			sum=0;
		}
}

/***** Matrix multiplication: sign=0 --> res=a1*a2, sign=1 --> res=-a1*a2, op = 0 --> res=a1*a2, op = 1 --> res=res+a1*a2  *****/
void mm(double ** a1,int x_1,int y_1, double ** a2, int x_2, int y_2 , double ** res, int x_3, int y_3, double ** tmp_res, int M, int N, int P, int sign)
{
	int i,j,k;
	double sum=0;

	for (i=0;i<M;i++)
		for (j=0;j<P;j++)
		{
			for (k=0;k<N;k++)
				sum=sum+a1[i+x_1][k+y_1]*a2[k+x_2][j+y_2];
			tmp_res[i+x_3][j+y_3]=sum;			
			sum=0;
		}
	
	for (i=0;i<M;i++)
		for (j=0;j<P;j++)
		{
			for (k=0;k<N;k++)
				sum=sum-a1[i+x_1][k+y_1]*a2[k+x_2][j+y_2];
			tmp_res[i+x_3][j+y_3]=sum;			
			sum=0;
		}
	
	if (sign==0)
	{
		for (i=0;i<M;i++)
			for (j=0;j<P;j++)
				res[i+x_3][j+y_3]=tmp_res[i+x_3][j+y_3];
	}
	else if (sign==1)
	{
		for (i=0;i<M;i++)
			for (j=0;j<P;j++)
				res[i+x_3][j+y_3]=-tmp_res[i+x_3][j+y_3];
	}
}


/**** Recursive matrix inversion of a lower triangular matrix *****/
void rectrtri_lower(double ** a, int x_s,int y_s, int N, int M, double ** tmp_res)
{
	int n,m;
	if (N==1 && M==1)
		a[x_s][y_s]=1.0/a[x_s][y_s];
	else {
		n=N/2;
		m=M/2;
		rectrtri_lower(a,x_s,y_s,n,m, tmp_res);
		rectrtri_lower(a,x_s+n,y_s+m,N-n,M-m, tmp_res);
		mm(a,x_s+n,y_s,a,x_s,y_s,a,x_s+n,y_s,tmp_res,N-n,m,m,0);
		mm(a,x_s+n,y_s+m,a,x_s+n,y_s,a,x_s+n,y_s,tmp_res,N-n,M-m,m,1);	
	}		
}

/***** Recursive matrix inversion of an upper triangular matrix *****/
void rectrtri_upper(double ** a, int x_s,int y_s, int N, int M, double ** tmp_res)
{
	int n,m;
	if (N==1 && M==1)
		a[x_s][y_s]=1.0/a[x_s][y_s];
	else
	{
		n=N/2;
		m=M/2;

		rectrtri_upper(a,x_s,y_s,n,m, tmp_res);
		rectrtri_upper(a,x_s+n,y_s+m,N-n,M-m, tmp_res);
		mm(a,x_s,y_s+m,a,x_s+n,y_s+m,a,x_s,y_s+m,tmp_res,n,M-m,M-m,0);
		mm(a,x_s,y_s,a,x_s,y_s+m,a,x_s,y_s+m,tmp_res, n,m,M-m,1);	
	} 
	
}

