#include <stdio.h>
#include <stdlib.h>

int main(int argc,char**argv){
	if(argc!=5){
		printf("\nusage: ./validate   file1  file2  N   error\n");
		return 0;//to N einai h diastash toy tetragwnikoy pinaka kai to error einai to epitrepto sfalma logw arithmwn double
	}
	char path1[100], path2[100];
	sprintf(path1,"./%s",argv[1]); 
	sprintf(path2,"./%s",argv[2]); 
//printf("path1 = %s\npath2 = %s\n",path1,path2);
	FILE *f1 = fopen(path1,"r");
	FILE *f2 = fopen(path2,"r");
	int N = atoi(argv[3]);
	float error = atof(argv[4]);
	int error_num=0;	//se posa diaforetika stoixeia toy pinaka diaferoyn
	float temp1,temp2;
	int i;
//printf("error = %f\n",error);
	for(i=0;i<N*N;i++){
		fscanf(f1,"%f",&temp1);
		fscanf(f2,"%f",&temp2);
		//if( (temp1-temp2 > error) || (temp1-temp2 < (-error)))
		if(  ((temp1>temp2)&&((temp1-temp2) > error)) || ((temp1<temp2)&&((temp2-temp1) > error)) ){
//			printf("%f - %f    =  %f\n",temp1,temp2,temp1-temp2);
			error_num++;
		}
	}
	printf("error_nums = %d\n",error_num);
}


