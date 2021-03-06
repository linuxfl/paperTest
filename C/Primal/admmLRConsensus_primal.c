#include <stdio.h>  
#include <stdlib.h> 
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>

#define P 8
#define NODE 64
//#define RHO 0.5

#define NUMDATASET (P*NODE)
#define CONNECTED 0
#define CYCLE 1
#define BUTTERFLY 2
#define ITERMAX 4000

typedef int VertexType;

/*every node infomation in Graph*/
typedef struct node
{
	VertexType vertex;
	VertexType *adjustV;
	int num;
}Node;

typedef struct graph{
	Node *vertex;
	int n;
}Graph;

typedef struct DataSet{
	gsl_matrix *A;
	gsl_matrix *b;
	gsl_matrix *x;
	gsl_matrix *y;
	gsl_matrix *oldx;
	gsl_matrix *sumx;
}DS;

void printMatirx(gsl_matrix *m){
	int i,j;
	for(i=0;i<m->size1;i++)
	{
		for(j = 0;j < m->size2;j++)
		{
			printf("%f ",gsl_matrix_get(m,i,j));
		}
		printf("\n");
	}
}

void getInverse(gsl_matrix *inverse,gsl_matrix *A)
{
	int n = A->size1;	
	int sign = 0;
	gsl_matrix *tmpA = gsl_matrix_alloc(n, n);
	gsl_matrix_memcpy(tmpA, A);
	gsl_permutation *p = gsl_permutation_alloc(n);
	gsl_linalg_LU_decomp(tmpA, p, &sign);
	gsl_linalg_LU_invert(tmpA, p, inverse);
	gsl_permutation_free(p);
	gsl_matrix_free(tmpA);
} 

float matrixDnrm2(gsl_matrix *a,gsl_matrix *b)
{
	float error = 0.0,solunorm = 0.0;
	int i;
	for(i=0;i<a->size1;i++){
		error += pow((gsl_matrix_get(a,i,0) - gsl_matrix_get(b,i,0)),2);
		solunorm += pow(gsl_matrix_get(b,i,0),2);
	}
	return sqrt(error)/sqrt(solunorm);
}

int admmTrain(float rho)
{
	int i,j,k = 0;
	double minerror;
	DS *data = (DS *)malloc(sizeof(DS)*NUMDATASET);
	FILE *fA,*fb,*fs;
	char strA[10],strb[10];
	gsl_matrix *eye;
	gsl_matrix *Atb;
	gsl_matrix *inverse;
	gsl_matrix *AtA;
	gsl_matrix *At;
	gsl_matrix *solution = gsl_matrix_calloc(3,1);
	fs = fopen("./data/solution.dat","r");
	gsl_matrix_fscanf(fs,solution);
	gsl_vector *error = gsl_vector_calloc(NUMDATASET);
	for(i = 0;i < NUMDATASET;i++){
		data[i].A = gsl_matrix_calloc(3,3);
		data[i].b = gsl_matrix_calloc(3,1);
		data[i].x = gsl_matrix_calloc(3,1);
		data[i].y = gsl_matrix_calloc(3,1);
		data[i].oldx = gsl_matrix_calloc(3,1);
		data[i].sumx = gsl_matrix_calloc(3,1);
		sprintf(strA,"./data/A%d.dat",i);
		sprintf(strb,"./data/b%d.dat",i);
		fA = fopen(strA,"r");
		fb = fopen(strb,"r");
		gsl_matrix_fscanf(fA,data[i].A);
		gsl_matrix_fscanf(fb,data[i].b);
	}
	eye = gsl_matrix_calloc(3,3);
	Atb = gsl_matrix_calloc(3,1);
	AtA = gsl_matrix_calloc(3,3);
	inverse = gsl_matrix_calloc(3,3);
	At = gsl_matrix_calloc(3,3);
	/* (At * A + 2 * size * rho * eye(n))\(At * b + rho * size * x +rho * sumx - y) */

	while(k < ITERMAX){
		for(i = 0;i < NUMDATASET;i++){

			gsl_matrix_set_identity(eye);
			gsl_matrix_transpose_memcpy(At,data[i].A);
			//printMatirx(At);
			/* (At * A + 2 * size * rho * eye(n))\*/

			gsl_matrix_scale(eye,2 * NUMDATASET * rho);
			gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,1.0,At,data[i].A,0,AtA);
			gsl_matrix_add(AtA,eye);
			getInverse(inverse,AtA);

			/*(At * b + rho * size * x +rho * sumx - y) */
			gsl_blas_dgemm(CblasTrans, CblasNoTrans,1.0,data[i].A,data[i].b,0,Atb); /*At * b*/
			gsl_matrix_scale(data[i].oldx,NUMDATASET * rho);/*rho * size * x*/
			
			gsl_matrix_scale(data[i].sumx,rho);
			gsl_matrix_add(Atb,data[i].oldx);
			gsl_matrix_sub(Atb,data[i].y);
			gsl_matrix_add(Atb,data[i].sumx);
			/*x update*/
			gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,1.0,inverse,Atb,0,data[i].x);
		}
		
		for(i=0;i<NUMDATASET;i++){
			gsl_matrix_set_zero(data[i].sumx);
			for(j=0;j<NUMDATASET;j++)
			{
				gsl_matrix_add(data[i].sumx,data[j].x);
			}
		}
		
		for(i=0;i<NUMDATASET;i++)
		{
			gsl_matrix_memcpy(data[i].oldx,data[i].x);
			/*y update:data[i].y = data[i].y + rho * (data[i].x * size - data[i].sumx)*/
			gsl_matrix_scale(data[i].x,NUMDATASET);
			gsl_matrix_sub(data[i].x,data[i].sumx);
			gsl_matrix_scale(data[i].x,rho);
			
			gsl_matrix_add(data[i].y,data[i].x);
			gsl_vector_set(error,i,matrixDnrm2(data[i].oldx,solution));
		}
		minerror = gsl_vector_min(error);
		
		//if( fabs(minerror - 10e-10) < 10e-10)
		if(minerror < 10e-9)
			break;
		//printf("iter:%d primal error :%0.15f\n",k,minerror);
		k++;
	}

	//release memory
	gsl_matrix_free(eye);
	gsl_matrix_free(Atb);
	gsl_matrix_free(AtA);
	gsl_matrix_free(At);
	gsl_matrix_free(solution);
	gsl_matrix_free(inverse);
	gsl_vector_free(error);

	for(i = 0;i < NUMDATASET;i++){
		gsl_matrix_free(data[i].A);
		gsl_matrix_free(data[i].b);
		gsl_matrix_free(data[i].x);
		gsl_matrix_free(data[i].y);
		gsl_matrix_free(data[i].oldx);
		gsl_matrix_free(data[i].sumx);
	}
	
	free(data);
	fclose(fA);
	fclose(fb);
	fclose(fs);
	return k;
}

int main(int argc,char **argv)
{
	char str[100];
	if(argc < 2)
	{
		printf("need two arguments!!!\n");
		return 0;
	}
	float rho = atof(argv[1]);
	//Graph *g = creatGraph(P,NODE,0);
	sprintf(str,"iter count :%d\n",admmTrain(rho));
	printf("%s",str);
	return 0;
}
