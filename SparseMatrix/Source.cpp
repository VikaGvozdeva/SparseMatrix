#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h> 
#include <stdio.h> 
#include <cmath>
#include <string.h>

#define MAX_LINE_LEN 1000

struct COOMATRIX
{
	int N;
	int NNZ;
	double *val;
	int *col_ind;
	int *row_ind;
	int *NNZ_row;
	COOMATRIX(int _NNZ, int _N)
	{
		int i;
		N = _N;
		NNZ = _NNZ;
		val = (double*)malloc(NNZ * sizeof(double));
		col_ind = (int*)malloc(NNZ * sizeof(int));
		row_ind = (int*)malloc(NNZ * sizeof(int));
		NNZ_row = (int*)malloc(N * sizeof(int));
		for (i = 0; i < N; i++)
			NNZ_row[i] = 0;
	}
	~COOMATRIX()
	{
		free(val);
		free(col_ind);
		free(row_ind);
		free(NNZ_row);
	}
	void COOMATRIX::PrintMatrix(int NNZ)
	{
		int i;
		for (i = 0; i < NNZ; i++)
		{
			printf("(%lf,%d,%d) , ", val[i], row_ind[i], col_ind[i]);
		}
		printf("\n");
	}
};

struct CRSMATRIX {

	int N;			//dimension
	int NNZ;		//numbers of NNZ elements
	double *val;	//array of NNZ elements
	int *col_ind;
	int *row_ptr;


	CRSMATRIX(int _NNZ, int _N)
	{
		int N = _N;
		int NNZ = _NNZ;
		val = (double*)malloc(NNZ * sizeof(double));
		col_ind = (int*)malloc(NNZ * sizeof(int));
		row_ptr = (int*)malloc((N + 1) * sizeof(int));
	}
	CRSMATRIX(const CRSMATRIX& Matrix)
	{
		N = Matrix.N;
		NNZ = Matrix.NNZ;
		int i;
		val = (double*)malloc(NNZ * sizeof(double));
		col_ind = (int*)malloc(NNZ * sizeof(int));
		row_ptr = (int*)malloc((N + 1) * sizeof(int));
	}

	~CRSMATRIX()
	{
		free(val);
		free(col_ind);
		free(row_ptr);
	}

	void CRSMATRIX::PrintCRSMatrix(int _NNZ, int _N)
	//CRSMATRIX PrintCRSMatrix(int _NNZ, int _N)
	{
		int N = _N;
		int NNZ = _NNZ;
		int i;
		printf("val:");
		for (i = 0; i < NNZ; i++)
			printf("%lf , ", val[i]);
		printf(" \n col_ind:");
		for (i = 0; i < NNZ; i++)
			printf("%d , ", col_ind[i]);
		printf(" \n row_ptr:");
		for (i = 0; i < N + 1; i++)
			printf("%d , ", row_ptr[i]);
		printf("\n");
	}
};

struct CCSMATRIX {

	int N;			//dimension
	int NNZ;		//numbers of NNZ elements
	double *val;	//array of NNZ elements
	int *row_ind;
	int *col_ptr;


	CCSMATRIX(int _NNZ, int _N)
	{
		int N = _N;
		int NNZ = _NNZ;
		val = (double*)malloc(NNZ * sizeof(double));
		row_ind = (int*)malloc(NNZ * sizeof(int));
		col_ptr = (int*)malloc((N + 1) * sizeof(int));
	}
	CCSMATRIX(const CCSMATRIX& Matrix)
	{
		N = Matrix.N;
		NNZ = Matrix.NNZ;
		int i;
		val = (double*)malloc(NNZ * sizeof(double));
		row_ind = (int*)malloc(NNZ * sizeof(int));
		col_ptr = (int*)malloc((N + 1) * sizeof(int));
	}

	~CCSMATRIX()
	{
		free(val);
		free(row_ind);
		free(col_ptr);
	}

	void CCSMATRIX::PrintCCSMatrix(int _NNZ, int _N)
	//CCSMATRIX PrintCCSMatrix(int _NNZ, int _N)
	{
		int N = _N;
		int NNZ = _NNZ;
		int i;
		printf("val:");
		for (i = 0; i < NNZ; i++)
			printf("%lf , ", val[i]);
		printf(" \n row_ind:");
		for (i = 0; i < NNZ; i++)
			printf("%d , ", row_ind[i]);
		printf(" \n col_ptr:");
		for (i = 0; i < N + 1; i++)
			printf("%d , ", col_ptr[i]);
		printf("\n");
	}
};


void ReadMatrixInfo(COOMATRIX &Matrix)
{
	FILE *f;
	int  i;
	char *line;
	char *p = NULL;
	char name[256] = "b1_ss.mtx";

	f = fopen(name, "r");
	if (f == NULL)
	{
		printf("%s- File Not Found!\n", name);
	}
	line = (char*)malloc((MAX_LINE_LEN) * sizeof(char));
	do
		fgets(line, MAX_LINE_LEN, f);
	while (line[0] == '%');
	p = strtok(line, " ");
	p = strtok(NULL, " ");
	p = strtok(NULL, " ");
	for (i = 0; i < Matrix.NNZ; i++)
	{
		fgets(line, MAX_LINE_LEN, f);
		p = strtok(line, " ");
		Matrix.row_ind[i] = atoi(p) - 1;
		Matrix.NNZ_row[Matrix.row_ind[i]]++;
		p = strtok(NULL, " ");
		Matrix.col_ind[i] = atoi(p) - 1;
		p = strtok(NULL, " ");
		Matrix.val[i] = atof(p);
	}
	free(line);
	fclose(f);
}
void ReadNumberForMatrix(int &N, int &NNZ)
{
	FILE *f;
	char *line;
	char *p = NULL;
	char name[256] = "b1_ss.mtx";

	f = fopen(name, "r");
	if (f == NULL)
		printf("%s- File Not Found!\n", name);
	line = (char*)malloc((MAX_LINE_LEN) * sizeof(char));
	do
		fgets(line, MAX_LINE_LEN, f);
	while (line[0] == '%');
	p = strtok(line, " ");
	N = atoi(p);
	p = strtok(NULL, " ");
	p = strtok(NULL, " ");
	NNZ = atoi(p);
	free(line);
}
int SearchMax(int *arr, int N)
{
	int i;
	int max_arr = arr[0];
	for (i = 1; i < N; i++)
	{
		if (max_arr < arr[i])
			max_arr = arr[i];
	}
	return max_arr;
}

CCSMATRIX* ConverterToCÑS(COOMATRIX &Matrix)
{
	int i = 0, j = 0, k = 1, NNZ_per_row = 0, NNZ = 0, N = 0;
	NNZ = Matrix.NNZ;
	N = Matrix.N;
	CCSMATRIX Mtx(NNZ, N);

	for (i = 0; i < NNZ; i++)
	{
		Mtx.val[i] = Matrix.val[i];
		printf("%lf , ", Mtx.val[i]);
		Mtx.row_ind[i] = Matrix.row_ind[i];
		printf("%d , ", Mtx.row_ind[i]);
	}

	Mtx.col_ptr[0] = 0;
	Mtx.col_ptr[N + 1] = NNZ;

}

int main()
{
	int N = 0;
	int NNZ = 0;
	ReadNumberForMatrix(N, NNZ);
	COOMATRIX Matrix(NNZ, N);
	ReadMatrixInfo(Matrix);
	Matrix.PrintMatrix(Matrix.NNZ);
	//MATRIX_ELLPACK* Matrix_E = ConverterInELLPACK(Matrix);
	//Matrix_E->PrintMatrix(Matrix_E->NNZ_max);
	system("pause");
}