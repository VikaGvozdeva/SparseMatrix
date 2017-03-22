#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h> 
#include <stdio.h> 
#include <cmath>
#include <string.h>
#include <iostream>
#include <algorithm>

using namespace std;

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

	void COOMATRIX::Sort(int NNZ)
	{
		int i, j;
		int tmp1, tmp2;
		double tmp3;
		for (i = 0; i < NNZ - 1; i++)
			for (j = i + 1; j < NNZ; j++)
			{
				if (row_ind[i] > row_ind[j])
				{
					tmp1 = row_ind[i];
					row_ind[i] = row_ind[j];
					row_ind[j] = tmp1;
					tmp2 = col_ind[i];
					col_ind[i] = col_ind[j];
					col_ind[j] = tmp2;
					tmp3 = val[i];
					val[i] = val[j];
					val[j] = tmp3;
				}
				if (row_ind[i] == row_ind[j])
				{
					if (col_ind[i] > col_ind[j])
					{
						tmp1 = row_ind[i];
						row_ind[i] = row_ind[j];
						row_ind[j] = tmp1;
						tmp2 = col_ind[i];
						col_ind[i] = col_ind[j];
						col_ind[j] = tmp2;
						tmp3 = val[i];
						val[i] = val[j];
						val[j] = tmp3;
					}
				}
			}
	}
};


struct CRSMATRIX {

	int N;			
	int NNZ;		
	double *val;	
	int *col_ind;
	int *row_ptr;


	CRSMATRIX(int _NNZ, int _N)
	{
		int i;
		N = _N;
		NNZ = _NNZ;
		val = (double*)malloc(NNZ * sizeof(double));
		col_ind = (int*)malloc(NNZ * sizeof(int));
		row_ptr = (int*)malloc((N + 1) * sizeof(int));
		for (i = 0; i < N + 1; i++)
		{
			row_ptr[i] = 0;
		}
	}
	
	CRSMATRIX(const CRSMATRIX& Matrix)
	{
		int i;
		N = Matrix.N;
		NNZ = Matrix.NNZ;
		val = (double*)malloc(NNZ * sizeof(double));
		col_ind = (int*)malloc(NNZ * sizeof(int));
		row_ptr = (int*)malloc((N + 1) * sizeof(int));
		for (i = 0; i < N + 1; i++)
		{
			row_ptr[i] = 0;
		}
	
	}

	~CRSMATRIX()
	{
		free(val);
		free(col_ind);
		free(row_ptr);
	}

	void CRSMATRIX::PrintCRSMatrix(int _NNZ, int _N)
	{
		int i;
		N = _N;
		NNZ = _NNZ;
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

	int N;
	int NNZ;
	double *val;
	int *row_ind;
	int *col_ptr;


	CCSMATRIX(int _NNZ, int _N)
	{
		int i, j;
		N = _N;
		NNZ = _NNZ;
		val = (double*)malloc(NNZ * sizeof(double));
		row_ind = (int*)malloc(NNZ * sizeof(int));
		col_ptr = (int*)malloc((N + 1) * sizeof(int));
		for (i = 0; i < N + 1; i++)
		{
			col_ptr[i] = 0;
		}
	}


	CCSMATRIX(const CCSMATRIX& Matrix)
	{
		int i;
		N = Matrix.N;
		NNZ = Matrix.NNZ;
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
	{
		int i;
		int N = _N;
		int NNZ = _NNZ;
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


double * Matrix_VectorMultiplicationInCRS(CRSMATRIX* Matrix, double *Vector, int N)
{
	int i, j;
	double tmp;
	double * result = (double*)malloc(N * sizeof(double));
	for (i = 0; i < N; i++)
	{
		result[i] = 0;
	}
	if (N == Matrix->N)
	{
		for (i = 0; i < N; i++)
		{
			tmp = 0;
			for (j = Matrix->row_ptr[i]; j < Matrix->row_ptr[i + 1]; j++)
			{
				tmp += Matrix->val[j] * Vector[Matrix->col_ind[j]];
			}
			result[i] = tmp;
		}
	}
	return result;
}

double * Matrix_VectorMultiplicationInCCS(CCSMATRIX* Matrix, double *Vector, int N)
{
	int i, j, k;

	double * result = (double*)malloc(N * sizeof(double));
	for (k = 0; k < N; k++)
	{
		result[k] = 0;
	}
	if (N == Matrix->N)
	{
		for (i = 0; i < N; i++)
		{
			for (j = Matrix->col_ptr[i]; j < Matrix->col_ptr[i + 1]; j++)
			{
				result[Matrix->row_ind[j]] += Vector[i] * Matrix->val[j];
			}
			/*printf("step is %d", i);
			for (int f = 0; f < N; f++)
				printf("%lf \n  ", result[f]);
*/
		}
	}
	return result;
}



struct COMPDIAGMATRIX{

	int B;
	double** value;

	void GetB(COMPDIAGMATRIX* Matrix);
	

};
struct JDIAGMATRIX {

	int NNZ;
	int N;
	int numbdiag;
	int *perm;
	double *jdiag;
	int *col_ind;
	int *jd_ptr;

	JDIAGMATRIX(int _NNZ, int _N)// int _numbdiag)
	{
		//numbdiag = _numbdiag;
		N = _N;
		NNZ = _NNZ;
		int i;
		jdiag = (double*)malloc(NNZ * sizeof(double));
		col_ind = (int*)malloc(NNZ * sizeof(int));
		jd_ptr = (int*)malloc((N ) * sizeof(int));
		for (i = 0; i < N + 1; i++)
		{
			jd_ptr[i] = 0;
		}
	}
	void JDIAGMATRIX::PrintJDIAGMatrix(int _NNZ, int _N)
	{
		int i;
		N = _N;
		NNZ = _NNZ;
		printf("val:");
		for (i = 0; i < NNZ; i++)
			printf("%lf , ", jdiag[i]);
		printf(" \n col_ind:");
		for (i = 0; i < NNZ; i++)
			printf("%d , ", col_ind[i]);
		printf(" \n diag_ptr:");
		for (i = 0; i < N + 1; i++)
			printf("%d , ", jd_ptr[i]);
		for (i = 0; i < N + 1; i++)
			printf("%d , ", perm[i]);
		printf("\n");
	}

	~JDIAGMATRIX()
	{
		free(col_ind);
		free(jd_ptr);
		free(jdiag);
		free(perm);
	}


};
struct BCRSMATRIX{
};
struct SKYLINEMATRIX{

	int N;
	int NNZ;
	double *adiag;
	double *altr;
	double *autr;
	int* jptr;
	int *iptr;

	SKYLINEMATRIX(int _N, int _NNZ)
	{
		N = _N;
		NNZ = _NNZ;
		adiag = (double*)malloc(N * sizeof(double));
		altr = (double*)malloc((NNZ-N)/2 * sizeof(double));
		autr = (double*)malloc((NNZ - N) / 2 * sizeof(double));
		jptr = (int*)malloc((NNZ - N) / 2 * sizeof(int));
		iptr = (int*)malloc((N + 1) * sizeof(int));

	}
	~SKYLINEMATRIX()
	{
		
		free(adiag);
		free(autr);
		free(altr);
		free(iptr);
		free(jptr);
	}

	void SKYLINEMATRIX::PrintSLMatrix(int _NNZ, int _N)
	{
		int i;
		N = _N;
		NNZ = _NNZ;
		printf("adiag:");
		for (i = 0; i < N; i++)
			printf("%lf , ", adiag[i]);
		printf(" \n autr:");
		for (i = 0; i < (NNZ-N)/2; i++)
			printf("%d , ", autr[i]);
		printf(" \n altr:");
		for (i = 0; i < (NNZ - N) / 2; i++)
			printf("%d , ", altr[i]);
		for (i = 0; i < (N + 1); i++)
			printf("%d , ", iptr[i]);
		printf(" \n jptr:");
		for (i = 0; i < (NNZ - N) / 2; i++)
			printf("%d , ", jptr[i]);
		
		printf("\n");
	}

};


//struct COMPR_DIAG_MATRIX {
//
//	int N;
//	int NNZ;
//	int B;
//	int block_rows, block_cols; //number of rows or columns blocks
//	int r, c; //rows and columns in each subblock
//	double **val;
//	double *val_c;
//	int **col_ind;
//	int *col_ind_c;
//	int *row_ind_c;
//
//
//	COMPR_DIAG_MATRIX(int NNZ, int N)
//	{
//		int i;
//		val = (double**)malloc(N * sizeof(double*));
//		for (i = 0; i<N; i++)
//			val[i] = (double*)malloc(NNZ_max * sizeof(double));

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

CCSMATRIX* ConverterToCCS(COOMATRIX &Matrix)
{
	int i = 0, j = 0, k = 0, NNZ = 0, N = 0, tmp_ind = 0;
	NNZ = Matrix.NNZ;
	N = Matrix.N;
	CCSMATRIX* Mtx= new CCSMATRIX(NNZ,N);

	for (i = 0; i < NNZ; i++)
	{
		Mtx->val[i] = Matrix.val[i];
		Mtx->row_ind[i] = Matrix.row_ind[i];
	}

	for (j = 0; j < NNZ; j++)
	{
		tmp_ind = Matrix.col_ind[j];
		Mtx->col_ptr[++tmp_ind]++;
	}
	while (j < NNZ);

	for (k = 2; k < (Matrix.N + 1 ); k++)
	{
		Mtx->col_ptr[k] += Mtx->col_ptr[k - 1];
	}

	return Mtx;
}

CRSMATRIX* ConverterToCRS(COOMATRIX &Matrix)
{
	int i = 0, j = 0, k = 0, NNZ = 0, N = 0, tmp_ind = 0;
	NNZ = Matrix.NNZ;
	N = Matrix.N;
	CRSMATRIX* Mtx = new CRSMATRIX(NNZ, N);

	for (i = 0; i < NNZ; i++)
	{
		Mtx->val[i] = Matrix.val[i];
		Mtx->col_ind[i] = Matrix.col_ind[i];
	}

	for (j = 0; j < NNZ; j++)
	{
		tmp_ind = Matrix.row_ind[j];
		Mtx->row_ptr[++tmp_ind]++;
	}
	while (j < NNZ);

	for (k = 2; k < (Matrix.N + 1); k++)
	{
		Mtx->row_ptr[k] += Mtx->row_ptr[k - 1];
	}

	return Mtx;
}


void qs(int* arr, int* _arr, int first, int last)
{
	int i = first, j = last, x = arr[(first + last) / 2];
	int temp;

	do {
		while (arr[i] < x) i++;
		while (arr[j] > x) j--;

		if (i <= j) {
			if (arr[i] > arr[j])
			{
				temp = arr[i];
				arr[i] = arr[j];
				arr[j] = arr[i];

				temp = _arr[i];
				arr[i] = _arr[j];
				arr[j] = _arr[i];

			}
			i++;
			j--;
		}
	} while (i <= j);

	if (i < last)
		qs(arr, _arr, i, last);
	if (first < j)
		qs(arr, _arr, first, j);
}


JDIAGMATRIX* ConverterToJDIAG(COOMATRIX &Matrix)
{
	int i = 0, j = 0, k = 0, NNZ = 0, N = 0, tmp_ind = 0, maxval = 0;
	NNZ = Matrix.NNZ;
	N = Matrix.N;
	int *row_len = new int[N];

	JDIAGMATRIX* Mtx = new JDIAGMATRIX(NNZ, N);

	//??
	///Mtx->numbdiag = Matrix.NNZ_row[0];
	for (i = 0; i < N; i++)
	{
		//row_len[i] = Matrix.row_ptr[i + 1] - Matrix.row_ptr[i];
		row_len[i] = Matrix.NNZ_row[i];
		Mtx->perm[i] = i;
	}

	qs(row_len, Mtx->perm, 0, N - 1);
	maxval = row_len[0];
	Mtx->numbdiag = maxval;
	CRSMATRIX* Test = ConverterToCRS(Matrix);
	for (i = 0; i < maxval; i++)
	{
		Mtx->jd_ptr[i] = k;
		for (j = 0; j < N; j++)
		{
			if (i < Mtx->perm[j])
			{
				Mtx->jdiag[k] = Matrix.val[Test->row_ptr[Mtx->perm[j]] + i];
				Mtx->col_ind[k] = Matrix.col_ind[Test->row_ptr[Mtx->perm[j]] + i];
				k++;
			}
		}
	}

	return Mtx;
}

SKYLINEMATRIX* ConverterToSL(COOMATRIX &Matrix)
{
	int i = 0, j = 0, k = 0, NNZ = 0, N = 0, tmp_ind = 0;
	NNZ = Matrix.NNZ;
	N = Matrix.N;
	SKYLINEMATRIX* Mtx = new SKYLINEMATRIX(NNZ, N);
	//sort coo(rows)
	//autr
	for (i = 0; i < NNZ; i++)
	{
		//altr
		if (Matrix.col_ind[i] < Matrix.row_ind[i])
		{
			Mtx->altr[j++] = Matrix.val[i];
			//Mtx->altr[j++] = Matrix.val[Matrix.col_ind[i]];
			Mtx->jptr[j++] = Matrix.col_ind[i];
		}
		else if (Matrix.col_ind[i] = Matrix.row_ind[i])
		{
			Mtx->adiag[Matrix.col_ind[i]] = Matrix.val[i];
			//Mtx->adiag[Matrix.col_ind[i]] = Matrix.val[Matrix.col_ind[i]];
		}

	}
	Matrix.Sort(Matrix.NNZ);
	//äåëàåì ñîðò
	for (i = 0; i < NNZ; i++)
	{//autr
		if (Matrix.col_ind[i] > Matrix.row_ind[i])
		{
			Mtx->autr[j++] = Matrix.val[Matrix.col_ind[i]];
		}
	}
		//sort
		//altr
	return Mtx;
}
	int main()
{
	int N = 0;
	int NNZ = 0;
	double *v, *res;
	ReadNumberForMatrix(N, NNZ);
	COOMATRIX Matrix(NNZ, N);

	ReadMatrixInfo(Matrix);
	Matrix.PrintMatrix(Matrix.NNZ);
	//for CRS
	//Matrix.Sort(Matrix.NNZ);
	//Matrix.PrintMatrix(Matrix.NNZ);
	//CCSMATRIX* Test = ConverterToCÑS(Matrix);
	//CRSMATRIX* Test = ConverterToCRS(Matrix);
	/*CCSMATRIX* test = ConverterToCCS(Matrix);*/
	//test->PrintCCSMatrix(test->NNZ, test->N);
	//JDIAGMATRIX* test = ConverterToJDIAG(Matrix);
	//test->PrintJDIAGMatrix(test->NNZ, test->N);
	SKYLINEMATRIX* test = ConverterToSL(Matrix);
	test->PrintSLMatrix(test->NNZ, test->N);
//	v = (double*)malloc(Matrix.N * sizeof(double));
	//for (int i = 0; i < Matrix.N; i++)
		/*v[0] = 2;
		v[1] = 0;
		v[2] = 3;
		v[3] = 10;
		v[4] = 0;
		v[5] = 1;
		v[6] = 3;*/
	//res=Matrix_VectorMultiplicationInCRS(test, v, N);
//	res = Matrix_VectorMultiplicationInCCS(test, v, N);

	//printf("Result: \n");
	//for (int i = 0; i < Matrix.N; i++)
	//	printf("%lf  ", res[i]);
	system("pause");
}