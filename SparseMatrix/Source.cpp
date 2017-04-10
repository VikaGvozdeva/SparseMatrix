#define _CRT_SECURE_NO_WARNINGS

#include <Windows.h>
#include <stdlib.h> 
#include <stdio.h> 
#include <cmath>
#include <string.h>
#include <iostream>
#include <algorithm>
#include <time.h>

using namespace std;

#define MAX_LINE_LEN 1000000

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

	void PrintMatrix(int NNZ)
	{
		int i;
		for (i = 0; i < NNZ; i++)
		{
			printf("(%lf,%d,%d) , ", val[i], row_ind[i], col_ind[i]);
		}
		printf("\n");
	}

	void Sort(int NNZ)
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
void DeleteCOOMATRIX(COOMATRIX &Matrix)
{
	free(Matrix.val);
	free(Matrix.col_ind);
	free(Matrix.row_ind);
	free(Matrix.NNZ_row);
}

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
	
	/*CRSMATRIX(const CRSMATRIX& Matrix)
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
	
	}*/

	void PrintCRSMatrix(int _NNZ, int _N)
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
void DeleteCRSMATRIX(CRSMATRIX *Matrix)
{
	free(Matrix->val);
	free(Matrix->col_ind);
	free(Matrix->row_ptr);
}

CRSMATRIX* ConvertToCRS(COOMATRIX &Matrix)
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

	void PrintCCSMatrix(int _NNZ, int _N)
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

void DeleteCCSMATRIX(CCSMATRIX *Matrix)
{
	free(Matrix->val);
	free(Matrix->row_ind);
	free(Matrix->col_ptr);
}


CCSMATRIX* ConvertToCCS(COOMATRIX &Matrix)
{
	int i = 0, j = 0, k = 0, NNZ = 0, N = 0, tmp_ind = 0;
	NNZ = Matrix.NNZ;
	N = Matrix.N;
	CCSMATRIX* Mtx = new CCSMATRIX(NNZ, N);

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

	for (k = 2; k < (Matrix.N + 1); k++)
	{
		Mtx->col_ptr[k] += Mtx->col_ptr[k - 1];
	}

	return Mtx;
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

struct COMPDIAGMATRIX {

	int N;
	int NNZ;

	int B;
	double* diag;
	double** val;

	COMPDIAGMATRIX(int _N, int _NNZ, int _B)
	{
		int i;
		int j;
		N = _N;
		NNZ = _NNZ;
		B = _B;

		diag = (double*)malloc(B * sizeof(double*));
		val = (double**)malloc(N * sizeof(double*));
		for (i = 0; i < B; i++)
			val[i] = (double*)malloc(N * sizeof(double*));
		/*val = (double**)malloc(B * sizeof(double*));
		for (i = 0; i < N; i++)
			val[i] = (double*)malloc(N * sizeof(double*));*/
		for (i = 0; i < B; i++)
			for (j = 0; j < N; j++)
				val[i][j] = 0;
	}
	
	~COMPDIAGMATRIX()
	{
		int i;
		for (int i = 0; i < N; i++)
			free(val[i]);
		free(val);
		free(diag);
	}

	void PrintMatrix(int _N, int _NNZ, int _B)
	{
		int i, j;
		int N = _N;
		int B = _B;
		int NNZ = _NNZ;
		printf("val:");
		for (i = 0; i < B; i++)
		{
			printf("diag %d \n", i);
			for (j = 0; j < N; j++)
				printf("%lf , ", val[i][j]);
		}
		printf(" \n diag:");
		for (i = 0; i < B; i++)
			printf("%d , ", diag[i]);
	
		printf("\n");
	}
};
void qs(int* s_arr,int* _s_arr, int first, int last)
{
	int i = first, j = last, x = s_arr[(first + last) / 2], tmp = 0;

	do {
		/*while (s_arr[i] < x) i++;
		while (s_arr[j] > x) j--;*/
		while (s_arr[i] > x) i++;
		while (s_arr[j] < x) j--;

		if (i <= j) {
			if (s_arr[i] < s_arr[j])
			{
				tmp = s_arr[i];
				s_arr[i] = s_arr[j];
				s_arr[j] = tmp;

				tmp = _s_arr[i];
				_s_arr[i] = _s_arr[j];
				_s_arr[j] = tmp;

			}
			i++;
			j--;
		}
	} while (i <= j);

	if (i < last)
		qs(s_arr, _s_arr, i, last);
	if (first < j)
		qs(s_arr, _s_arr, first, j);
}



COMPDIAGMATRIX* ConvertToCompDiag(COOMATRIX &Matrix)
{
	int i = 0, j = 0, k = 0, l = 0, NNZ = 0, N = 0, diag_ind = 0, B = 0, tmp_ind = 0;
	NNZ = Matrix.NNZ;
	N = Matrix.N;

	int diag_numb = 0;
	diag_numb = Matrix.NNZ_row[0];
	printf("%d , ", diag_numb);
	for (i = 1; i < N; i++)
	{
		if (Matrix.NNZ_row[i] > diag_numb)
		{
			diag_numb = Matrix.NNZ_row[i];
			//where is the max of nnz (row)
			diag_ind = i;
		}
	}
	B = diag_numb;
	//printf("%d , ", diag_ind);
	COMPDIAGMATRIX* Mtx = new COMPDIAGMATRIX(N, NNZ, B);
	Matrix.Sort(Matrix.NNZ);
	//fill diag 0,diag_ind
	for (i = 0; i < NNZ; i++)
	{
		if (Matrix.row_ind[i] == diag_ind)
		{
			Mtx->val[diag_ind][j] = Matrix.val[i];
			Mtx->diag[j] = Matrix.col_ind[i] - Matrix.row_ind[i];
			j++;
		}
	}
	for (i = 0; i < NNZ; i++)
	{
		tmp_ind = Matrix.col_ind[i] - Matrix.row_ind[i];
		for (j = 0; j < diag_numb; j++)
		{
			if (tmp_ind == Mtx->diag[j])
				Mtx->val[Matrix.row_ind[i]][j] = Matrix.val[i];

		}
	}

	return Mtx;
}

double * Matrix_VectorMultiplicationInCompDiag(COMPDIAGMATRIX* Matrix, double *Vector, int N)
{
	int i, j;
	int tmp;
	int max;
	int B = Matrix->B;
	double * result = (double*)malloc(N * sizeof(double));
	for (i = 0; i < N; i++)
	{
		result[i] = 0;
	}

	for (j = 0; j < B; j++) {
		tmp = Matrix->diag[j];
		if (0 >(0 - tmp))
			i = 0;
		else i = 0 - tmp;

		if (0 > tmp)
			max = Matrix->N - 0;
		else max = Matrix->N - tmp;

		for (i; i < max; i++)
			result[i] += Matrix->val[i][j] * Vector[tmp + i];

	}
	return result;
}

struct JDIAGMATRIX {

	int NNZ;
	int N;
	int numbdiag;
	int *perm;
	double *jdiag;
	int *col_ind;
	int *jd_ptr;

	JDIAGMATRIX(int _N, int _NNZ)// int _numbdiag)
	{
		//numbdiag = _numbdiag;
		int i;
		N = _N;
		NNZ = _NNZ;
		jdiag = (double*)malloc(NNZ * sizeof(double));
		perm = (int*)malloc(N * sizeof(int));
		col_ind = (int*)malloc(NNZ * sizeof(int)); 
		jd_ptr = (int*)malloc((N ) * sizeof(int));
		for (i = 0; i < N + 1; i++)
		{
			jd_ptr[i] = 0;
		}
	}
	void PrintMatrix(int _NNZ, int _N)
	{
		int i;
		N = _N;
		NNZ = _NNZ;
		printf(" \n jdiag:");
		for (i = 0; i < NNZ; i++)
			printf("%d , ", jdiag[i]);
		printf(" \n col_ind:");
		for (i = 0; i < NNZ; i++)
			printf("%d , ", col_ind[i]);
		printf(" \n diag_ptr:");
		for (i = 0; i < N + 1; i++)
			printf("%d , ", jd_ptr[i]);
		printf(" \n perm:");
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

JDIAGMATRIX* ConvertToJDIAG(COOMATRIX &Matrix)
{//dosort
	int i = 0, j = 0, k = 0, NNZ = 0, N = 0, tmp_ind = 0, maxval = 0;
	NNZ = Matrix.NNZ;
	N = Matrix.N;
	int *row_len = new int[N];

	JDIAGMATRIX* Mtx = new JDIAGMATRIX(N, NNZ);

	///Mtx->numbdiag = Matrix.NNZ_row[0];
	for (i = 0; i < N; i++)
	{
		//row_len[i] = Matrix.row_ptr[i + 1] - Matrix.row_ptr[i];
		row_len[i] = Matrix.NNZ_row[i];
		Mtx->perm[i] = i;
	}
	/*printf(" \n perm----------");
	for (i = 0; i < N; i++)
		printf("%d , ", Mtx->perm[i]);

	printf(" \n row_len----------");
	for (i = 0; i < N; i++)
		printf("%d , ", row_len[i]);
	qs(row_len, Mtx->perm, 0, N-1);
	printf(" \n perm----------");
	for (i = 0; i < N; i++)
		printf("%d , ", Mtx->perm[i]);

	printf(" \n row_len----------");
	for (i = 0; i < N; i++)
		printf("%d , ", row_len[i]);*/
	maxval = row_len[0];
	Mtx->numbdiag = maxval;
	double** val = (double**)malloc(N * sizeof(double*));
	for (i = 0; i < N; i++)
		val[i] = (double*)malloc(maxval * sizeof(double*));
	//for (i = 0; i < N; i++)
	//	for (j = 0; j < Matrix.NNZ_row[i]; j++)
	//	{
	//		val[i][j]=
	//	}
	/*for (i=0; i<N; i++)
		Mtx->perm[i]=*/
	/*CRSMATRIX* Test = ConverterToCRS(Matrix);
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
*/
	return Mtx;
}
double*  Matrix_VectorMultiplicationInJDIAG(JDIAGMATRIX* Matrix, double *Vector, int N)
{
	int i = 0, j = 0, k = 0, NNZ = 0, tmp_ind = 0, maxval = 0, upper = 0;
	NNZ = Matrix->NNZ;

	double * result = (double*)malloc(N * sizeof(double));
	for (i = 0; i < N; i++)
	{
		result[i] = 0;
	}


	for (int j = 0; j < Matrix->numbdiag; j++)
	{
		for (int i = Matrix->jd_ptr[j]; i < Matrix->jd_ptr[j + 1]; i++)
		{
			tmp_ind = Matrix->col_ind[i];
			result[i] += Matrix->jdiag[i] * Vector[tmp_ind];
		}
	}
	return result;
}

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

	SKYLINEMATRIX(int _NNZ, int _N)
	{
		int i;
		N = _N;
		NNZ = _NNZ;
		adiag = (double*)malloc(N * sizeof(double));
		altr = (double*)malloc((NNZ - N)/2 * sizeof(double));
		autr = (double*)malloc((NNZ - N) / 2 * sizeof(double));
		jptr = (int*)malloc((NNZ - N) / 2 * sizeof(int));
		iptr = (int*)malloc((N + 1) * sizeof(int));
		for (i = 0; i < (NNZ - N) / 2; i++)
		{
			altr[i] = 0;
			autr[i] = 0;
			jptr[i] = 0;
		}
		for (i = 0; i < N + 1; i++)
			iptr[i] = 0;

	}
	bool IsSymmetric(COOMATRIX &M)
	{
		int N = M.N;
		int NNZ = M.NNZ;
		int i = 0, l = 0;
		double* row = (double*)malloc((NNZ - N) / 2 * sizeof(double));
		double* col = (double*)malloc((NNZ - N) / 2 * sizeof(double));
		bool* check = (bool*)malloc((NNZ - N) / 2 * sizeof(bool));
		for (i = 0; i < NNZ; i++)
		{
			if (M.row_ind[i] > M.col_ind[i])
			{
				row[l] = M.row_ind[i];
				col[l] = M.col_ind[i];
				l++;
			}

		}
		l = 0;
		for (i = 0; i < NNZ; i++)
		{
			if (M.row_ind[i] < M.col_ind[i])
			{
				if ((M.col_ind[i] == col[l]) && (M.row_ind[i] == row[l]))
				{
					check[l] = true;
					l++;
				}
				else
				{
					check[l] = false;
					l++;
				}
			}
	
		}
		for (i = 0; i < (NNZ - N) / 2; i++)
		{
			if (check[i] == false)
				return false;
		}
		return true;

	}

	void PrintMatrix(int _NNZ, int _N)
	{
		int i;
		N = _N;
		NNZ = _NNZ;
		printf("adiag:");
		for (i = 0; i < N; i++)
			printf("%lf , ", adiag[i]);
		printf(" \n autr:");
		for (i = 0; i < (NNZ-N)/2; i++)
			printf("%lf , ", autr[i]);
		printf(" \n altr:");
		for (i = 0; i < (NNZ - N) / 2; i++)
			printf("%lf , ", altr[i]);
		printf(" \n iptr:");
		for (i = 0; i < (N + 1); i++)
			printf("%d , ", iptr[i]);
		printf(" \n jptr:");
		for (i = 0; i < (NNZ - N) / 2; i++)
			printf("%d , ", jptr[i]);
		printf("\n");
	}

};

void DeleteSKYLINEMATRIX(SKYLINEMATRIX *Matrix)
{

	free(Matrix->adiag);
	free(Matrix->autr);
	free(Matrix->altr);
	free(Matrix->iptr);
	free(Matrix->jptr);
}

SKYLINEMATRIX* ConvertToSL(COOMATRIX &Matrix)
{
	int i = 0, j = 0, k = 0, l = 0, NNZ = 0, N = 0, tmp_ind = 0, ad = 0;
	NNZ = Matrix.NNZ;
	N = Matrix.N;
	SKYLINEMATRIX* Mtx = new SKYLINEMATRIX(NNZ, N);

	for (i = 0; i < NNZ; i++)
	{
		if (Matrix.col_ind[i] == Matrix.row_ind[i])
		{
			Mtx->adiag[Matrix.col_ind[i]] = Matrix.val[i];
		}
		else if (Matrix.col_ind[i] > Matrix.row_ind[i])
		{
			Mtx->autr[j] = Matrix.val[i];
			Mtx->jptr[l] = Matrix.row_ind[i];
			//Mtx->iptr[Matrix.row_ind[i+1]]++;
			//Mtx->iptr[Matrix.row_ind[i]]++;
			j++;
			l++;
		}

	}
	Matrix.Sort(Matrix.NNZ);
	for (i = 0; i < NNZ; i++)
	{
		if (Matrix.col_ind[i] < Matrix.row_ind[i])
		{
			Mtx->altr[k] = Matrix.val[i];
			k++;
			Mtx->iptr[Matrix.row_ind[i]+1]++;
		}
	}

	for (k = 2; k < (Matrix.N + 1); k++)
	{
		Mtx->iptr[k] += Mtx->iptr[k - 1];
	}
	return Mtx;
}

double * Matrix_VectorMultiplicationInSL(SKYLINEMATRIX* Matrix, double *Vector, int N)
{
	int i, j;
	double tmp;
	double * result = (double*)malloc(N * sizeof(double));
	for (i = 0; i < N; i++)
	{
		result[i] = 0;
	}
	for (i = 0; i < N; i++)
	{
		result[i] = Vector[i] * Matrix->adiag[i];
	}
	for (i = 0; i < N; i++)
		for (j = Matrix->iptr[i]; j < Matrix->iptr[i + 1]; j++)
		{
			result[i] += Vector[Matrix->jptr[j]] * Matrix->altr[j];
			result[Matrix->jptr[j]] += Vector[i] * Matrix->autr[j];
		}

	return result;
}
void ReadMatrixInfo(COOMATRIX& Matrix, char *name)//������ ������� �� �����
{
	FILE* f;
	int i;
	char* line;
	char* p = NULL;
	//char name[256] = "ck656.mtx";

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

void ReadNumberForMatrix(int& N, int& NNZ, char *name) //������ ������� � ���-�� ��������� ���������
{
	FILE* f;
	char* line;
	char* p = NULL;
	//char name[256] = "ck656.mtx";

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
double SearchMax_double(double *arr, int N) // ïîèñê ìàêñèìóìà äëÿ double
{
	int i;
	double max_arr = arr[0];
	for (i = 1; i < N; i++)
	{
		if (max_arr < arr[i])
			max_arr = arr[i];
	}
	return max_arr;
}
double CheckCorrectness(double * my_mult, double * mkl_mult, int N)  //ïðîâåðêà êîððåêòíîñòè óìíîæåíèÿ,íóæíî åùå ïîëó÷èòü ðåçóëüòàò â MKL
{
	int i;
	double res;
	double * arr_abs;
	arr_abs = (double*)malloc(N * sizeof(double));
	for (i = 0; i < N; i++)
		arr_abs[i] = abs(my_mult[i] - mkl_mult[i]);
	res = SearchMax_double(arr_abs, N);
	free(arr_abs);
	return res;
}

int main(int argc, char** argv)
{





	char *fileName;
	double timer;
	int i;
	double* v;
	int N = 0;
	int NNZ = 0;
	FILE *fp = fopen("data.dt", "w");
	LARGE_INTEGER frequency, start, finish;

	fileName = (char*)malloc(256 * sizeof(char));
	strcpy(fileName, argv[1]);

	fprintf(fp, "Matrix file name: %s\n\n", fileName);

	QueryPerformanceFrequency(&frequency);
	QueryPerformanceCounter(&start);
	ReadNumberForMatrix(N, NNZ, fileName);
	COOMATRIX Matrix(NNZ, N);
	ReadMatrixInfo(Matrix, fileName);
	QueryPerformanceCounter(&finish);
	timer = (finish.QuadPart - start.QuadPart) * 1000.0f / frequency.QuadPart;
	fprintf(fp, "Time Read Matrix Info: \t%lf\n\n", timer);

	//Matrix.Sort(Matrix.NNZ);

	QueryPerformanceFrequency(&frequency);
	QueryPerformanceCounter(&start);
	CCSMATRIX* Mtx_CCS = ConvertToCCS(Matrix);
	QueryPerformanceCounter(&finish);
	timer = (finish.QuadPart - start.QuadPart) * 1000.0f / frequency.QuadPart;
	fprintf(fp, "Time Convert in \n\nCCS: \t\t%lf\n", timer);

	QueryPerformanceFrequency(&frequency);
	QueryPerformanceCounter(&start);
	SKYLINEMATRIX* Mtx_SL = ConvertToSL(Matrix);
	QueryPerformanceCounter(&finish);
	timer = (finish.QuadPart - start.QuadPart) * 1000.0f / frequency.QuadPart;
	fprintf(fp, "Time Convert in SL \t%lf\n", timer);

	Matrix.Sort(Matrix.NNZ);

	QueryPerformanceFrequency(&frequency);
	QueryPerformanceCounter(&start);
	CRSMATRIX* Mtx_CRS = ConvertToCRS(Matrix);
	QueryPerformanceCounter(&finish);
	timer = (finish.QuadPart - start.QuadPart) * 1000.0f / frequency.QuadPart;
	fprintf(fp, "Time Convert in CRS \t%lf\n", timer);

	v = (double*)malloc(Matrix.N * sizeof(double));
	for (i = 0; i < Matrix.N; i++)
		v[i] = 1;

	QueryPerformanceFrequency(&frequency);
	QueryPerformanceCounter(&start);
	Matrix_VectorMultiplicationInCCS(Mtx_CCS, v, Mtx_CCS->N);
	QueryPerformanceCounter(&finish);
	timer = (finish.QuadPart - start.QuadPart) * 1000.0f / frequency.QuadPart;
	fprintf(fp, "\nTime Matrix-Vector multiplication in \n\nCCS: \t\t%lf\n", timer);

	QueryPerformanceFrequency(&frequency);
	QueryPerformanceCounter(&start);
	Matrix_VectorMultiplicationInSL(Mtx_SL, v, Mtx_SL->N);
	QueryPerformanceCounter(&finish);
	timer = (finish.QuadPart - start.QuadPart) * 1000.0f / frequency.QuadPart;
	fprintf(fp, "Time Matrix-Vector multiplication in SL \t%lf\n", timer);

	QueryPerformanceFrequency(&frequency);
	QueryPerformanceCounter(&start);
	Matrix_VectorMultiplicationInCRS(Mtx_CRS, v, Mtx_CRS->N);
	QueryPerformanceCounter(&finish);
	timer += (finish.QuadPart - start.QuadPart) * 1000.0f / frequency.QuadPart;
	fprintf(fp, "Time Matrix-Vector multiplication in CRS \t%lf\n", timer);



	DeleteCOOMATRIX(Matrix);
	DeleteCOOMATRIX(Matrix);
	DeleteCRSMATRIX(Mtx_CRS);
	DeleteCCSMATRIX(Mtx_CCS);
	DeleteSKYLINEMATRIX(Mtx_SL);
	fclose(fp);
}







//	int i;
//	int N = 0;
//	int NNZ = 0;
//	double *res_CCS;
//	double *res_CRS;
//	double *res_SL;
//	double *v;
//	clock_t start, finish;
//
//	start = clock();
//	ReadNumberForMatrix(N, NNZ);
//	COOMATRIX Matrix(NNZ, N);
//	ReadMatrixInfo(Matrix);
//	finish = clock();
//	printf(" Time Read Matrix Info = %lf\n", (long double)(finish - start) / CLOCKS_PER_SEC);
//	//if (N > 20)
//	//	Matrix.PrintMatrix(Matrix.NNZ);
//	////Matrix.Sort(Matrix.NNZ);
//	//if (N > 20)
//	//	Matrix.PrintMatrix(Matrix.NNZ);
//
//	v = (double*)malloc(Matrix.N * sizeof(double));
//	for (i = 0; i < Matrix.N; i++)
//		v[i] = 1;
//
//	start = clock();
//	CCSMATRIX* Mtx_CCS = ConvertToCCS(Matrix);
//	finish = clock();
//	printf(" Time Convert in CCS = %lf\n", (long double)(finish - start) / CLOCKS_PER_SEC);
//
//	start = clock();
//	res_CCS = Matrix_VectorMultiplicationInCCS(Mtx_CCS, v, Mtx_CCS->N);
//	finish = clock();
//	printf(" Time Matrix-Vector multiplication in CCS = %lf\n", (long double)(finish - start) / CLOCKS_PER_SEC);
//
//
//	start = clock();
//	SKYLINEMATRIX* Mtx_SL = ConvertToSL(Matrix);
//	finish = clock();
//	printf(" Time Convert in SL = %lf\n", (long double)(finish - start) / CLOCKS_PER_SEC);
//
//	start = clock();
//	res_SL = Matrix_VectorMultiplicationInSL(Mtx_SL, v, Mtx_SL->N);
//	finish = clock();
//	printf(" Time Matrix-Vector multiplication in SL = %lf\n", (long double)(finish - start) / CLOCKS_PER_SEC);
//
//	Matrix.Sort(Matrix.NNZ);
//
//	start = clock();
//	CRSMATRIX* Mtx_CRS = ConvertToCRS(Matrix);
//	finish = clock();
//	printf(" Time Convert in CRS = %lf\n", (long double)(finish - start) / CLOCKS_PER_SEC);
//
//	start = clock();
//	res_CRS = Matrix_VectorMultiplicationInCRS(Mtx_CRS, v, Mtx_CRS->N);
//	finish = clock();
//	printf(" Time Matrix-Vector multiplication in CRS = %lf\n", (long double)(finish - start) / CLOCKS_PER_SEC);
//
//	DeleteCOOMATRIX(Matrix);
//	DeleteCRSMATRIX(Mtx_CRS);
//	DeleteCCSMATRIX(Mtx_CCS);
//	DeleteSKYLINEMATRIX(Mtx_SL);
//	system("pause");
//
//}
//
//
//	int N = 0;
//	int NNZ = 0;
//	double *v, *res;
//	ReadNumberForMatrix(N, NNZ);
//	COOMATRIX Matrix(NNZ, N);
//
//	ReadMatrixInfo(Matrix);
//	Matrix.PrintMatrix(Matrix.NNZ);
//	//for CRS
//	//Matrix.Sort(Matrix.NNZ);
//	//Matrix.PrintMatrix(Matrix.NNZ);
//	//CCSMATRIX* Test = ConverterToCСS(Matrix);
//	//CRSMATRIX* Test = ConverterToCRS(Matrix);
//	/*CCSMATRIX* test = ConverterToCCS(Matrix);*/
//	//Test->PrintCRSMatrix(Test->NNZ, Test->N);
//	JDIAGMATRIX* test = ConverterToJDIAG(Matrix);
//	test->PrintMatrix(test->NNZ, test->N);
//	//SKYLINEMATRIX* test = ConverterToSL(Matrix);
//	//test->PrintMatrix(test->NNZ, test->N);
////	v = (double*)malloc(Matrix.N * sizeof(double));
////	res = (double*)malloc(Matrix.N * sizeof(double));
////	for (int i = 0; i < Matrix.N; i++)
////	{
////		v[i] = 1;
////	}
////
////	res=Matrix_VectorMultiplicationInJDIAG(test, v, N);
//////	res = Matrix_VectorMultiplicationInCCS(test, v, N);
////
////	printf("Result: \n");
////	for (int i = 0; i < Matrix.N; i++)
////		printf("%lf  ", res[i]);
//	system("pause");
///}