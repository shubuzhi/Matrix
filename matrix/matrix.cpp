// matrix.cpp: 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
using namespace std;


typedef struct {
	double **m;//the indiviual entries in the matrix
	int l, c;//row and column
} dmatrix_t;

void nrerror(char error_text[])

{
	fprintf(stderr, "Run-time error...\n");
	fprintf(stderr, "%s\n", error_text);

	exit(1);
}

dmatrix_t *dmat_alloc(int l, int c)

{
	//Copies the content of one matrix1 into matrix2
	int j;
	dmatrix_t *new_matrix;
	new_matrix = (dmatrix_t *)malloc(sizeof(dmatrix_t));

	new_matrix->c = c;
	new_matrix->l = l;
	new_matrix->m = (double **)malloc(new_matrix->l * sizeof(double));
	if (!new_matrix->m) {
		nrerror("MATRIX.cpp: allocation failure");
	}
	new_matrix->m -= 1;

	for (j = 0; j < new_matrix->l; j++) {
		new_matrix->m[j] = (double *)malloc(new_matrix->c * sizeof(double));
		if (!new_matrix->m[j]) {
			nrerror("MATRIX.H: allocation failure");
		}
		new_matrix->m[j] -= 1;
	}

	return new_matrix;
}

void write_dmatrix(dmatrix_t *M)

{
	int i, j;
	printf("\n");
	for (i = 0; i < M->l; i++)
	{
		printf("\t\t");
		for (j = 0; j < M->c; j++)
		{
			printf("%7.4f\t", M->m[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

void delete_dmatrix(dmatrix_t *A)

{
	int i;
	for (i = 0; i < A->l; i++)
	{
		free(A->m[i]);
	}
	free(A->m);
	free(A);
}

dmatrix_t *dmat_init(dmatrix_t *A, double a)

{
	int i, j;

	for (i = 0; i < A->l; i++) {
		for (j = 0; j < A->c; j++) {
			A->m[i][j] = a;
		}
	}
	return(A);
}

void dmat_add(dmatrix_t *result,dmatrix_t *A, dmatrix_t *B)

{
	int i, j;

	if ((*A).l != (*B).l || (*A).c != (*B).c) {
		nrerror("MATRIX.cpp: incompatible matrix sizes");
	}
	
	for (i = 0; i < A->l; i += 1)
	{
		for (j = 0; j < A->c; j += 1)
		{
			result->m[i][j] = A->m[i][j] + B->m[i][j];
		}
	}
}

void dmat_sub(dmatrix_t *result, dmatrix_t *A, dmatrix_t *B)

{
	int i, j;

	if ((*A).l != (*B).l || (*A).c != (*B).c) {
		nrerror("MATRIX.cpp: incompatible matrix sizes");
	}

	for (i = 0; i < A->l; i += 1)
	{
		for (j = 0; j < A->c; j += 1)
		{
			result->m[i][j] = A->m[i][j] - B->m[i][j];
		}
	}
}

void dmat_scalar_mult(dmatrix_t *result, dmatrix_t *A, double a)

{
	int i, j;

	for (i = 0; i < A->l; i += 1) {
		for (j = 0; j < A->c; j += 1) {
			result->m[i][j] = A->m[i][j] * a;
		}
	}
}

dmatrix_t *dmat_mult(const dmatrix_t *A, const dmatrix_t *B)

{
	dmatrix_t *C;
	double sum = 0;
	int i, j, k;
	
	if ((*A).c != (*B).l) {
		nrerror("MATRIX.cpp: incompatible matrix sizes");
	}
	C = dmat_alloc((*A).l, (*B).c);
	
	for (i = 0; i < A->l; i += 1) {
		for (j = 0; j < B->c; j += 1) {
			for (k = 0; k < A->c; k += 1) {
				sum += A->m[i][k] * B->m[k][j];
			}
			C->m[i][j] = sum;
		}
	}
	return(C);
}

dmatrix_t *dmat_transpose(dmatrix_t *A)

{
	dmatrix_t *B = dmat_alloc(A->l, A->c);
	int i, j;

	for (i = 0; i < A->l; i += 1) {
		for (j = 0; j < A->c; j += 1) {
			B->m[j][i] = A->m[i][j];
		}
	}
	return(B);
}

double ddot_product(dmatrix_t *A, dmatrix_t *B)

{
	dmatrix_t *C;

	C = (dmatrix_t *)malloc(sizeof(dmatrix_t));

	if ((*A).c == (*B).c && (*A).l == 1 && (*B).l == 1) {
		C = dmat_mult(A, dmat_transpose(B));
	}
	else if ((*A).c == (*B).l && (*A).l == 1 && (*B).c == 1) {
		C = dmat_mult(A, B);
	}
	else if ((*A).l == (*B).c && (*A).c == 1 && (*B).l == 1) {
		C = dmat_mult(B, A);
	}
	else if ((*A).l == (*B).l && (*A).c == 1 && (*B).c == 1) {
		C = dmat_mult(dmat_transpose(A), B);
	}
	else nrerror("MATRIX.cpp: Incompatible matrix sizes");
	return((*C).m[1][1]);
}

dmatrix_t *sub_matrix(dmatrix_t *A, int row, int col)

{
	int i, j, k, l;
	dmatrix_t *B;

	if (row < 1 || row > A->l || col < 1 || col > A->c || A->c < 2 || A->l < 2) {
		nrerror("MATRIX.cpp: erroneous indices");
	}
	B = dmat_alloc(A->l - 1, A->c - 1);

	for (i = 0, k = 0; i < A->l; i++) {
		for (j = 0, l = 0; j < A->c; j++) {
			if (j != col && i != row) {
				B->m[k][l] = A->m[i][j];
			}
			if (j != col) l++;
		}
		if (i != row) k++;
	}
	return B;
}

double determinant(dmatrix_t *A)

{
	int i;
	double det;

	if ((*A).l < 1 || (*A).c < 1) {
		nrerror("MATRIX.cpp: erroneous matrix size");
	}
	else if ((*A).l != (*A).c) {
		nrerror("MATRIX.cpp: not a square matrix");
	}
	else if ((*A).l == 1) {
		det = (*A).m[1][1];
	}
	else {
		det = 0;
		for (i = 0; i < A->c; i++) {
			det += pow(-1.0, i + 1)*A->m[1][i] * determinant(sub_matrix(A, 1, i));
		}
	}
	return det;
}

dmatrix_t *cofactor(dmatrix_t *A)

{
	int i, j;
	dmatrix_t *B;

	B = dmat_alloc(A->l, A->c);

	for (i = 0; i < A->l; i++) {
		for (j = 0; j < A->c; j++) {
			B->m[i][j] = pow(-1.0, i + j)*determinant(sub_matrix(A, i, j));
		}
	}
	return B;
}

dmatrix_t *dmat_inverse(dmatrix_t *A)

{
	dmatrix_t *B = dmat_alloc(A->l, A->c);
	dmat_scalar_mult(B, dmat_transpose(cofactor(A)), 1.0 / determinant(A));
	return B;
}

dmatrix_t *to_homogeneous(dmatrix_t *A, double l)

{
	int i, j;
	dmatrix_t *B;

	if ((*A).l <= 0 || (*A).c <= 0) {
		nrerror("MATRIX.cpp: erroneous matrix size");
	}

	if (A->c == 1) {
		B = dmat_alloc(A->l + 1, 1);
		for (i = 0; i < A->l; i++) {
			B->m[i][0] = A->m[i][0];
		}
		B->m[B->l - 1][0] = l;
	}
	else if (A->l == 1) {
		B = dmat_alloc(1, A->l + 1);
		for (i = 0; i < A->c; i++) {
			B->m[0][i] = A->m[0][i];
		}
		B->m[0][B->c - 1] = l;
	}
	else {
		B = dmat_alloc(A->l + 1, A->c + 1);
		B = dmat_init(B, 0);

		for (i = 0; i < A->l; i ++) {
			for (j = 0; j < A->c; j ++) {
				B->m[i][j] = A->m[i][j];
			}
		}
		B->m[B->l - 1][B->c - 1] = l;
	}
	return B;
}

dmatrix_t *from_homogeneous(dmatrix_t *A)

{
	int i, j;
	dmatrix_t *B;

	if ((*A).l < 1 || (*A).c < 1) {
		nrerror("MATRIX.cpp: erroneous matrix size");
	}

	if (A->c == 1) {
		B = dmat_alloc(A->l - 1, 1);
		for (i = 0; i < B->l; i++) {
			B->m[i][0] = A->m[i][0];
		}
	}
	else if (A->l == 1) {
		B = dmat_alloc(1, A->c - 1);
		for (i = 0; i < B->c; i++) {
			B->m[0][i] = A->m[1][i];
		}
	}
	else {
		B = dmat_alloc(A->l - 1, A->c - 1);
		for (i = 0; i < B->l; i++) {
			for (j = 0; j < A->c; j++) {
				B->m[i][j] = A->m[i][j];
			}
		}
	}
	return B;
}

int main()
{
	dmatrix_t *matrix1 = dmat_alloc(4, 4);
	dmatrix_t *matrix2 = dmat_alloc(4, 4);
	dmatrix_t *result = dmat_alloc(4, 4);
	dmatrix_t *result1;
	dmatrix_t *tohomo;
	char option;
	double dot;
	int i, j;

	for (i = 0; i<4; i++)
	{
		for (j = 0; j < 4; j++)
		{
			matrix1->m[i][j] = rand() % 5;
		}
	}
	printf("\n\tMatrix1 is:\n");
	write_dmatrix(matrix1);

	for (i = 0; i<4; i++)
	{
		for (j = 0; j < 4; j++)
		{
			matrix2->m[i][j] = rand() % 5;
		}
	}
	printf("\n\tMatrix2 is:\n");
	write_dmatrix(matrix2);

	//switch option list start

	printf("\n\tA. Addition\n");
	printf("\n\tB. Subtraction\n");
	printf("\n\tC. Scalar Multiplication\n");
	printf("\n\tD. Multiplication\n");
	printf("\n\tE. Transpose\n");
	printf("\n\tF. Ddot Product\n");
	printf("\n\tG. To & From Homogeneous\n");

	printf("\n\tQ. Quit\n");
	printf("Your option is:");
	cin >> option;

	//switch option list end

	//switch statement

	while (option != 'q')
	{
		switch (option) {
			case 'a': {
				dmat_add(result, matrix1, matrix2);

				printf("\n\tAddition is:\n");
				write_dmatrix(result);
				break;
			}
			case 'b': {
				dmat_sub(result, matrix1, matrix2);

				printf("\n\tSubtraction is:\n");
				write_dmatrix(result);
				break;
			}
			case 'c': {
				dmat_scalar_mult(result, matrix1, 2);

				printf("\n\tMatrix1*2 is:\n");
				write_dmatrix(result);
				break;
			}
			case 'd': {
				result1 = dmat_mult(matrix1, matrix2);

				printf("\n\tMutiplication is:\n");
				write_dmatrix(result1);
				break;
			}
			case 'e': {
				result = dmat_transpose(matrix1);

				printf("\n\tTranspose of Matrix1 is:\n");
				write_dmatrix(result);
				break;
			}
			case 'f': {
				dot = ddot_product(matrix1, matrix2);

				printf("\n\tDdot Product is:\n");
				printf("%lf", dot);
				break;
			}
			case 'g': {
				double num = 3;
				tohomo = to_homogeneous(matrix1, num);
				printf("\n\tMatrix1 to homo is:\n");
				write_dmatrix(tohomo);

				result1 = from_homogeneous(tohomo);
				printf("\n\tHomo_matrix1 from homo is:\n");
				write_dmatrix(result1);

				break;
			}
			case 'q': {
				break;
			}
		}
		printf("Your option is:");
		cin >> option;

	}


}