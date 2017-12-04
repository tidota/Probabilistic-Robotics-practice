// MVn_gsl.hpp

// to compile:
// $ g++ <source file> -lgsl -lm -lgslcblas

// general infor of gsl matrix: https://www.gnu.org/software/gsl/manual/html_node/Matrices.html
// decomposition sample: http://www.macapp.net/pmwiki/pmwiki.php?n=Main.InvertMatrix
// assignment sample: https://www.gnu.org/software/gsl/manual/html_node/Example-programs-for-matrices.html

// Ver. 0.5.0-141207: file created and codes are derived from testing codes of gsl features
// Ver. 0.6.0-141208: det() added
// Ver. 0.7.0-141210: rankm() and other features have been added

/* sample of implementation
	double vec_result[3];
	double vec[3];
	gsl_matrix *Mat = gsl_matrix_alloc(3,3);
	gsl_matrix_set(Mat,0,0,0);gsl_matrix_set(Mat,0,1,9);gsl_matrix_set(Mat,0,2,0);
	gsl_matrix_set(Mat,1,0,0);gsl_matrix_set(Mat,1,1,0);gsl_matrix_set(Mat,1,2,8);
	gsl_matrix_set(Mat,2,0,7);gsl_matrix_set(Mat,2,1,0);gsl_matrix_set(Mat,2,2,0);
	vec[0] = 1; vec[1] = 2; vec[2] = 3;
	
	cout << "test Mat*Vec" << endl;
	MatVec(vec_result,Mat,vec);
	cout << "result: " << vec_result[0] << " " << vec_result[1] << " " << vec_result[2] << endl;
	
	cout << "test Vec*Mat" << endl;
	VecMat(vec_result,vec,Mat);
	cout << "result: " << vec_result[0] << " " << vec_result[1] << " " << vec_result[2] << endl;
	
	gsl_matrix_free(Mat);
	
	// Define the dimension n of the matrix
	// and the signum s (for LU decomposition)
	int n = 3;

	// Define all the used matrices
	gsl_matrix * m = gsl_matrix_alloc (n, n);
	gsl_matrix * inverse = gsl_matrix_alloc (n, n);
	gsl_matrix * result = gsl_matrix_alloc(n, n);

	// Fill the matrix m
	gsl_matrix_set(m,0,0,9.5);
	gsl_matrix_set(m,0,1,1.5);
	gsl_matrix_set(m,0,2,2.5);
	gsl_matrix_set(m,1,0,3.5);
	gsl_matrix_set(m,1,1,7.5);
	gsl_matrix_set(m,1,2,4.5);
	gsl_matrix_set(m,2,0,6.5);
	gsl_matrix_set(m,2,1,8.5);
	gsl_matrix_set(m,2,2,5.5);

	cout << "matrix: " << endl;
	printMat(m);
	
	inv(inverse,m);
	
	cout << "inverse: " << endl;
	printMat(inverse);
	
	// operation
	MatMat(result,inverse,m);
	
	// display
	cout << "result: " << endl;
	printMat(result);
	
	// memory release
	gsl_matrix_free(m);
	gsl_matrix_free(inverse);
	gsl_matrix_free(result);

*/


#ifndef _MVN_GSL_HPP
#define _MVN_GLS_HPP

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

// operation of matrix*vector
void MatVec(double *result, gsl_matrix *m, double *vec);

// operation of matrix*matrix
void MatMat(gsl_matrix *result, gsl_matrix *m1, gsl_matrix *m2);

// operation of vector*matrix
void VecMat(double *result, double *vec, gsl_matrix *m);

// generate an inverse of matrix m
// the inverse is stored to "result"
// result needs to have assinged memory space
void inv(gsl_matrix *result, gsl_matrix *m);

// calculate the determinant of matrix m
double det(gsl_matrix *m);

// calculate the rank of matrix m
int rankm(gsl_matrix *m);

// display matrix
void printMat(gsl_matrix *m);



#endif //_MVN_GSL_HPP
