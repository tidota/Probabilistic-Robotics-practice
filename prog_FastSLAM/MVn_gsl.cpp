// MVn_gsl.cpp

#include "MVn_gsl.hpp"

void MatVec(double *result, gsl_matrix *m, double *vec){
	int n = (m->size1 < m->size2)? m->size1: m->size2;
	for(int i = 0; i < n; i++){
		result[i] = 0;
		for(int j = 0; j < n; j++){
			result[i] += gsl_matrix_get(m,i,j)*vec[j];
		}
	}
}

void MatMat(gsl_matrix *result, gsl_matrix *m1, gsl_matrix *m2){
	int n1 = (m1->size1 < m1->size2)? m1->size1: m1->size2;
	int n2 = (m2->size1 < m2->size2)? m2->size1: m2->size2;
	int n = (n1 < n2)? n1: n2;
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			double buff = 0;
			for(int k = 0; k < n; k++){
				buff += gsl_matrix_get(m1,i,k)*gsl_matrix_get(m2,k,j);
			}
			gsl_matrix_set(result,i,j,buff);
		}
	}
}

void VecMat(double *result, double *vec, gsl_matrix *m){
	int n = (m->size1 < m->size2)? m->size1: m->size2;
	for(int i = 0; i < n; i++){
		result[i] = 0;
		for(int j = 0; j < n; j++){
			result[i] += vec[j] * gsl_matrix_get(m,j,i);
		}
	}
}

void inv(gsl_matrix *result, gsl_matrix *m){
	int n = (m->size1 < m->size2)? m->size1: m->size2;
	int s;
	gsl_matrix *buff = gsl_matrix_alloc(n,n);
	gsl_matrix_memcpy(buff,m);
	gsl_permutation * perm = gsl_permutation_alloc (n);
	
	// Make LU decomposition of matrix m
	gsl_linalg_LU_decomp (buff, perm, &s);
	// Invert the matrix m
	gsl_linalg_LU_invert (buff, perm, result);
	
	gsl_matrix_free(buff);
	gsl_permutation_free(perm);
}

double det(gsl_matrix *m){
	int n = (m->size1 < m->size2)? m->size1: m->size2;
	if(rankm(m) < n) return 0;
	
	int s;
	gsl_matrix *buff = gsl_matrix_alloc(n,n);
	gsl_matrix_memcpy(buff,m);
	gsl_permutation * perm = gsl_permutation_alloc (n);
	
	// Make LU decomposition of matrix m
	gsl_linalg_LU_decomp (buff, perm, &s);
	// calculate determinant
	double det = gsl_linalg_LU_det(buff,s);
	
	gsl_matrix_free(buff);
	gsl_permutation_free(perm);
	
	return det;
}

int rankm(gsl_matrix *m){
	int small_dim = (m->size1 < m->size2)? m->size1: m->size2;
	gsl_matrix *U = gsl_matrix_alloc(m->size1,m->size1);
	gsl_matrix *V = gsl_matrix_alloc(m->size2,m->size2);
	gsl_vector *S = gsl_vector_alloc(small_dim);
	gsl_vector *work = gsl_vector_alloc(small_dim);
	
	gsl_matrix_memcpy(U,m);
	
	gsl_linalg_SV_decomp(U,V,S,work);
	int rank = 0;
	for(int ip = 0; ip < small_dim; ip++){
		if(isnan(gsl_vector_get(S,ip))){
			rank = 0;
			break;
		}else if(gsl_vector_get(S,ip) != 0){
			rank++;
		}
	}
	
	gsl_matrix_free(U);
	gsl_matrix_free(V);
	gsl_vector_free(S);
	gsl_vector_free(work);
	
	return rank;
}

void printMat(gsl_matrix *m){
	for(int i = 0; i < (int)m->size1; i++){
		for(int j = 0; j < (int)m->size2; j++){
			cout << fixed << setw(5) << setprecision(2) << gsl_matrix_get(m,i,j) << " ";
		}
		cout << endl;
	}
}
