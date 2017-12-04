// MV3.cpp
// 
// This code provides the definitions of matrix3X3 and vector, and their calculations.
// The general information is available in hpp file.
//
// Ver.0.0.0-140909: file created
// Ver.0.9.0-140909: functions implemented
// Ver.1.0.0-140909: inverse and multiplication debugged
// Ver.1.1.0-140910: norm() functions added

#include "MV3.hpp"

void initM33(double m1[][3],double val){
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			m1[i][j] = val;	
}
void cpM33(double m1[][3], double m2[][3]){
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			m1[i][j] = m2[i][j];
}
void cpV3(double* v1, double* v2){
	for(int i = 0; i < 3; i++)
		v1[i] = v2[i];
}

//=== matrix ===//
void Mat33::operator=(Mat33 m2){
	cpM33(val,m2.val);
}
Mat33 Mat33::operator+(Mat33 m2){
	Mat33 buff;
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			buff.val[i][j] = val[i][j] + m2.val[i][j];
	return buff;
}
Mat33 Mat33::operator-(Mat33 m2){
	Mat33 buff;
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			buff.val[i][j] = val[i][j] - m2.val[i][j];
	return buff;
}
Mat33 Mat33::operator*(double scl){
	Mat33 buff;
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			buff.val[i][j] = val[i][j] * scl;
	return buff;
}
Mat33 Mat33::operator/(double scl){
	Mat33 buff;
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			buff.val[i][j] = val[i][j] / scl;
	return buff;
}
Mat33 Mat33::operator*(Mat33 m2){
	Mat33 buff;
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++){
			buff.val[i][j] = 0;
			for(int k = 0; k < 3; k++)
				buff.val[i][j] += val[i][k] * m2.val[k][j];
		}
	return buff;
}
Vec3T Mat33::operator*(Vec3T v){
	Vec3T vec;
	for(int i = 0; i < 3; i++){
		vec.val[i] = 0;
		for(int j = 0; j < 3; j++)
			vec.val[i] += val[i][j] * v.val[j];
	}
	return vec;
}
double Mat33::subdet(int r, int c){
	double a[2][2];
	int ia = 0; int ja = 0;
	int r_updated = 0;
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++)
			if(i != r && j != c){
				a[ia][ja] = val[i][j];
				ja++;
				r_updated = 1;
			}
		if(r_updated){
			ia++;
			ja=0;
			r_updated = 0;
		}
	}
	return a[0][0]*a[1][1]-a[1][0]*a[0][1];
}
double Mat33::det(){
	double det_buff = 0;
	double sig = 1;
	for(int i = 0; i < 3; i++){
		det_buff += sig * val[0][i] * subdet(0,i);
		sig *= -1;
	}

	return det_buff;
}
Mat33 Mat33::inv(){
	Mat33 buff;
	initM33(buff.val,0);
	
	double det_buff = det();

	if(det_buff == 0){
		cout << "WARNING: det = 0" << endl;
		for(int i = 0; i < 3; i++)
			for(int j = 0; j < 3; j++)
				buff.val[i][j] = 0;	
	}else{
		for(int i = 0; i < 3; i++)
			for(int j = 0; j < 3; j++){
				buff.val[i][j] = subdet(j,i)/det_buff;
				if((i+j)%2)
					buff.val[i][j] *= -1;
			}
	}

	return buff;
}
Mat33 Mat33::trns(){
	Mat33 buff;
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			buff.val[i][j] = val[j][i];
	return buff;
}
void Mat33::print(){
	cout << "==== 3x3 matrix ====" << endl;
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++)
			cout << setw(10) << setprecision(2) << val[i][j] << "\t";
		cout << endl;
	}
	cout << endl;
}
Mat33 operator*(double scl, Mat33 m){
	Mat33 buff;
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			buff.val[i][j] = scl * m.val[i][j];	
	return buff;
}

//=== horizontal vector ===//
void Vec3::operator=(Vec3 v2){
	cpV3(val,v2.val);
}
Vec3 Vec3::operator+(Vec3 v2){
	Vec3 buff;
	for(int i = 0; i < 3; i++)
		buff.val[i] = val[i] + v2.val[i];
	return buff;
}
Vec3 Vec3::operator-(Vec3 v2){
	Vec3 buff;
	for(int i = 0; i < 3; i++)
		buff.val[i] = val[i] - v2.val[i];
	return buff;
}
Vec3 Vec3::operator*(double scl){
	Vec3 buff;
	for(int i = 0; i < 3; i++)
		buff.val[i] = val[i]*scl;
	return buff;
}
Vec3 Vec3::operator/(double scl){
	Vec3 buff;
	for(int i = 0; i < 3; i++)
		buff.val[i] = val[i]/scl;
	return buff;
}
Vec3 Vec3::operator*(Mat33 m){
	Vec3 buff;
	for(int i = 0; i < 3; i++){
		buff.val[i] = 0;
		for(int j = 0; j < 3; j++)
			buff.val[i] += val[j]*m.val[j][i];
	}
	return buff;
}
double Vec3::operator*(Vec3T v2){
	double buff = 0;
	for(int i = 0; i < 3; i++)
		buff += val[i]*v2.val[i];	
	return buff;
}
Vec3T Vec3::trns(){
	Vec3T buff;
	for(int i = 0; i < 3; i++)
		buff.val[i] = val[i];
	return buff;
}
double Vec3::norm(){
	return sqrt(val[0]*val[0]+val[1]*val[1]+val[2]*val[2]);
}
void Vec3::print(){
	cout << "==== horizontal vector ====" << endl;
	for(int i = 0; i < 3; i++)
		cout << setw(10) << setprecision(2) << val[i] << "\t";
	cout << endl << endl;
}
Vec3 operator*(double scl, Vec3 v){
	Vec3 buff;
	for(int i = 0; i < 3; i++)
		buff.val[i] = scl*v.val[i];
	return buff;
}

//=== vertical vector ===//
void Vec3T::operator=(Vec3T v2){
	cpV3(val,v2.val);
}
Vec3T Vec3T::operator+(Vec3T v2){
	Vec3T buff;
	for(int i = 0; i < 3; i++)
		buff.val[i] = val[i] + v2.val[i];
	return buff;
}
Vec3T Vec3T::operator-(Vec3T v2){
	Vec3T buff;
	for(int i = 0; i < 3; i++)
		buff.val[i] = val[i] - v2.val[i];
	return buff;
}
Vec3T Vec3T::operator*(double scl){
	Vec3T buff;
	for(int i = 0; i < 3; i++)
		buff.val[i] = val[i]*scl;
	return buff;
}
Vec3T Vec3T::operator/(double scl){
	Vec3T buff;
	for(int i = 0; i < 3; i++)
		buff.val[i] = val[i]/scl;
	return buff;
}
Mat33 Vec3T::operator*(Vec3 v2){
	Mat33 buff;
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			buff.val[i][j] = val[i] * v2.val[j];
	return buff;
}
Vec3 Vec3T::trns(){
	Vec3 buff;
	for(int i = 0; i < 3; i++)
		buff.val[i] = val[i];
	return buff;
}
double Vec3T::norm(){
	return sqrt(val[0]*val[0]+val[1]*val[1]+val[2]*val[2]);
}
void Vec3T::print(){
	cout << "==== vertical vector ====" << endl;
	for(int i = 0; i < 3; i++)
		cout << setw(10) << setprecision(2) << val[i] << endl;
	cout << endl << endl;
}
Vec3T operator*(double scl, Vec3T v){
	Vec3T buff;
	for(int i = 0; i < 3; i++)
		buff.val[i] = scl * v.val[i];
	return buff;
}
