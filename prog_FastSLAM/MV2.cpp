// MV2.cpp
// 
// This code is derived from MV2.cpp (Ver.1.1.0-140910)
// It provides the definitions of matrix2X2 and vector, and their calculations.
// The general information is available in hpp file.
// The version of this file is integrated with that of hpp file for simplicity.

#include "MV2.hpp"

void initM22(double m1[][2],double val){
	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 2; j++)
			m1[i][j] = val;	
}
void cpM22(double m1[][2], double m2[][2]){
	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 2; j++)
			m1[i][j] = m2[i][j];
}
void cpV2(double* v1, double* v2){
	for(int i = 0; i < 2; i++)
		v1[i] = v2[i];
}

//=== matrix ===//
void Mat22::operator=(Mat22 m2){
	cpM22(val,m2.val);
}
Mat22 Mat22::operator+(Mat22 m2){
	Mat22 buff;
	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 2; j++)
			buff.val[i][j] = val[i][j] + m2.val[i][j];
	return buff;
}
Mat22 Mat22::operator-(Mat22 m2){
	Mat22 buff;
	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 2; j++)
			buff.val[i][j] = val[i][j] - m2.val[i][j];
	return buff;
}
Mat22 Mat22::operator*(double scl){
	Mat22 buff;
	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 2; j++)
			buff.val[i][j] = val[i][j] * scl;
	return buff;
}
Mat22 Mat22::operator/(double scl){
	Mat22 buff;
	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 2; j++)
			buff.val[i][j] = val[i][j] / scl;
	return buff;
}
Mat22 Mat22::operator*(Mat22 m2){
	Mat22 buff;
	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 2; j++){
			buff.val[i][j] = 0;
			for(int k = 0; k < 2; k++)
				buff.val[i][j] += val[i][k] * m2.val[k][j];
		}
	return buff;
}
Vec2T Mat22::operator*(Vec2T v){
	Vec2T vec;
	for(int i = 0; i < 2; i++){
		vec.val[i] = 0;
		for(int j = 0; j < 2; j++)
			vec.val[i] += val[i][j] * v.val[j];
	}
	return vec;
}
double Mat22::det(){
	return val[0][0]*val[1][1]-val[1][0]*val[0][1];
}
Mat22 Mat22::inv(){
	Mat22 buff;
	initM22(buff.val,0);
	
	double det_buff = det();

	if(det_buff == 0){
		cout << "WARNING: det = 0" << endl;
		for(int i = 0; i < 2; i++)
			for(int j = 0; j < 2; j++)
				buff.val[i][j] = 0;	
	}else{
		buff.val[0][0] = val[1][1];
		buff.val[0][1] = -1*val[0][1];
		buff.val[1][0] = -1*val[1][0];
		buff.val[1][1] = val[0][0];
		buff = buff / det_buff;
	}

	return buff;
}
Mat22 Mat22::trns(){
	Mat22 buff;
	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 2; j++)
			buff.val[i][j] = val[j][i];
	return buff;
}
void Mat22::print(){
	cout << "==== 2x2 matrix ====" << endl;
	for(int i = 0; i < 2; i++){
		for(int j = 0; j < 2; j++)
			cout << setw(10) << setprecision(2) << val[i][j] << "\t";
		cout << endl;
	}
	cout << endl;
}
Mat22 operator*(double scl, Mat22 m){
	Mat22 buff;
	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 2; j++)
			buff.val[i][j] = scl * m.val[i][j];	
	return buff;
}

//=== horizontal vector ===//
void Vec2::operator=(Vec2 v2){
	cpV2(val,v2.val);
}
Vec2 Vec2::operator+(Vec2 v2){
	Vec2 buff;
	for(int i = 0; i < 2; i++)
		buff.val[i] = val[i] + v2.val[i];
	return buff;
}
Vec2 Vec2::operator-(Vec2 v2){
	Vec2 buff;
	for(int i = 0; i < 2; i++)
		buff.val[i] = val[i] - v2.val[i];
	return buff;
}
Vec2 Vec2::operator*(double scl){
	Vec2 buff;
	for(int i = 0; i < 2; i++)
		buff.val[i] = val[i]*scl;
	return buff;
}
Vec2 Vec2::operator/(double scl){
	Vec2 buff;
	for(int i = 0; i < 2; i++)
		buff.val[i] = val[i]/scl;
	return buff;
}
Vec2 Vec2::operator*(Mat22 m){
	Vec2 buff;
	for(int i = 0; i < 2; i++){
		buff.val[i] = 0;
		for(int j = 0; j < 2; j++)
			buff.val[i] += val[j]*m.val[j][i];
	}
	return buff;
}
double Vec2::operator*(Vec2T v2){
	double buff = 0;
	for(int i = 0; i < 2; i++)
		buff += val[i]*v2.val[i];	
	return buff;
}
Vec2T Vec2::trns(){
	Vec2T buff;
	for(int i = 0; i < 2; i++)
		buff.val[i] = val[i];
	return buff;
}
double Vec2::norm(){
	return sqrt(val[0]*val[0]+val[1]*val[1]);
}
void Vec2::print(){
	cout << "==== horizontal vector ====" << endl;
	for(int i = 0; i < 2; i++)
		cout << setw(10) << setprecision(2) << val[i] << "\t";
	cout << endl << endl;
}
Vec2 operator*(double scl, Vec2 v){
	Vec2 buff;
	for(int i = 0; i < 2; i++)
		buff.val[i] = scl*v.val[i];
	return buff;
}

//=== vertical vector ===//
void Vec2T::operator=(Vec2T v2){
	cpV2(val,v2.val);
}
Vec2T Vec2T::operator+(Vec2T v2){
	Vec2T buff;
	for(int i = 0; i < 2; i++)
		buff.val[i] = val[i] + v2.val[i];
	return buff;
}
Vec2T Vec2T::operator-(Vec2T v2){
	Vec2T buff;
	for(int i = 0; i < 2; i++)
		buff.val[i] = val[i] - v2.val[i];
	return buff;
}
Vec2T Vec2T::operator*(double scl){
	Vec2T buff;
	for(int i = 0; i < 2; i++)
		buff.val[i] = val[i]*scl;
	return buff;
}
Vec2T Vec2T::operator/(double scl){
	Vec2T buff;
	for(int i = 0; i < 2; i++)
		buff.val[i] = val[i]/scl;
	return buff;
}
Mat22 Vec2T::operator*(Vec2 v2){
	Mat22 buff;
	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 2; j++)
			buff.val[i][j] = val[i] * v2.val[j];
	return buff;
}
Vec2 Vec2T::trns(){
	Vec2 buff;
	for(int i = 0; i < 2; i++)
		buff.val[i] = val[i];
	return buff;
}
double Vec2T::norm(){
	return sqrt(val[0]*val[0]+val[1]*val[1]);
}
void Vec2T::print(){
	cout << "==== vertical vector ====" << endl;
	for(int i = 0; i < 2; i++)
		cout << setw(10) << setprecision(2) << val[i] << endl;
	cout << endl << endl;
}
Vec2T operator*(double scl, Vec2T v){
	Vec2T buff;
	for(int i = 0; i < 2; i++)
		buff.val[i] = scl * v.val[i];
	return buff;
}
