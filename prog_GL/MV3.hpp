// MV3.hpp
// 
// This code provides the definitions of matrix3X3 and vector, and their calculations.
//
// Ver.0.0.0-140909: file created
// Ver.0.0.1-140909: temporary definitions
// Ver.0.9.0-140909: members in classes defined
// Ver.1.0.0-140909: function list completed at this point
// Ver.1.1.0-140910: norm() functions added

#ifndef _MV3_HPP
#define _MV3_HPP

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

class Mat33; // 3x3 matrix
class Vec3;  // horizontal vector with 3 params
class Vec3T; // vertical vector with 3 params

//=== matrix ===//
class Mat33{
public:
	double val[3][3];
	void operator=(Mat33 m2);
	Mat33 operator+(Mat33 m2);
	Mat33 operator-(Mat33 m2);
	Mat33 operator*(double scl);
	Mat33 operator/(double scl);
	Mat33 operator*(Mat33 m2);
	Vec3T operator*(Vec3T v);
	double subdet(int r, int c);
	double det();
	Mat33 inv();
	Mat33 trns();
	void print();
};
Mat33 operator*(double scl, Mat33 m);

//=== horizontal vector ===//
class Vec3{
public:
	double val[3];
	void operator=(Vec3 v2);
	Vec3 operator+(Vec3 v2);
	Vec3 operator-(Vec3 v2);
	Vec3 operator*(double scl);
	Vec3 operator/(double scl);
	Vec3 operator*(Mat33 m);
	double operator*(Vec3T v2);
	Vec3T trns();
	double norm();
	void print();
};
Vec3 operator*(double scl, Vec3 v);

//=== vertical vector ===//
class Vec3T{
public:
	double val[3];
	void operator=(Vec3T v2);
	Vec3T operator+(Vec3T v2);
	Vec3T operator-(Vec3T v2);
	Vec3T operator*(double scl);
	Vec3T operator/(double scl);
	Mat33 operator*(Vec3 v2);
	Vec3 trns();
	double norm();
	void print();
};
Vec3T operator*(double scl, Vec3T v);

	
#endif //_MV3_HPP
