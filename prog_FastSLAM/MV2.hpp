// MV2.hpp
// 
// This code is derived from MV3.hpp (Ver.1.1.0-140910)
// It provides the definitions of matrix2X2 and vector, and their calculations.
//
// Ver.0.0.0-141202: file derived from MV3.hpp
// Ver.0.1.0-141202: definition of matrix and vector renamed so as to indicate they deal with 2D
// Ver.0.1.1-141210: inverse debugged

#ifndef _MV2_HPP
#define _MV2_HPP

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

class Mat22; // 2x2 matrix
class Vec2;  // horizontal vector with 2 params
class Vec2T; // vertical vector with 2 params

//=== matrix ===//
class Mat22{
public:
	double val[2][2];
	void operator=(Mat22 m2);
	Mat22 operator+(Mat22 m2);
	Mat22 operator-(Mat22 m2);
	Mat22 operator*(double scl);
	Mat22 operator/(double scl);
	Mat22 operator*(Mat22 m2);
	Vec2T operator*(Vec2T v);
	double subdet(int r, int c);
	double det();
	Mat22 inv();
	Mat22 trns();
	void print();
};
Mat22 operator*(double scl, Mat22 m);

//=== horizontal vector ===//
class Vec2{
public:
	double val[3];
	void operator=(Vec2 v2);
	Vec2 operator+(Vec2 v2);
	Vec2 operator-(Vec2 v2);
	Vec2 operator*(double scl);
	Vec2 operator/(double scl);
	Vec2 operator*(Mat22 m);
	double operator*(Vec2T v2);
	Vec2T trns();
	double norm();
	void print();
};
Vec2 operator*(double scl, Vec2 v);

//=== vertical vector ===//
class Vec2T{
public:
	double val[3];
	void operator=(Vec2T v2);
	Vec2T operator+(Vec2T v2);
	Vec2T operator-(Vec2T v2);
	Vec2T operator*(double scl);
	Vec2T operator/(double scl);
	Mat22 operator*(Vec2 v2);
	Vec2 trns();
	double norm();
	void print();
};
Vec2T operator*(double scl, Vec2T v);

	
#endif //_MV2_HPP
