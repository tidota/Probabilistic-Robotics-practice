// env_xyz.hpp
//
// this file defines the motion/measurement model in addition to the roobt's state.
//
// Ver.0.0.0-140909: file created
// Ver.0.1.0-140910: functions added
// Ver.0.9.0-140911: reorganized setting section
// Ver.1.0.0-140911: option for enabling the correction step and # of sensors
// Ver.1.0.1-140916: optional floating effect added
// Ver.1.1.0-140919: odometry added
// Ver.1.1.1-140920: motion model updated (SIGM added to matrix M)
// Ver.1.1.2-140921: renamed SIGM to SIGF
//

#ifndef _ENV_XYZ_HPP
#define _ENV_XYZ_HPP

#include "MV3.hpp"
#include <cmath>
#include <chrono>
#include <random>

// ====================================================
//               SETTINGS
// ====================================================
// Delta t (sec)
#define DELTAT 1.0

// Control Noise
#define ALPHA1 0.2
#define ALPHA2 0.2
#define ALPHA3 0.2

// Measurement Noise
#define SIGR 10

// Odometry Noise (mm per time step)
#define SIGO 5*DELTAT

// Floating Noise (mm per time step)
#define SIGF 50*DELTAT

// maximum number of sensors (this corresponds m[][3] in env_xyz.cpp) 
#define MAX_N_SENS 8
// =====================================================
// =====================================================


#ifndef M_PI /* c++11 seems to eliminate M_PI */
#define M_PI 3.14159265358979
#endif

// === initialization of random number seriese === //
void initRand();

// === real world === //
Vec3T getRealNextState(Vec3T prev, Vec3T u, bool en_float);

// === real measurement === //
// location of sensor
Vec3T getLocSens(int indx_sens);
// to get measurements given the location of the robot
void getRealMeasure(Vec3T loc, double* r, int n_sens);

// === real mesurement (odometry) === //
// to get odometry
Vec3T getRealOdometry();

// === EKF === //
class EKFPARAMS{
public:
	Vec3T mu; // mean of the distribution
	Mat33 Sigma; // covariance
	double pz; // belief for sensor data (as output from the algorithm)
};
// function for EKF algorithm
//  mu: current estimated position, Sigma: covariances
//  u: input, z: list of measurement, n_sens: # of sensors
//  odm: odometry in x, y, and z
//  en_cor: enable the correction step if it is true
EKFPARAMS EKFlocalization(Vec3T mu, Mat33 Sigma, Vec3T u, double *z, int n_sens, bool en_cor, Vec3T odm, bool en_odm);

#endif //_ENV_XYZ_HPP

