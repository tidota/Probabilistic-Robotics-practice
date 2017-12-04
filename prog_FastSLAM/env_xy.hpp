// env_xy.hpp
//
// This file was derived from env_xyz.hpp (Ver.1.1.2-140921)
// It defines the motion/measurement model in addition to the roobt's state.
//
// Ver.0.0.0-141202: file derived from env_xyz.hpp
// Ver.0.1.0-141202: some variables renamed to be compatible with MV2.cpp/hpp (EKF codes not removed yet)
// Ver.0.1.1-141202: debugged and compiled correctly
// Ver.0.2.0-141208: direction of sensor reading added (now sensor reading is 2D)
// Ver.0.9.0-141210: mostly done
//

#ifndef _ENV_XY_HPP
#define _ENV_XY_HPP

#include "MV2.hpp"
#include "MVn_gsl.hpp"
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

// Measurement Noise
#define SIGR 10
#define SIGR2 (3.0/180.0*M_PI)

// Odometry Noise (mm per time step)
#define SIGO 5*DELTAT

// Floating Noise (mm per time step)
#define SIGF 50*DELTAT

// maximum number of sensors (this corresponds m[][2] in env_xyz.cpp) 
#define MAX_N_SENS 8
// =====================================================
// =====================================================


#ifndef M_PI /* c++11 seems to eliminate M_PI */
#define M_PI 3.14159265358979
#endif

// === initialization of random number seriese === //
void initRand();

// === real world === //
Vec2T getRealNextState(Vec2T prev, Vec2T u, bool en_float);

// === real measurement === //
// location of sensor
Vec2T getLocSens(int indx_sens);
// to get measurements given the location of the robot
void getRealMeasure(Vec2T loc, Vec2T* r, int n_sens);
//void getRealMeasure(Vec2T loc, double* r, int n_sens);

// === real mesurement (odometry) === //
// to get odometry
Vec2T getRealOdometry();

// === FastSLAM1.0 === //
// # of particles
#define N_PART 10000

// default importance weight
#define P_0 1

class FEAT{      // feature class
public:
	Vec2T mu;     // estimated position
	Mat22 Sigma;  // Gaussian matrix
};

class PARTICLE{ // contains robot's pose and a list of features
public:
	Vec2T pose; // robot's pose
	FEAT feat[MAX_N_SENS]; // estimation of features
	double omega; // importance weight
public:
	PARTICLE(Vec2T init_pose);
	PARTICLE();
};

//=== main estimation function ===
// Y is a list of particles
// n_part is # of particles
// z is a list of sensor readings
// n_sens is # of sensor readings
// u is input
// init is a flag to indicate this interation is the first
// * correspondence of sensors is supposed to be known
void estimate(PARTICLE *Y, Vec2T *z, int n_sens, Vec2T u, bool init, bool enable);

void initPart(PARTICLE *Y, Vec2T vec);
Vec2T getAve(PARTICLE *Y);

#endif //_ENV_XY_HPP

