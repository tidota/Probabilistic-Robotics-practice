// env_xyz.hpp
//
// this file defines the motion/measurement model in addition to the roobt's state.
//
// Ver.0.0.0-140927: derived from EKF simulator
// Ver.0.2.0-140929: implementation of GL
// Ver.0.3.0-140930: updates and bug fixes, inflation of noise, macros implemented
// Ver.0.4.0-141001: initialization of the first estimation implemented
// Ver.0.4.1-141003: inflation of noises as variables
// Ver.0.4.2-141003: minor updates in the algorithm to speed up (e.g., constant of log(2pi))
// Ver.0.4.3-141003: inflation of noises as constants (variable version was slow)
// Ver.0.4.4-141003: minor updates in the algorithm to speed up (eliminated function calls)
// Ver.0.5.0-141005: odometry implemented
// Ver.0.5.1-141005: GSL random number generator implemented

#ifndef _ENV_XYZ_HPP
#define _ENV_XYZ_HPP

#include "MV3.hpp"
#include <cmath>

// =====================================================
//              Random number generation
// =====================================================
//***define either one to choose standard library or GNU scientific library
//#define RAND_STD
#define RAND_GSL

#ifdef RAND_STD
#include <chrono>
#include <random>
#endif

#ifdef RAND_GSL
#include <time.h>
#include <gsl/gsl_randist.h>
#endif

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
#define SIGR 30.0

// Odometry Noise (mm per time step)
#define SIGO 5.0*DELTAT

// Floating Noise (mm per time step)
#define SIGF 50.0*DELTAT

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
double getMotionModelP(Vec3T next, Vec3T u, Vec3T prev, bool en_float);

// === real measurement === //
// location of sensor
Vec3T getLocSens(int indx_sens);
// to get measurements given the location of the robot
void getRealMeasure(Vec3T loc, double* r, int n_sens);
double getBeacModelP(Vec3T loc, double* r, int n_sens);

// === real mesurement (odometry) === //
// to get odometry
Vec3T getRealOdometry();
double getOdomModelP(Vec3T next, Vec3T prev);

// === Grid Localization === //
// grids[x][y][z]
//    each grid contains a probability that a true location is in it.
// parameters:
//    interval of grids (mm)
#define Unit_GRID 800.0
//    number of grids
#define N_GRIDS_X 12
#define N_GRIDS_Y 12
#define N_GRIDS_Z 12
//    location of the origin grid (grid[0][0][0])
#define ORG_GRID_X (-N_GRIDS_X/2*Unit_GRID)
#define ORG_GRID_Y (-N_GRIDS_Y/2*Unit_GRID)
#define ORG_GRID_Z (-N_GRIDS_Z/2*Unit_GRID)

// inflation of noises
// for motion model
#define INF_M (Unit_GRID*1.1)
// for measurement model (beacon)
#define INF_R (Unit_GRID*1.1)
// for odometry
#define INF_O (Unit_GRID*1.1)

// Grid Localization class
// grids containing probabilities
// this provides actual localization routine and output data
class RESULTPARAMS{
public:
	Vec3T mu; Mat33 Sigma;
};
class GL_MAIN{
private:
	double ***grids;
	double ***grids_old;
	double ***grids_buf;
	
public:
	GL_MAIN();
	~GL_MAIN();
	
	// to set an exact location to the estimation
	void setpos(Vec3T loc);
	
	// main algorithm for GL
	void localize(Vec3T u, double *z, int n_sens, bool en_cor, Vec3T odm, bool en_odm, bool en_float);
	
	// to get results of estimation
	RESULTPARAMS getResults();
};


#endif //_ENV_XYZ_HPP

