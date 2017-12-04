// env_xyz.cpp
//
// this file defines the motion/measurement model in addition to the roobt's state.
//
//
//

#include "env_xyz.hpp"

// === random number seriese === //
#ifdef RAND_STD
default_random_engine generator;
//normal_distribution<double> motion_noise(0.0,SIGF);
normal_distribution<double> measur_noise(0.0,SIGR*DELTAT);
normal_distribution<double> odom_noise(0.0,SIGO);
#endif
#ifdef RAND_GSL
gsl_rng *generator;
#endif

void initRand(){
#ifdef RAND_STD
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine g (seed);
	generator = g;
#endif
#ifdef RAND_GSL
	unsigned seed = time(NULL);
	generator = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(generator,seed);
#endif
}

// === real world === //
static Vec3T delta_mu; // motion at each time step (it is used for odometry)

Vec3T getRealNextState(Vec3T prev, Vec3T u, bool en_float){
	Vec3T a;

#ifdef RAND_STD
	normal_distribution<double> motion_noise_x(0.0,sqrt(ALPHA1)*abs(u.val[0]));
	normal_distribution<double> motion_noise_y(0.0,sqrt(ALPHA2)*abs(u.val[1]));
	normal_distribution<double> motion_noise_z(0.0,sqrt(ALPHA3)*abs(u.val[2]));
	a.val[0] = motion_noise_x(generator);
	a.val[1] = motion_noise_y(generator);
	a.val[2] = motion_noise_z(generator);
#endif
#ifdef RAND_GSL
	a.val[0] = gsl_ran_gaussian(generator,sqrt(ALPHA1)*abs(u.val[0]));
	a.val[1] = gsl_ran_gaussian(generator,sqrt(ALPHA1)*abs(u.val[1]));
	a.val[2] = gsl_ran_gaussian(generator,sqrt(ALPHA1)*abs(u.val[2]));
#endif

	if(en_float){ // floating effect is enabled
		Vec3T flt;
#ifdef RAND_STD
		normal_distribution<double> floating_noise(0.0,SIGF);
		flt.val[0] = floating_noise(generator);
		flt.val[1] = floating_noise(generator);
		flt.val[2] = floating_noise(generator);
#endif
#ifdef RAND_GSL
		flt.val[0] = gsl_ran_gaussian(generator,SIGF);
		flt.val[1] = gsl_ran_gaussian(generator,SIGF);
		flt.val[2] = gsl_ran_gaussian(generator,SIGF);
#endif
		delta_mu = (u+a)*DELTAT + flt;
	}else{
		delta_mu = (u+a)*DELTAT;
	}
	return prev + delta_mu;
}

// === real measurement === //
Vec3T m[MAX_N_SENS] = {\
	{{ 10000, 10000, 10000}},\
	{{-10000, 10000, 10000}},\
	{{-10000,-10000,10000}},\
	{{ 10000,-10000,10000}},\
	{{ 10000, 10000,-10000}},\
	{{-10000, 10000,-10000}},\
	{{-10000,-10000,-10000}},\
	{{ 10000,-10000,-10000}}};

// location of sensor
Vec3T getLocSens(int indx_sens){
	return m[indx_sens];
}
// to get measurements given the location of the robot
void getRealMeasure(Vec3T loc, double* r, int n_sens){
	for(int i = 0; i < n_sens; i++){
#ifdef RAND_STD
		r[i] = (m[i] - loc).norm() + measur_noise(generator);
#endif
#ifdef RAND_GSL
		r[i] = (m[i] - loc).norm() + gsl_ran_gaussian(generator,SIGR*DELTAT);
#endif
	}
}

// === real mesurement (odometry) === //
// to get odometry
Vec3T getRealOdometry(){
	Vec3T odm;
#ifdef RAND_STD
	odm.val[0] = delta_mu.val[0] + odom_noise(generator);
	odm.val[1] = delta_mu.val[1] + odom_noise(generator);
	odm.val[2] = delta_mu.val[2] + odom_noise(generator);
#endif
#ifdef RAND_GSL
	odm.val[0] = delta_mu.val[0] + gsl_ran_gaussian(generator,SIGO);
	odm.val[1] = delta_mu.val[1] + gsl_ran_gaussian(generator,SIGO);
	odm.val[2] = delta_mu.val[2] + gsl_ran_gaussian(generator,SIGO);
#endif
	return odm;
}

// Grid Localization class
GL_MAIN::GL_MAIN(){
	
	grids = new double**[N_GRIDS_X];
	for(int ix = 0; ix < N_GRIDS_X; ix++){
		grids[ix] = new double*[N_GRIDS_Y];
		for(int iy = 0; iy < N_GRIDS_Y; iy++){
			grids[ix][iy] = new double[N_GRIDS_Z];
			for(int iz = 0; iz < N_GRIDS_Z; iz++)
				grids[ix][iy][iz] = log(1.0 / N_GRIDS_X / N_GRIDS_Y / N_GRIDS_Z);
		}
	}
	
	grids_old = new double**[N_GRIDS_X];
	for(int ix = 0; ix < N_GRIDS_X; ix++){
		grids_old[ix] = new double*[N_GRIDS_Y];
		for(int iy = 0; iy < N_GRIDS_Y; iy++)
			grids_old[ix][iy] = new double[N_GRIDS_Z];
	}
	
	grids_buf = new double**[N_GRIDS_X];
	for(int ix = 0; ix < N_GRIDS_X; ix++){
		grids_buf[ix] = new double*[N_GRIDS_Y];
		for(int iy = 0; iy < N_GRIDS_Y; iy++)
			grids_buf[ix][iy] = new double[N_GRIDS_Z];
	}
}
GL_MAIN::~GL_MAIN(){
	for(int ix = 0; ix < N_GRIDS_X; ix++){
		for(int iy = 0; iy < N_GRIDS_Y; iy++){
			delete grids[ix][iy];
		}
		delete grids[ix];
	}
	delete grids;
	for(int ix = 0; ix < N_GRIDS_X; ix++){
		for(int iy = 0; iy < N_GRIDS_Y; iy++){
			delete grids_old[ix][iy];
		}
		delete grids_old[ix];
	}
	delete grids_old;
	for(int ix = 0; ix < N_GRIDS_X; ix++){
		for(int iy = 0; iy < N_GRIDS_Y; iy++){
			delete grids_buf[ix][iy];
		}
		delete grids_buf[ix];
	}
	delete grids_buf;
}

// set the initial estimation
//  loc has coordinates of x, y, and z
void GL_MAIN::setpos(Vec3T loc){
	// find the nearest grid
	int x = (loc.val[0] - ORG_GRID_X + Unit_GRID/2)/Unit_GRID;
	int y = (loc.val[1] - ORG_GRID_Y + Unit_GRID/2)/Unit_GRID;
	int z = (loc.val[2] - ORG_GRID_Z + Unit_GRID/2)/Unit_GRID;
	
	// initalize grids so that only the specified grid has 1.
	for(int ix = 0; ix < N_GRIDS_X; ix++)
	 for(int iy = 0; iy < N_GRIDS_Y; iy++)
	  for(int iz = 0; iz < N_GRIDS_Z; iz++)
	  	if(x == ix && y == iy && z == iz)
	  		grids[ix][iy][iz] = 1.0;
	  	else
	  		grids[iz][iy][iz] = 0.0;
	
}


// probability based on motion model
// refer to p.123 of Probabilistic Robotics
/*
#define getMotionModelLogP(logp, next, u, prev, en_float) do{ \
	Vec3T u_est;\
	u_est.val[0] = (next.val[0] - prev.val[0])/DELTAT;\
	u_est.val[1] = (next.val[1] - prev.val[1])/DELTAT;\
	u_est.val[2] = (next.val[2] - prev.val[2])/DELTAT;\
	\
	logp = 0;\
	double alpha[3] = {ALPHA1, ALPHA2, ALPHA3};\
	for(int i = 0; i < 3; i++){\
		double variance = alpha[i]*u.val[i]*u.val[i] + INF_M*INF_M;\
		logp += \
			-1.0/2*(u_est.val[i]-u.val[i])*(u_est.val[i]-u.val[i])/variance +\
			log(1.0/sqrt(2*M_PI*variance));\
	}\
}while(0)
*/
/*
#define getMotionModelLogP(logp, next, u, prev, en_float) do{ \
	Vec3T u_est = ((next) - (prev))/DELTAT;\
	\
	(logp) = 0;\
	const double alpha[3] = {ALPHA1, ALPHA2, ALPHA3};\
	for(int i = 0; i < 3; i++){\
		double variance = alpha[i]*u.val[i]*u.val[i] + INF_M*INF_M;\
		logp += \
			-1.0/2*(u_est.val[i]-(u).val[i])*(u_est.val[i]-(u).val[i])/variance\
			-(log(variance) + LOG_2PI)/2.0;\
	}\
}while(0)
*/

// log(2*pi)
#define LOG_2PI 1.8378770664093454835606594728112
#define getMotionModelLogP(logp, next, u, prev, en_float) do{ \
	Vec3T u_est;\
	u_est.val[0] = ((next).val[0] - (prev).val[0])/DELTAT;\
	u_est.val[1] = ((next).val[1] - (prev).val[1])/DELTAT;\
	u_est.val[2] = ((next).val[2] - (prev).val[2])/DELTAT;\
	\
	(logp) = 0;\
	const double alpha[3] = {ALPHA1, ALPHA2, ALPHA3};\
	for(int i = 0; i < 3; i++){\
		double variance = alpha[i]*u.val[i]*u.val[i] + INF_M*INF_M;\
		logp += \
			-1.0/2*(u_est.val[i]-(u).val[i])*(u_est.val[i]-(u).val[i])/variance\
			-(log(variance) + LOG_2PI)/2.0;\
	}\
}while(0)

// probability based on beacon measurements
// refer to p.179 in Probabilistic Robotics
#define getBeacModelLogP(logp, loc, r, n_sens) do{\
	(logp) = 0;\
	\
	for(int i = 0; i < (n_sens); i++){\
		double leng = ((loc) - m[i]).norm();\
		logp += \
			-1.0/2*(leng - (r)[i])*(leng - (r)[i])/(INF_R*INF_R)\
			-(log((INF_R*INF_R)) + LOG_2PI)/2.0;\
	}\
}while(0)

#define getOdomModelLogP(logp, next, prev, odm) do{\
	Vec3T odm_est;\
	odm_est.val[0] = ((next).val[0] - (prev).val[0])/DELTAT;\
	odm_est.val[1] = ((next).val[1] - (prev).val[1])/DELTAT;\
	odm_est.val[2] = ((next).val[2] - (prev).val[2])/DELTAT;\
	\
	(logp) = 0;\
	for(int i = 0; i < 3; i++){\
		(logp) += \
			-1.0/2*(odm_est.val[i]-(odm).val[i])*(odm_est.val[i]-(odm).val[i])/(INF_O*INF_O)\
			-(log(INF_O*INF_O) + LOG_2PI)/2.0;\
	}\
}while(0)


void GL_MAIN::localize(Vec3T u, double *z, int n_sens, bool en_cor, Vec3T odm, bool en_odm, bool en_float){
	// copy old beliefs
	for(int ix = 0; ix < N_GRIDS_X; ix++)
	 for(int iy = 0; iy < N_GRIDS_Y; iy++)
	  for(int iz = 0; iz < N_GRIDS_Z; iz++)
	   grids_old[ix][iy][iz] = grids[ix][iy][iz];
	
	
	double max_logp = 0; // maximum probability in grids
	bool f_max_logp_updated = false;
	for(int ix = 0; ix < N_GRIDS_X; ix++)
	 for(int iy = 0; iy < N_GRIDS_Y; iy++)
	  for(int iz = 0; iz < N_GRIDS_Z; iz++){
	  	// ============ motion model ============ //
	  	Vec3T next;
	  	next.val[0] = ORG_GRID_X + Unit_GRID * ix;
	  	next.val[1] = ORG_GRID_Y + Unit_GRID * iy;
	  	next.val[2] = ORG_GRID_Z + Unit_GRID * iz;
	  	// store lograrithm of probability for each grid
	  	double max_logp_lcl; // local maximum probability in the process for a grid
	  	bool f_max_logp_lcl_updated = false;
	  	double buff_logp;
	  	for(int jx = 0; jx < N_GRIDS_X; jx++)
	  	 for(int jy = 0; jy < N_GRIDS_Y; jy++)
	  	  for(int jz = 0; jz < N_GRIDS_Z; jz++){
	  	  	Vec3T prev;
	  	  	prev.val[0] = ORG_GRID_X + Unit_GRID * jx;
	  	  	prev.val[1] = ORG_GRID_Y + Unit_GRID * jy;
	  	  	prev.val[2] = ORG_GRID_Z + Unit_GRID * jz;
	  	  	
	  	  	getMotionModelLogP(buff_logp, next, u, prev, en_float);
	  	  	grids_buf[jx][jy][jz] = grids_old[jx][jy][jz] + buff_logp;
	  	  	if(f_max_logp_lcl_updated == false || max_logp_lcl < grids_buf[jx][jy][jz]){
	  	  		max_logp_lcl = grids_buf[jx][jy][jz];
	  	  		f_max_logp_lcl_updated = true;
	  	  	}
	  	  }
	  	
	  	// calculate sum divided by the maximum probability
	  	double sum_p = 0;
	  	for(int jx = 0; jx < N_GRIDS_X; jx++)
	  	 for(int jy = 0; jy < N_GRIDS_Y; jy++)
	  	  for(int jz = 0; jz < N_GRIDS_Z; jz++)
	  	  	sum_p += exp(grids_buf[jx][jy][jz] - max_logp_lcl);
	  	
	  	// calculate logarithm of the exact probability
	  	double logp = log(sum_p) + max_logp_lcl;
	  	
	  	// ============ mesurement model (beacons) ============ //
	  	if(en_cor){
	  		getBeacModelLogP(buff_logp, next, z, n_sens);
	  		logp += buff_logp;
	  	}
	  	
	  	// ============ mesurement model (odometry) ============ //
	  	if(en_odm){
	  		f_max_logp_lcl_updated = false;
	  		for(int jx = 0; jx < N_GRIDS_X; jx++)
	  		 for(int jy = 0; jy < N_GRIDS_Y; jy++)
	  		  for(int jz = 0; jz < N_GRIDS_Z; jz++){
	  		  	Vec3T prev;
	  		  	prev.val[0] = ORG_GRID_X + Unit_GRID * jx;
	  		  	prev.val[1] = ORG_GRID_Y + Unit_GRID * jy;
	  		  	prev.val[2] = ORG_GRID_Z + Unit_GRID * jz;
	  		  	getOdomModelLogP(buff_logp, next, prev, odm);
	  		  	grids_buf[jx][jy][jz] = grids_old[jx][jy][jz] + buff_logp;
	  		  	if(f_max_logp_lcl_updated == false || max_logp_lcl < grids_buf[jx][jy][jz]){
	  		  		max_logp_lcl = grids_buf[jx][jy][jz];
	  		  		f_max_logp_lcl_updated = true;
	  		  	}
	  		  }
	  		// calculate sum divided by the maximum probability
	  		sum_p = 0;
	  		for(int jx = 0; jx < N_GRIDS_X; jx++)
	  		 for(int jy = 0; jy < N_GRIDS_Y; jy++)
	  		  for(int jz = 0; jz < N_GRIDS_Z; jz++)
	  		  	sum_p += exp(grids_buf[jx][jy][jz] - max_logp_lcl);
	  		
	  		// calculate logarithm of the exact probability
	  		logp += log(sum_p) + max_logp_lcl;
	  	}
	
	  	grids[ix][iy][iz] = logp;
	  	
	  	// set maximum
	  	if(f_max_logp_updated == false || max_logp < logp){
	  		max_logp = logp;
	  		f_max_logp_updated = true;
	  	}
	  }
	
	// normalization
	double total_p = 0;
	for(int ix = 0; ix < N_GRIDS_X; ix++)
	 for(int iy = 0; iy < N_GRIDS_Y; iy++)
	  for(int iz = 0; iz < N_GRIDS_Z; iz++){
	  	grids[ix][iy][iz] -= max_logp;
	  	total_p += exp(grids[ix][iy][iz]);
	  }
	if(total_p != 0)
		for(int ix = 0; ix < N_GRIDS_X; ix++)
		 for(int iy = 0; iy < N_GRIDS_Y; iy++)
		  for(int iz = 0; iz < N_GRIDS_Z; iz++)
		  	grids[ix][iy][iz] -= log(total_p);
}
RESULTPARAMS GL_MAIN::getResults(){
	// calculate mean
	Vec3T mu = {{0,0,0}};
	
	for(int ix = 0; ix < N_GRIDS_X; ix++)
	 for(int iy = 0; iy < N_GRIDS_Y; iy++)
	  for(int iz = 0; iz < N_GRIDS_Z; iz++){
	  	mu.val[0] += exp(grids[ix][iy][iz])*(ORG_GRID_X+Unit_GRID*ix);
	   	mu.val[1] += exp(grids[ix][iy][iz])*(ORG_GRID_Y+Unit_GRID*iy);
	   	mu.val[2] += exp(grids[ix][iy][iz])*(ORG_GRID_Z+Unit_GRID*iz);
	  }
	
	// calculate covariances
	Mat33 Sig = {{{0,0,0},{0,0,0},{0,0,0}}};
	for(int ix = 0; ix < N_GRIDS_X; ix++)
	 for(int iy = 0; iy < N_GRIDS_Y; iy++)
	  for(int iz = 0; iz < N_GRIDS_Z; iz++){
		double buf[3];
		buf[0] = ORG_GRID_X+Unit_GRID*ix;
		buf[1] = ORG_GRID_Y+Unit_GRID*iy;
		buf[2] = ORG_GRID_Z+Unit_GRID*iz;
		for(int ir = 0; ir < 3; ir++)
			for(int ic = 0; ic < 3; ic++)
				Sig.val[ir][ic] += exp(grids[ix][iy][iz])*(buf[ir] - mu.val[ir])*(buf[ic] - mu.val[ic]);
	  }
	
	RESULTPARAMS results;
	results.mu = mu;
	results.Sigma = Sig;
	
	return results;
}



