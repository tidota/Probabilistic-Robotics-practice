// env_xyz.cpp
//
// this file defines the motion/measurement model in addition to the roobt's state.
//
// Ver.0.0.0-140910: file created
// Ver.0.9.0-140910: mostly done
// Ver.1.0.0-140911: correction step can be disabled by giving an argument
// Ver.1.1.0-140916: motion model updated
// Ver.1.1.1-140916: optional floating effect added
// Ver.1.2.0-140919: odometry added
// Ver.1.2.1-140920: motion model updated (SIGF added to matrix M)
// Ver.1.2.2-140921: renamed SIGM to SIGF
//
//

#include "env_xyz.hpp"

// === random number seriese === //
default_random_engine generator;
//normal_distribution<double> motion_noise(0.0,SIGF);
normal_distribution<double> measur_noise(0.0,SIGR*DELTAT);
normal_distribution<double> odom_noise(0.0,SIGO);
void initRand(){
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine g (seed);
	generator = g;
}

// === real world === //
static Vec3T delta_mu; // motion at each time step (it is used for odometry)

Vec3T getRealNextState(Vec3T prev, Vec3T u, bool en_float){
	normal_distribution<double> motion_noise_x(0.0,sqrt(ALPHA1)*abs(u.val[0]));
	normal_distribution<double> motion_noise_y(0.0,sqrt(ALPHA2)*abs(u.val[1]));
	normal_distribution<double> motion_noise_z(0.0,sqrt(ALPHA3)*abs(u.val[2]));
	Vec3T a;
	a.val[0] = motion_noise_x(generator);
	a.val[1] = motion_noise_y(generator);
	a.val[2] = motion_noise_z(generator);
	if(en_float){ // floating effect is enabled
		normal_distribution<double> floating_noise(0.0,SIGF);
		Vec3T flt;
		flt.val[0] = floating_noise(generator);
		flt.val[1] = floating_noise(generator);
		flt.val[2] = floating_noise(generator);
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
		r[i] = (m[i] - loc).norm() + measur_noise(generator);
	}
}

// === real mesurement (odometry) === //
// to get odometry
Vec3T getRealOdometry(){
	Vec3T odm;
	odm.val[0] = delta_mu.val[0] + odom_noise(generator);
	odm.val[1] = delta_mu.val[1] + odom_noise(generator);
	odm.val[2] = delta_mu.val[2] + odom_noise(generator);
	return odm;
}


// === EKF === //
EKFPARAMS EKFlocalization(Vec3T mu, Mat33 Sigma, Vec3T u, double *z, int n_sens, bool en_cor, Vec3T odm, bool en_odm){
	Vec3T mu_prev = mu;
	
	Mat33 Gt = {{	{1,0,0},\
			{0,1,0},\
			{0,0,1} }};
	Mat33 Vt = {{	{DELTAT,0,0},\
			{0,DELTAT,0},\
			{0,0,DELTAT} }};
	Mat33 Mt = {{	{ALPHA1*u.val[0]*u.val[0]+SIGF*SIGF,0,0},\
			{0,ALPHA2*u.val[1]*u.val[1]+SIGF*SIGF,0},\
			{0,0,ALPHA3*u.val[2]*u.val[2]+SIGF*SIGF} }};
	//estimation step
	mu = mu + u*DELTAT;
	Sigma = Gt*Sigma*Gt.trns() + Vt*Mt*Vt.trns();
	
	// === correction by beacon signals === //
	double St[n_sens];
	double z_hat[n_sens];
	if(en_cor){
		Vec3 Ht[n_sens];
		//correction step
		double Qt = SIGR*SIGR;
		for(int i = 0; i < n_sens; i++){
			z_hat[i] = (m[i]-mu).norm();
			Ht[i] = (-1*(m[i]-mu)/z_hat[i]).trns();
			St[i] = Ht[i]*Sigma*Ht[i].trns()+Qt;
			Vec3T Kt = Sigma*Ht[i].trns()/St[i];
			mu = mu + Kt*(z[i] - z_hat[i]);
			Mat33 I = {{{1,0,0},{0,1,0},{0,0,1}}};
			Sigma = (I - Kt*Ht[i])*Sigma;
		}
	}
	
	// === correction by odometry === //
	Mat33 S2t;
	Vec3T odm_hat;
	if(en_odm){
		Mat33 Q2t = {{ {SIGO*SIGO,0,0},{0,SIGO*SIGO,0},{0,0,SIGO*SIGO} }};
		Mat33 H2t = {{ {1,0,0},{0,1,0},{0,0,1} }};
		odm_hat = mu - mu_prev;
		S2t = H2t*Sigma*H2t.trns()+Q2t;
		Mat33 K2t = Sigma*H2t.trns()*S2t.inv();
		mu = mu + K2t*(odm - odm_hat);
		Mat33 I = {{ {1,0,0},{0,1,0},{0,0,1} }};
		Sigma = (I - K2t*H2t)*Sigma;
	}
	
	EKFPARAMS results;
	results.mu = mu;
	results.Sigma = Sigma;
	results.pz = 1;

	if(en_cor){
		for(int i = 0; i < n_sens; i++){
			z_hat[i] = (m[i]-mu).norm();
			results.pz *= exp(-1/2*(z[i]-z_hat[i])/St[i]*(z[i]-z_hat[i]))/sqrt(2*M_PI*St[i]);
		}
	}
	if(en_odm){
		odm_hat = mu - mu_prev;
		results.pz *= exp(-1/2*(odm-odm_hat).trns()*S2t.inv()*(odm-odm_hat))/sqrt(2*M_PI*S2t.det());
	}
	
	return results;
}

