// env_xy.cpp
//
// This file was derived from env_xyz.cpp (Ver.1.2.2-140921)
// It defines the motion/measurement model in addition to the roobt's state.
// The version of this file is integrated with that of hpp file for simplicity.
//

#include "env_xy.hpp"

// === random number seriese === //
default_random_engine generator;
//normal_distribution<double> motion_noise(0.0,SIGF);
normal_distribution<double> measur_noise(0.0,SIGR*DELTAT);
normal_distribution<double> measur_noise2(0.0,SIGR2*DELTAT);
normal_distribution<double> odom_noise(0.0,SIGO);
default_random_engine uni_generator;
uniform_real_distribution<double> uni_distribution(0.0,1.0);
void initRand(){
	unsigned seed;
	seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine g (seed);
	generator = g;
	seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine uni_g (seed);
	uni_generator = uni_g;
}

// === real world === //
static Vec2T delta_mu; // motion at each time step (it is used for odometry)

Vec2T getRealNextState(Vec2T prev, Vec2T u, bool en_float){
	normal_distribution<double> motion_noise_x(0.0,sqrt(ALPHA1)*abs(u.val[0]));
	normal_distribution<double> motion_noise_y(0.0,sqrt(ALPHA2)*abs(u.val[1]));
	Vec2T a;
	a.val[0] = motion_noise_x(generator);
	a.val[1] = motion_noise_y(generator);
	if(en_float){ // floating effect is enabled
		normal_distribution<double> floating_noise(0.0,SIGF);
		Vec2T flt;
		flt.val[0] = floating_noise(generator);
		flt.val[1] = floating_noise(generator);
		delta_mu = (u+a)*DELTAT + flt;
	}else{
		delta_mu = (u+a)*DELTAT;
	}
	return prev + delta_mu;
}

// === real measurement === //
Vec2T m[MAX_N_SENS] = {
	{{ 10000, 10000}},
	{{     0, 10000}},
	{{-10000, 10000}},
	{{ 10000,     0}},
	{{-10000,     0}},
	{{ 10000,-10000}},
	{{     0,-10000}},
	{{-10000,-10000}}};


// location of sensor
Vec2T getLocSens(int indx_sens){
	return m[indx_sens];
}
// to get measurements given the location of the robot
void getRealMeasure(Vec2T loc, Vec2T* r, int n_sens){
	for(int i = 0; i < n_sens; i++){
		Vec2T buff = m[i] - loc;
		r[i].val[0] = buff.norm() + measur_noise(generator);
		r[i].val[1] = atan2(buff.val[1],buff.val[0]) + measur_noise2(generator);
	}
}

// === real mesurement (odometry) === //
// to get odometry
Vec2T getRealOdometry(){
	Vec2T odm;
	odm.val[0] = delta_mu.val[0] + odom_noise(generator);
	odm.val[1] = delta_mu.val[1] + odom_noise(generator);
	return odm;
}


// === FastSLAM1.0 === //
PARTICLE::PARTICLE(Vec2T init_pose){
	pose = init_pose;
	for(int i = 0; i < MAX_N_SENS; i++){
		feat[i].mu = (Vec2T){{0,0}};
		feat[i].Sigma = (Mat22){{{10000000,0},{0,10000000}}};
	}
}
PARTICLE::PARTICLE(): PARTICLE((Vec2T){{0,0}}){}

Vec2T sample(Vec2T prev, Vec2T u){
	normal_distribution<double> motion_noise_x(0.0,sqrt(ALPHA1)*abs(u.val[0]));
	normal_distribution<double> motion_noise_y(0.0,sqrt(ALPHA2)*abs(u.val[1]));
	Vec2T a;
	a.val[0] = motion_noise_x(generator);
	a.val[1] = motion_noise_y(generator);
	return prev + (u+a)*DELTAT;
}

// calculate a sensor reading based on a feature location and robot's pose
Vec2T h(Vec2T mu, Vec2T pose){
	Vec2T buff = mu - pose;
	double r = buff.norm();
	double theta = atan2(buff.val[1],buff.val[0]);
	return (Vec2T){{r,theta}};
}

// calculate the feature location based on a sensor reading and the robot's pose
Vec2T h_inv(Vec2T z, Vec2T pose){
	double r = z.val[0];
	double theta = z.val[1];
	Vec2T buff = pose;
	buff.val[0] = r*cos(theta);
	buff.val[1] = r*sin(theta);
	return buff;
}
// calculate a Jacobian of sensor reading wrt robot's pose
Mat22 h_dash(Vec2T mu, Vec2T pose){
	// * h is a function of locations of feature and robot to generate a sensor reading
	// h(mu_true,pose)~=h(mu,pose)+h'(pose,mu)*(m_true-mu) = z_hat + H*(m_true - mu)
	//
	// x ---h---> mu
	//        <- m_true ->
	// given mu, m_true moves
	// so at the point of mu, function h is linearized
	//
	// if m_true grows, h grows
	// if x grows, h decreases
	// if mu grows, h decreases
	//
	double sqrt_q = (mu-pose).norm();
	double q = sqrt_q*sqrt_q;
	Mat22 H;
	H.val[0][0] = (mu.val[0] - pose.val[0])/sqrt_q;
	H.val[0][1] = (mu.val[1] - pose.val[1])/sqrt_q;
	H.val[1][0] = -1*(mu.val[1] - pose.val[1])/q;
	H.val[1][1] = (mu.val[0] - pose.val[0])/q;
	
	return H;
}

void copy_part(PARTICLE *dest, PARTICLE *src){
	dest->pose = src->pose;
	dest->omega = src->omega;
	for(int iz = 0; iz < MAX_N_SENS; iz++){
		dest->feat[iz] = src->feat[iz];
	}
}

#define IdenMat ((Mat22){{{1,0},{0,1}}})

void estimate(PARTICLE *Y, Vec2T *z, int n_sens, Vec2T u, bool init, bool enable){
	// for all particle
	for(int ip = 0; ip < N_PART; ip++){
		
		// === sampling === //
		PARTICLE *p = &Y[ip];
		p->pose = sample(p->pose, u); // sample pose
		
		//debug
		/*
		if(ip == 0){
			for(int iz = 0; iz < 1; iz++){
				cout << "iz = " << iz << ": " << z[iz].val[0]*cos(z[iz].val[1])-p->feat[iz].mu.val[0]
				     << ", " << z[iz].val[0]*sin(z[iz].val[1])-p->feat[iz].mu.val[1] << endl;
			}
		}*/
		// === update of features === //
		// initialize a jointed Q matrix
		gsl_matrix *Q_all = gsl_matrix_alloc(n_sens*2,n_sens*2);
		gsl_matrix_set_zero(Q_all); // initialize the matrix with 0
		Mat22 Q_t = (Mat22){{{SIGR*SIGR,0},{0,SIGR2*SIGR2}}};
		// for all features
		for(int iz = 0; iz < n_sens; iz++){
			if(init){// if the beginning
				p->feat[iz].mu = h_inv(z[iz], p->pose); // initialize mean
				Mat22 H = h_dash(p->feat[iz].mu, p->pose); // Jacobian
				p->feat[iz].Sigma = H.inv()*Q_t*(H.inv().trns()); // initialize covariance
			}else{
				Vec2T z_hat = h(p->feat[iz].mu, p->pose); // measurement prediction
				Mat22 H = h_dash(p->feat[iz].mu, p->pose); // calculate Jacobian
				Mat22 Q = H*(p->feat[iz].Sigma)*(H.trns())+Q_t; // measurement covariance
				if(Q.det() != 0){
					Mat22 K = p->feat[iz].Sigma*(H.trns())*(Q.inv()); // calculate Kalman gain
					p->feat[iz].mu = p->feat[iz].mu + K*(z[iz]-z_hat); // update mean
					p->feat[iz].Sigma = (IdenMat-K*H)*(p->feat[iz].Sigma); // update covariance
					
					//debug
					/*
					if(ip == 0 && iz < 1){
						cout << "z[" << iz << "]-z_hat: " << endl;
						(z[iz]-z_hat).trns().print();
						cout << "correction: " << endl;
						(K*(z[iz]-z_hat)).trns().print();
					}*/
				}
				
				// set Q
				gsl_matrix_set(Q_all,iz*2  ,iz*2  ,Q.val[0][0]);
				gsl_matrix_set(Q_all,iz*2+1,iz*2  ,Q.val[1][0]);
				gsl_matrix_set(Q_all,iz*2  ,iz*2+1,Q.val[0][1]);
				gsl_matrix_set(Q_all,iz*2+1,iz*2+1,Q.val[1][1]);
			}
		}
		
		// === importance weight === //
		if(init){
			p->omega = P_0;
		}else{
			double z_all[MAX_N_SENS*2];
			for(int iz = 0; iz < n_sens; iz++){
				Vec2T buff = p->feat[iz].mu - p->pose;
				z_all[iz*2] = buff.norm();
				z_all[iz*2+1] = atan2(buff.val[1],buff.val[2]);
			}
			double coeff = det(Q_all);
			
			if(coeff != 0){
				for(int iz = 0; iz < 2*n_sens; iz++) coeff*=2*M_PI;
				coeff = 1.0/sqrt(coeff);
				double exp_buff = 0;
				gsl_matrix *inverse = gsl_matrix_alloc(n_sens*2,n_sens*2);
				inv(inverse,Q_all);
				for(int iz = 0; iz < 2*n_sens; iz++){
					for(int jz = 0; jz < 2*n_sens; jz++){
						exp_buff += z_all[jz]*gsl_matrix_get(inverse,jz,iz)*z_all[iz];
					}
				}
				exp_buff *= -1/2;
				// update importance weight
				p->omega = coeff*exp(exp_buff);
				
				gsl_matrix_free(inverse);
			}else{
				bool match = true;
				for(int iz = 0; iz < 2*n_sens; iz++){
					if(z_all[iz] != 0) match = false;
				}
				if(match)
					p->omega = 1000000000000;
				else
					p->omega = 0;
			}
		}
		
		gsl_matrix_free(Q_all);
		
	}
	
	// === resampling ===
	if(enable){
		
		// normalize omegas
		double omega_sum = 0;
		for(int ip = 0; ip < N_PART; ip++){
			omega_sum += Y[ip].omega;
		}
		// if total is 0, assign default values
		if(omega_sum == 0){
			for(int ip = 0; ip < N_PART; ip++){
				Y[ip].omega = 1.0 / N_PART;
			}
		}else{
			for(int ip = 0; ip < N_PART; ip++){
				Y[ip].omega /= omega_sum;
			}
		}
		
		// copy particles to buffer
		PARTICLE *Y_buff = new PARTICLE[N_PART];
		for(int ip = 0; ip < N_PART; ip++){
			copy_part(&Y_buff[ip],&Y[ip]);
		}
		
		// resample
		for(int ip = 0; ip < N_PART; ip++){
			// get a random number
			double rand_key = uni_distribution(uni_generator);
			// get a corresponding one
			double buff = 0;
			int indx2pick = 0;
			while(indx2pick < N_PART-1){
				buff += Y_buff[indx2pick].omega;
				if(buff > rand_key)
					break;
				indx2pick++;
			}
			// set the chosen particle
			copy_part(&Y[ip],&Y_buff[indx2pick]);
		}
		delete Y_buff;
	}
	
}

void initPart(PARTICLE *Y, Vec2T vec){
	for(int ip = 0; ip < N_PART; ip++){
		Y[ip].pose = vec;
	}
}

Vec2T getAve(PARTICLE *Y){
	Vec2T buff = {{0,0}};
	for(int ip = 0; ip < N_PART; ip++){
		buff = buff + Y[ip].pose;
	}
	buff = buff / N_PART;
	return buff;
}
