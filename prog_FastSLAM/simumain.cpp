// simumain.cpp
//
// Simulation of FastSLAM1.0
//
// author: Testuya Idota
//
// Ver.0.0.1-141202: file derived from simumain.cpp (Ver.1.3.1-140921) for EKF
// Ver.0.1.0-141202: eliminate the portion of EKF
// Ver.0.9.0-141210: mostly done?
// Ver.0.9.1-141211: correlation matrix debugged
//
// compilation:
//    $ make
// to run:
//    $ ./FastSLAM1.0 [--dis_cor] [--time TIMESTEP] [--hide_table] [FILENAME]
//
// options:
//   --dis_cor: disables the resampling step
//   --time TIMESTEP: keeps simulation by the time of TIMESTEP. If no input data is available, zeros are taken.
//   --hide_table: hides the table of results. The program displays basic information besides the correlation matrix.
//
// notes:
//   given no filename, the simulator loads data from input.txt.
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include "MV2.hpp"
#include "env_xy.hpp"


using namespace std;

int main(int argc, char** argv){

	// setup for environment
	Vec2T realstate = {{0,0}};
	int n_sens = MAX_N_SENS;
	char fname[64] = "input.txt";

	// setup for estimation
	bool f_enable = true; // enable the resampling step
	int t_limit = 20;
	bool f_hide_table = false; // hide the table of resutls
	
	// processing arguments
	for(int iarg = 1; iarg < argc; iarg++){
		if(argv[iarg][0] == '-' && argv[iarg][1] == '-'){ // if argument is a parameter
			if(strcmp(argv[iarg] + 2, "dis_cor") == 0){ // to disable the correction step
				f_enable = false;
			}else if(strcmp(argv[iarg] + 2, "time") == 0){ // to set a limit of time
				iarg++;
				t_limit = atoi(argv[iarg]);
			}else if(strcmp(argv[iarg] + 2, "hide_table") == 0){ // to hide the table of results
				f_hide_table = true;
			}else{
				cout << "INVALID PARAMETERS: " << argv[iarg] << endl;
				return -1;
			}
		}else{ // supposed to be an file name
			strcpy(fname,argv[iarg]);
		}
	}

	ifstream ifs;
	ifs.open(fname);
	if(!ifs.good()){
		cout << "FILE OPEN ERROR:" << fname << endl;
		return -1;
	}
	
	// initialize the random number seriese
	initRand();
	
	// === header of the results === //
	PARTICLE *Y = new PARTICLE[N_PART];
	initPart(Y,(Vec2T){{0,0}});
	Vec2T est = getAve(Y);
	
	if(f_hide_table){
		cout << "processing";
	}else{
		if(f_enable)
			cout << "resampling enabled" << endl;
		else
			cout << "resampling disabled" << endl;
		cout << "# of sensors: " << n_sens << endl;
		cout << "time ||=input(mm/sec)==||=true path(mm)==||=estimated(mm)==||===error(mm)====||  error" << endl;
		cout << "(sec)||   dx      dy   ||   x       y    ||    x       y   ||    x       y   ||" << endl;
		cout << fixed;
		cout << setw(5) << setprecision(1) << 0 << "||";
		cout << "                ||";
		cout << setw(7) << setprecision(0) << realstate.val[0] << "  " \
		     << setw(7) << setprecision(0) << realstate.val[1] << "||";
		cout << setw(7) << setprecision(0) << est.val[0] << "  " \
		     << setw(7) << setprecision(0) << est.val[1] << "||";
		cout << setw(7) << setprecision(0) << abs(realstate.val[0] - est.val[0]) << "  " \
		     << setw(7) << setprecision(0) << abs(realstate.val[1] - est.val[1]) << "||";
		cout << setw(7) << setprecision(0) << (realstate - est).norm() << "  ";
		cout << endl;

		cout << "-----||----------------" << \
		             "||----------------" << \
		             "||----------------" << \
		             "||----------------" << \
		             "||----------------" << endl;
	}
	
	// === main iteration of FastSLAM1.0 === //
	int count = 1;
	while(1){
		if(count > t_limit) break;
		
		// get input data
		Vec2T u;
		ifs >> u.val[0] >> u.val[1];
		
		if(ifs.eof()){
			if(count <= t_limit)
				u = (Vec2T){{0,0}};
			else
				break;
		}
		
		// update for real world
		realstate = getRealNextState(realstate, u, false);
		
		// get sensor readings
		Vec2T r[n_sens];
		getRealMeasure(realstate,r,n_sens);
		//Vec2T odom = getRealOdometry();
		
		// FastSLAM1.0
		estimate(Y, r, n_sens, u, (count == 1), f_enable);
		est = getAve(Y);
		
		// output
		if(f_hide_table){
			if(count % 10 == 0)
				cout << "|" << flush;
			else
				cout << "." << flush;
		}else{
			cout << fixed;
			cout << setw(5) << setprecision(1) << count*DELTAT << "||";
			cout << setw(7) << setprecision(0) << u.val[0] << "  " \
			     << setw(7) << setprecision(0) << u.val[1] << "||";
			cout << setw(7) << setprecision(0) << realstate.val[0] << "  " \
			     << setw(7) << setprecision(0) << realstate.val[1] << "||";
			cout << setw(7) << setprecision(0) << est.val[0] << "  " \
			     << setw(7) << setprecision(0) << est.val[1] << "||";
			cout << setw(7) << setprecision(0) << abs(realstate.val[0] - est.val[0]) << "  " \
			     << setw(7) << setprecision(0) << abs(realstate.val[1] - est.val[1]) << "||";
			cout << setw(7) << setprecision(0) << (realstate - est).norm() << "  ";
			cout << endl;

			
			if(count % 10 == 0)
				cout << "-----||----------------" << \
				             "||----------------" << \
				             "||----------------" << \
				             "||----------------" << \
				             "||----------------" << endl;
		}
		
		count++;
	}
	
	ifs.close();
	
	
	// === calculation of  correlation matrix ===
	// calculate means
	Vec2T ave_pose = {{0,0}};
	Vec2T ave_feat[n_sens];
	for(int iz = 0; iz < n_sens; iz++) ave_feat[iz] = {{0,0}};
	for(int ip = 0; ip < N_PART; ip++){
		ave_pose = ave_pose + Y[ip].pose;
		for(int iz = 0; iz < n_sens; iz++) ave_feat[iz] = ave_feat[iz] + Y[ip].feat[iz].mu;
	}
	ave_pose = ave_pose / N_PART;
	for(int iz = 0; iz < n_sens; iz++) ave_feat[iz] = ave_feat[iz] / N_PART;
	
	// create a covariance matrix of joint space over the robot's pose and features
	gsl_matrix* Sigma_all = gsl_matrix_alloc(2*(n_sens+1),2*(n_sens+1));
	gsl_matrix_set_zero(Sigma_all); //initialize with zeros
	
	double value1, value2;
	
	for(int ip = 0; ip < N_PART; ip++){
		for(int i = 0; i < 2*(n_sens+1); i++){
			if(i < 2)
				value1 = Y[ip].pose.val[i%2] - ave_pose.val[i%2];
			else{
				value1 = Y[ip].feat[i/2-1].mu.val[i%2] - ave_feat[i/2-1].val[i%2];
			}
			for(int j = 0; j < 2*(n_sens+1); j++){
				if(j < 2)
					value2 = Y[ip].pose.val[j%2] - ave_pose.val[j%2];
				else{
					value2 = Y[ip].feat[j/2-1].mu.val[j%2] - ave_feat[j/2-1].val[j%2];
				}
				
				gsl_matrix_set(Sigma_all,i,j,gsl_matrix_get(Sigma_all,i,j)+value1*value2/N_PART);
			}
		}
	}
	
	// === operation with the elements ===
	gsl_matrix *diag = gsl_matrix_alloc(2*(n_sens+1),2*(n_sens+1));
	gsl_matrix_set_zero(diag);
	for(int i = 0; i < 2*(n_sens+1); i++){
		gsl_matrix_set(diag,i,i,1.0/sqrt(gsl_matrix_get(Sigma_all,i,i)));
	}
	
	// buffer of result
	gsl_matrix *sub_results = gsl_matrix_alloc(2*(n_sens+1),2*(n_sens+1));
	gsl_matrix *results = gsl_matrix_alloc(2*(n_sens+1),2*(n_sens+1));
	
	// calculate correlation matrix
	MatMat(sub_results,diag,Sigma_all);
	MatMat(results,sub_results,diag);
	
	// display the results
	cout << endl << "Correlation Matrix" << endl;
	printMat(results);

	// release the memory
	gsl_matrix_free(Sigma_all);
	gsl_matrix_free(diag);
	gsl_matrix_free(sub_results);
	gsl_matrix_free(results);
	
	delete Y;
	
	cout << "simulation end" << endl;
	
	return 0;
}
