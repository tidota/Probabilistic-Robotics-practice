// simumain.cpp
//
// Main code of EKF simulation.
// Arguments are the input filename and parameters. if no file name is given, the program loads "input.txt"
// The input file format is a table of velocities for x, y, and z (mm/sec) at each time step.
//     e.g.,
//          	0	0	0
//          	0	0	0
//          	0	0	0
//          	100	100	100
//          	100	100	100
//          	100	100	100
//          	100	100	100
//          	...
//
// parameters are as follows:
//	--dis_cor: disables the correction step in EKF
//	--n_sens: sets the number of sensors to be used (needs to be equal to or less than MAX_N_SENS)
//	           this parameter must be followed by an exact number.
//	--init_stat: set initial state, followed by three parameters(xyz). (e.g., "--init_stat 1.4 1.0 0.0")
//	--init_mu: set intial estimation, followed by three parameters(xyz). (e.g., "--init_mu 1.0 1.2 2.5")
//	--init_sig: set initial covariances, followed by 9 parameters(a11, a12, a13, a21, ..., a31, a32, a33) 
//	             e.g., "--init_sig 30.0 0 0 0 25.0 0 0 0 35"
//	--kidnap: forcefully set a state at a specific time step, followed by time step, x, y, and z 
//	           e.g., "--kidnap 10 100 200 300" to set x=100, y=200, z=300 at t=10 * Delta t
//	--float: enables the floating effect
//
// As output, the program displays a table of input, true path, estimated path by EKF, error, and measurement likelihood.
//
// For details of settings of the simulator, refer to env_xyz.hpp/cpp
//
// to compile the codes:
//      make
// to run the program:
//	./EKFsim <parameters>
//        or
//      ./EKFsim <input file name> <parameters>
//
//	e.g., 
// 	   ./EKFsim input.txt --n_sens 4
//	   ./EKFsim input.txt --dis_cor
//
//
// author: Testuya Idota
//
// Ver.0.0.1-140909: file created
// Ver.0.9.0-140910: mostly done including output
// Ver.0.9.1-140910: argument added
// Ver.0.9.2-140911: minor updates of output format + some documents
// Ver.1.0.0-140911: arguments added
// Ver.1.1.0-140913: arguments to set states added
// Ver.1.1.1-140913: bug around arugments fixed
// Ver.1.2.0-140916: motion model modified, optional floating effect added
// Ver.1.2.1-140916: erro of each dimension added to output
// Ver.1.3.0-140919: odometry added
// Ver.1.3.1-140921: minor updates of measurement belief

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include "MV3.hpp"
#include "env_xyz.hpp"


using namespace std;

int main(int argc, char** argv){

	// setup for the initial state
	//   the initial estimation is supposed to be identical to the initial true values
	Vec3T realstate = {{0,0,0}};
	Vec3T mu = realstate;
	Mat33 Sigma = {{{50,0,0},{0,50,0},{0,0,50}}};
	bool en_cor = true;
	int n_sens = MAX_N_SENS;
	char fname[64] = "input.txt";
	bool f_kidn = false; //if true, kidnapped state is inserted at a specified time step
	int t_kidn = 0; // time step to insert kidnapped state
	Vec3T stat_kidn = {{0,0,0}}; // kidnapped state
	bool f_flt = false; // if true, the motion model is affected by the floating effects.
	bool f_odm = false; // if true, the odometry correction is enabled.
	
	// processing arguments
	for(int iarg = 1; iarg < argc; iarg++){
		if(argv[iarg][0] == '-' && argv[iarg][1] == '-'){ // if argument is a parameter
			if(strcmp(argv[iarg] + 2, "dis_cor") == 0){ // to disable the correction step
				en_cor = false;
			}else if(strcmp(argv[iarg] + 2, "n_sens") == 0){ // to set the # of sensors to be used
				n_sens = atoi(argv[iarg+1]);
				if(n_sens > MAX_N_SENS){
					cout << "the number of sensors needs not to exceed " << MAX_N_SENS << "." << endl;
					cout << MAX_N_SENS << " was set instead." << endl;
					n_sens = MAX_N_SENS;
				}else if(n_sens < 1){
					cout << "the number of sensors needs to be more than 0." << endl;
					cout << MAX_N_SENS << " was set instead." << endl;
					n_sens = MAX_N_SENS;
				}
				iarg++;
			}else if(strcmp(argv[iarg] + 2, "init_stat") == 0){ // to set initial state
				realstate.val[0] = atof(argv[iarg+1]);
				realstate.val[1] = atof(argv[iarg+2]);
				realstate.val[2] = atof(argv[iarg+3]);
				iarg += 3;
			}else if(strcmp(argv[iarg] + 2, "init_mu") == 0){ // to set initial estimation
				mu.val[0] = atof(argv[iarg+1]);
				mu.val[1] = atof(argv[iarg+2]);
				mu.val[2] = atof(argv[iarg+3]);
				iarg += 3;
			}else if(strcmp(argv[iarg] + 2, "init_sig") == 0){ // to set initial covariances
				for(int iparam = 0; iparam < 9; iparam++)
					Sigma.val[(int)(iparam/3)][iparam%3] = atof(argv[iarg+1+iparam]);
				iarg+=9;
			}else if(strcmp(argv[iarg] + 2, "kidnap") == 0){ // to set kidnapped state
				f_kidn = true;
				t_kidn = atoi(argv[iarg+1]);
				stat_kidn.val[0] = atof(argv[iarg+2]);
				stat_kidn.val[1] = atof(argv[iarg+3]);
				stat_kidn.val[2] = atof(argv[iarg+4]);
				iarg+=4;
			}else if(strcmp(argv[iarg] + 2, "float") == 0){ // to enable the floating effect
				f_flt = true;
			}else if(strcmp(argv[iarg] + 2, "en_odm") == 0){ // to enable odometry correction
				f_odm = true;
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
	if(en_cor)
		cout << "# of sensors: " << n_sens << endl;
	else
		cout << "correction step in EKF algorithm is disabled." << endl;
	cout << "time ||===== input(mm/sec) =====||==== true path(mm) ======||===== estimated(mm) =====||======= error(mm) =======||  error ,   pz" << endl;
	cout << "(sec)||   dx       dy      dz   ||   x        y        z   ||    x       y        z   ||    x       y        z   ||" << endl;
	cout << fixed;
	cout << setw(5) << setprecision(1) << 0 << "||";
	cout << "                         ||";
	cout << setw(7) << setprecision(0) << realstate.val[0] << "  " \
	     << setw(7) << setprecision(0) << realstate.val[1] << "  " \
	     << setw(7) << setprecision(0) << realstate.val[2] << "||";
	cout << setw(7) << setprecision(0) << mu.val[0] << "  " \
	     << setw(7) << setprecision(0) << mu.val[1] << "  " \
	     << setw(7) << setprecision(0) << mu.val[2] << "||";
	cout << setw(7) << setprecision(0) << abs(realstate.val[0] - mu.val[0]) << "  " \
	     << setw(7) << setprecision(0) << abs(realstate.val[1] - mu.val[1]) << "  " \
	     << setw(7) << setprecision(0) << abs(realstate.val[2] - mu.val[2]) << "||";
	cout << setw(7) << setprecision(0) << (realstate - mu).norm() << "  ";
	cout << endl;

	cout << "-----||-------------------------" << \
	             "||-------------------------" << \
	             "||-------------------------" << \
	             "||-------------------------" << \
	             "||--------------------" << endl;

	// === main iteration of EKF === //
	int count = 1;
	while(1){
		if(f_kidn && count == t_kidn)
			realstate = stat_kidn;
		
		// get input data
		Vec3T u;
		ifs >> u.val[0] >> u.val[1] >> u.val[2];

		if(ifs.eof())
			break;
		
		// update for real world
		realstate = getRealNextState(realstate, u, f_flt);
		double r[n_sens];
		getRealMeasure(realstate,r,n_sens);
		Vec3T odom = getRealOdometry();
		
		// EKF
		EKFPARAMS results = EKFlocalization(mu,Sigma,u,r,n_sens,en_cor,odom,f_odm);
		mu = results.mu;
		Sigma = results.Sigma;
		
		// output
		cout << fixed;
		cout << setw(5) << setprecision(1) << count*DELTAT << "||";
		cout << setw(7) << setprecision(0) << u.val[0] << "  " \
		     << setw(7) << setprecision(0) << u.val[1] << "  " \
		     << setw(7) << setprecision(0) << u.val[2] << "||";
		cout << setw(7) << setprecision(0) << realstate.val[0] << "  " \
		     << setw(7) << setprecision(0) << realstate.val[1] << "  " \
		     << setw(7) << setprecision(0) << realstate.val[2] << "||";
		cout << setw(7) << setprecision(0) << mu.val[0] << "  " \
		     << setw(7) << setprecision(0) << mu.val[1] << "  " \
		     << setw(7) << setprecision(0) << mu.val[2] << "||";
		cout << setw(7) << setprecision(0) << abs(realstate.val[0] - mu.val[0]) << "  " \
		     << setw(7) << setprecision(0) << abs(realstate.val[1] - mu.val[1]) << "  " \
		     << setw(7) << setprecision(0) << abs(realstate.val[2] - mu.val[2]) << "||";
		cout << setw(7) << setprecision(0) << (realstate - mu).norm() << "  ";
		if(en_cor || f_odm)
			cout << setw(9) << setprecision(6) << 100.0*results.pz << "%";
		else
			cout << "    NA    ";
		cout << endl;

		if(count % 10 == 0)
			cout << "-----||-------------------------" << \
			             "||-------------------------" << \
			             "||-------------------------" << \
			             "||-------------------------" << \
			             "||--------------------" << endl;

		count++;
	}
	
	ifs.close();
	
	return 0;
}
