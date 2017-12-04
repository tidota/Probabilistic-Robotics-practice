// simumain.cpp
//
// Main code of Grid Localization simulation.
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

//************************** to edit *****************************
// parameters are as follows:
//	--dis_cor: disables the correction step
//	--en_odm: enables the odometry
//	--n_sens: sets the number of sensors to be used (needs to be equal to or less than MAX_N_SENS)
//	           this parameter must be followed by an exact number.
//	--init_stat: set initial state, followed by three parameters(xyz). (e.g., "--init_stat 1.4 1.0 0.0")
//	--kidnap: forcefully set a state at a specific time step, followed by time step, x, y, and z 
//	           e.g., "--kidnap 10 100 200 300" to set x=100, y=200, z=300 at t=10 * Delta t
//	--init_GL: set initial estimation, followed by x, y, and z (the nearest grid is assined with 1 and others with 0)
//	--float: enables the floating effect (141001: not guaranteed to work in this version)
//
// As output, the program displays a table of input, true path, estimated path by EKF, error, and measurement likelihood.
//
// For details of settings of the simulator, refer to env_xyz.hpp/cpp
//
// to compile the codes:
//      make
// to run the program:
//	./GLsim <parameters>
//        or
//      ./GLsim <input file name> <parameters>
//
//	e.g., 
// 	   ./GLsim input.txt --n_sens 4
//	   ./GLsim input.txt --dis_cor
//
// ***********************************************************

//
// author: Testuya Idota
//
// Ver.0.0.1-140927: derived from EKF simulator
// Ver.0.5.0-140930: GL implemented and updated format of output
// Ver.0.6.0-141001: position tracking problem implemented
// Ver.0.7.0-141002: floating effects added
// Ver.0.7.1-141003: inflation of noises as variables
// Ver.0.7.2-141003: inflation of noises as constants
// Ver.0.8.0-141005: odometry implemented

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
	
	bool en_cor = true;
	int n_sens = MAX_N_SENS;
	char fname[64] = "input.txt";
	bool f_kidn = false; //if true, kidnapped state is inserted at a specified time step
	int t_kidn = 0; // time step to insert kidnapped state
	Vec3T stat_kidn = {{0,0,0}}; // kidnapped state
	bool f_flt = false; // if true, the motion model is affected by the floating effects.
	bool f_odm = false;// if true, the odometry correction is enabled.
	
	
	Vec3T initest; // initial estimation
	bool f_initest = false;
	
	
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
			}else if(strcmp(argv[iarg] + 2, "kidnap") == 0){ // to set kidnapped state
				f_kidn = true;
				t_kidn = atoi(argv[iarg+1]);
				stat_kidn.val[0] = atof(argv[iarg+2]);
				stat_kidn.val[1] = atof(argv[iarg+3]);
				stat_kidn.val[2] = atof(argv[iarg+4]);
				iarg+=4;
			}else if(strcmp(argv[iarg] + 2, "init_GL") == 0){ // to set the initial estimation
				initest.val[0] = atof(argv[iarg+1]);
				initest.val[1] = atof(argv[iarg+2]);
				initest.val[2] = atof(argv[iarg+3]);
				f_initest = true;
				iarg+=3;
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
	
	// initalize random number series
	initRand();
	
	ifstream ifs;
	ifs.open(fname);
	if(!ifs.good()){
		cout << "FILE OPEN ERROR:" << fname << endl;
		return -1;
	}

	// allocation for grids
	GL_MAIN *gl_main = new GL_MAIN();
	
	if(f_initest) // initial estimation is specified
		gl_main->setpos(initest);

	// === header of the results === //
	if(en_cor)
		cout << "# of sensors: " << n_sens << endl;
	else
		cout << "correction step in EKF algorithm is disabled." << endl;

	cout << " interval of grids: " << Unit_GRID << endl;
	cout << " range of field " << endl;
	cout << " x(mm): " << ORG_GRID_X << " to " << ORG_GRID_X + (N_GRIDS_X-1)*Unit_GRID << endl;
	cout << " y(mm): " << ORG_GRID_Y << " to " << ORG_GRID_Y + (N_GRIDS_Y-1)*Unit_GRID << endl;
	cout << " z(mm): " << ORG_GRID_Z << " to " << ORG_GRID_Z + (N_GRIDS_Z-1)*Unit_GRID << endl;
	
	cout << "-----||--------- input ---------" << \
	             "||------- true path -------" << \
	             "||------- estimated -------" << \
	             "||--------- error ---------" << \
	             "||---- det of cov ----" << endl;
	cout << "     ||   x        y       z    " << \
	             "||   x        y       z    " << \
	             "||   x        y       z    " << \
	             "||   x        y       z    " << \
	             "||" << endl;
	
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
		
		// Grid Localization
		gl_main->localize(u, r, n_sens, en_cor, odom, f_odm, f_flt);
		
		// output
		RESULTPARAMS params = gl_main->getResults();
		
		cout << fixed;
		cout << setw(5) << setprecision(1) << count*DELTAT << "||";
		cout << setw(7) << setprecision(0) << u.val[0] << "  " \
		     << setw(7) << setprecision(0) << u.val[1] << "  " \
		     << setw(7) << setprecision(0) << u.val[2] << "||";
		cout << setw(7) << setprecision(0) << realstate.val[0] << "  " \
		     << setw(7) << setprecision(0) << realstate.val[1] << "  " \
		     << setw(7) << setprecision(0) << realstate.val[2] << "||";
		cout << setw(7) << setprecision(0) << params.mu.val[0] << "  " \
		     << setw(7) << setprecision(0) << params.mu.val[1] << "  " \
		     << setw(7) << setprecision(0) << params.mu.val[2] << "||";
		cout << setw(7) << setprecision(0) << abs(realstate.val[0]-params.mu.val[0]) << "  " \
		     << setw(7) << setprecision(0) << abs(realstate.val[1]-params.mu.val[1]) << "  " \
		     << setw(7) << setprecision(0) << abs(realstate.val[2]-params.mu.val[2]) << "||";
		cout << setw(7) << setprecision(0) << params.Sigma.det() << endl;
		
		if(count % 10 == 0)
			cout << "-----||-------------------------" << \
			             "||-------------------------" << \
			             "||-------------------------" << \
			             "||-------------------------" << \
			             "||--------------------" << endl;
		
		count++;
	}
	ifs.close();
	
	// memory reliese of grids
	delete gl_main;
	
	return 0;
}
