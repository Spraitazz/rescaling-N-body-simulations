#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf.h>
//#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>
#include <fftw3.h>

//my definitions
#define pi 3.14159265358979323846264338327

#define X 0
#define Y 1
#define Z 2

#define REAL_SPACE 2221
#define REDSHIFT_SPACE 2222

#define NGP 3331
#define CIC 3332

#define NORMAL 4441
#define REVERSED 4442

#define NORMAL_BIN 5551
#define LOG_BIN 5552

#define F4 1111
#define F5 1112
#define F6 1113
#define GR 1114

//my headers
#include "structs.h"
#include "main.h"

//my code
#include "parameters.c"
#include "memory.c"
#include "functions.c"
#include "populate_haloes.c"
#include "power_spectrum.c"
#include "catalogue_input.c"
#include "fold.c"
#include "runs.c"
#include "covariance_matrices.c"
#include "CubicSpline_double.c"
#include "rescaling.c"
#include "RSD_Pk.c"
#include "halo_model.c"
#include "Dplus.c"
#include "rescaling_functions.c"
#include "ZA.c"
#include "bias.c"
//#include "HOD.c"
//#include "tests.c"

//http://arxiv.org/pdf/astro-ph/0005010v2.pdf rescaling
//http://arxiv.org/pdf/1408.1047v2.pdf rescaling
//http://arxiv.org/pdf/1308.5183v3.pdf rescaling
//http://arxiv.org/pdf/1202.5559v4.pdf HOD
//http://arxiv.org/pdf/1005.2413v2.pdf HOD parameters


int main(int argc, char **argv) {	
	//unused parameters
	argc = argc;
	argv = argv;	
	//run time
	clock_t start_time = clock();	

	//set home directory!!
	//home_directory = "/shome/jonasp/Work/Testing_GR/code/Jonas";	
	home_directory = "/home/jonas/Testing_GR/code/Jonas";
	homeDirChk();
	sprintf(temp_directory, "%s/output/tmp", home_directory);
			
	//for checking if a catalogue was input, set this to a 0
	particle_no = 0;
	halo_no = 0;
	satellite_no = 0;
	central_no = 0;
	
	//define volume for particles
	volume_limits[0] = 256.0;  // h^-1 Mpc 
	volume_limits[1] = 256.0;
	volume_limits[2] = 256.0;
	volume = volume_limits[0]*volume_limits[1]*volume_limits[2];	
	
	//define grid size, EVEN NUMBERS ONLY TO AVOID PROBLEMS, 2^n best for FFT	
	cells[0] = 128;
	cells[1] = 128;
	cells[2] = 128;	
	
	//allocate grid memory
	initGrid();	
	
	//grid size for ZA displacements
	cells_displ[0] = 128;
	cells_displ[1] = 128;
	cells_displ[2] = 128;	
		
	//spectrum bins
	spectrum_size = 30;	
		
	//start random numb generator
	init_rng();			
	
	//HOD RANDS STUFF	
	set_params();
	
	//covariance_sequence();
	//HOD_main();
	//biased_displacements();
	//ZA_RSD();
	//ST_bias();

	//rescale_testRun();
	//rescale_model();
	rescale_catalogue_to_model();
	//rescale_catalogue_to_catalogue();
	exit(0);
	
	
	
	

	
	// GET DATA FROM THE MOCKS HERE
	char in[100];
	char out[100];
	char in_final[100];	
	char out_final[100];
	char out_filename[100];		
	
	char storage[] = "/home/jonas/Testing_GR/code/Jonas/data";
	char type[] = "GR";
	bool redshift_space = true;
	double fold_factor = 32.0;
	//char tmp0[100] = "/home/jonas/Testing_GR/code/Jonas/data/GR/"
	//char tmp1[100];
	//char tmp2[100]
	
	//double foldPoints = {0.04, 0.1, 0.21, 0.6}

	if (redshift_space) {
		sprintf(out_filename, "%s_redshift_%s_%dfold.dat", type, "%d", (int)fold_factor);
	} else {
		sprintf(out_filename, "%s_%s_%dfold.dat", type, "%d",  (int)fold_factor);
	}
		
	sprintf(in, "%s/AllHalo_%s_1500Mpc_a0.7_%s.txt", mock_directory, type, "%d");
	sprintf(out, "%s/%s/%s", storage, type, out_filename);
	
	overdensity_allocateMemory();
	for (int i = 1; i < 7; i++) {
		sprintf(out_final, out, i);
		sprintf(in_final, in, i);
		printf("run %d, in file: %s, out file: %s \n", i, in_final, out_final);
		read_full_xyz(in_final, 1, mocks_format);
		if (redshift_space) toRedshift(0.7);		
		//haloes_measure_Pk(out_final, fold_factor, NGP);	
		
		
	}
	overdensity_freeMemory();	
	
	
	//fftw_cleanup();
	clock_t end_time = clock();
	double time_spent_seconds = (double)(end_time - start_time)/CLOCKS_PER_SEC;
	printf("\nDone in %.1lf seconds\n", time_spent_seconds);

	return 0;
	
	
	
	
	
	
	//char galaxy_catalogue_centrals[] = "/home/jonas/Testing_GR/code/Jonas/output/HOD/centrals_positions.dat";
	//char galaxy_catalogue_satellites[] = "/home/jonas/Testing_GR/code/Jonas/output/HOD/satellites_positions.dat";	
	///populate_halos(galaxy_catalogue_centrals, galaxy_catalogue_satellites);
			
	/*
	char out_allpower[] =  "/home/jonas/Testing_GR/code/Jonas/output/HOD/full_power.dat";
	char out_power_mmin[] = "/home/jonas/Testing_GR/code/Jonas/output/HOD/treshold_power.dat";
	char out_power_centrals[] = "/home/jonas/Testing_GR/code/Jonas/output/HOD/centrals_power.dat";
	char out_power_satellites[] = "/home/jonas/Testing_GR/code/Jonas/output/HOD/satellites_power.dat";
	
	char halos_mmin[] = "/home/jonas/Testing_GR/code/Jonas/output/HOD/halo_positions_mmin.dat";
	char sim_in[100];
	sprintf(sim_in, "%s%s", mock_directory, "AllHalo_GR_1500Mpc_a0.7_1.txt");
	pre_covariance_run(sim_in, mocks_format, out_allpower, true, false, 1);	*/
	
	//pre_covariance_run(galaxy_catalogue_satellites, galaxies_format, out_power_satellites, true, false, 0);
	//prep_variances();
	
	//select_halos_mmin(sim_in, halos_mmin);
	
	//pre_covariance_run(halos_mmin, mocks_format, out_power_mmin, true, false, 0);
	//memCheck();
	
	
	//combine_multiple(in_uf, in_f, out, 0.02, 305, 1);
	//generate_covariance(out, cov_out, 305, 1);
	//integrals_test();
	
	/*
	char in[100];
	char out_unfolded[100];
	char out_folded[100];
	char ou[100];
	char of[100];
	char cov_in[100];
	char cov_out[100];
	char storage[] = "/home/jonas/Testing_GR/code/Jonas/output/";
	char type[] = "GR";
		
	char format[] = "%s%s/%s%s_power_%s_%d.txt";
	sprintf(out_unfolded, format, storage, type, type, "", "%d", spectrum_size);
	sprintf(out_folded, format, storage, type, type, "_folded", "%d", spectrum_size);
	sprintf(in, "%sAllHalo_%s_1500Mpc_a0.7_%s.txt", mock_directory, type, "%d");
	sprintf(cov_in, "%s%s/%s_combined_%s_%d.txt", storage, type, type, "%d", spectrum_size);
	sprintf(cov_out, "%s%s/%s_covariance_%s_%d.txt", storage, type, type, "%d", spectrum_size);
	*/
	
/*
	for (int i = 1; i < 7; i++) {
		sprintf(in, "%sAllHalo_GR_1500Mpc_a0.7_%d.txt", mock_directory, i);
		sprintf(ou, out_unfolded, i);
		printf("sim %d to %s\n", i, out_unfolded);
		//better keep settings like this in a specific struct, pass that struct then
		pre_covariance_run(in, ou, true, false);	
	}
	
	for (int i = 1; i < 7; i++) {
		sprintf(in, "%sAllHalo_GR_1500Mpc_a0.7_%d.txt", mock_directory, i);
		sprintf(of, out_folded, i);
		printf("sim %d\n", i);
		//better keep settings like this in a specific struct, pass that struct then
		pre_covariance_run(in, of, true, true);	
	}
*/
	
	
	
	//combine_multiple(out_unfolded, out_folded, cov_in, 0.02, 6, 1);	

	//generate_covariance(cov_in, cov_out, 6, 1);
	
	
	
}
