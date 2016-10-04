#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <fftw3.h>

//my definitions
#define pi 3.14159265358979323846264338327
#define F4 0
#define F5 1
#define F6 2
#define GR 3

//my headers
#include "main.h"

//my code
#include "parameters.c"
#include "memory.c"
#include "functions.c"
#include "power_spectrum.c"
#include "catalogue_input.c"
#include "fold.c"
#include "runs.c"
#include "covariance_matrices.c"
#include "CubicSpline_double.c"
#include "rescaling.c"
#include "halo_model.c"
#include "Dplus.c"
#include "rescaling_functions.c"
#include "ZA.c"
//#include "halomodel.c"
//#include "mockGalaxyCats.c"
//#include "jonas.c"

//#include "minimisation.c"
//#include "halomodel_Mike.c"
//#include "tests.c"

int main(int argc, char **argv) {	
	//unused parameters
	argc = argc;
	argv = argv;	
	//run time
	clock_t start_time = clock();	

	//set parameters for current cosmology
	//set_params();
	
	//set home directory!!
	//home_directory = "/shome/jonasp/Work/Testing_GR/code/Jonas";
	home_directory = "/home/jonas/Testing_GR/code/Jonas";
	//homeDirChk();
			
	//for checking if a catalogue was input, set this to a 0
	particle_no = 0;
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
		
	//spectrum bins, particles
	spectrum_size = 50;
	particle_no = 1000000;	
	initParticles();
		
	//start random numb generator
	init_rng();			
	
	
	
	//HOD RANDS STUFF	
	
	
	ZA_test();
	//ZA_test_Mikes();
	exit(0);
	
	
	
	
	
	//Variable to hold spline information for splinting
	SplineInfo Pk_spline_temp;// = input_spline_file(Pk_in_model, Pk_model_format, false);	
	SplineInfo Pk_spline_temp2 = Pk_spline_temp;
	SplineInfo Pk_spline_temp3 = Pk_spline_temp; 
	SplineInfo Pk_spline_temp4 = Pk_spline_temp;
	SplineInfo Pk_spline_temp5 = Pk_spline_temp;	
	//also adds regression by powerlaw for small and high k
	double z_test = 5.0;
	k_bins = 300;
	z_bins = 200;
	Dplus_bins = 1000;
	int R_bins = 500;
	k_min = 0.001;
	k_max = 11.35;
	z_min = 0.0;
	z_max = 10.0;
	rescaling_k_bin_info = prep_bins(k_min, k_max, k_bins, true);	//logbin
	rescaling_z_bin_info = prep_bins(z_min, z_max, z_bins, false);
	
	prep_Pk_constz(GR, 0.0, &Pk_spline_temp, &rescaling_k_bin_info);
	prep_Pk_constz(GR, 2.0, &Pk_spline_temp2, &rescaling_k_bin_info);
	prep_Pk_constz(GR, 4.0, &Pk_spline_temp3, &rescaling_k_bin_info);
	prep_Pk_constz(GR, 6.0, &Pk_spline_temp4, &rescaling_k_bin_info);
	prep_Pk_constz(GR, 8.0, &Pk_spline_temp5, &rescaling_k_bin_info);
	
	SplineInfo_Extended f1;
	SplineInfo_Extended f2;
	SplineInfo_Extended f3;
	SplineInfo_Extended f4;
	SplineInfo_Extended f5;	

	prep_Pk_model_spline(&Pk_spline_temp, &f1);
	prep_Pk_model_spline(&Pk_spline_temp2, &f2);
	prep_Pk_model_spline(&Pk_spline_temp3, &f3);
	prep_Pk_model_spline(&Pk_spline_temp4, &f4);
	prep_Pk_model_spline(&Pk_spline_temp5, &f5);
	
	FILE *ftmp = fopen("/home/jonas/Testing_GR/code/Jonas/output/tmp/test_Pk_model_spline.dat", "w");	
	
	
	double k_cur, Pk_cur, Pk_cur2, Pk_cur3, Pk_cur4, Pk_cur5;
	double logkmin = log(k_min);
	double logkmax = log(k_max);
	double logdk = (logkmax - logkmin) / (double)(k_bins - 1);
	double logk;
	for (int i = 0; i < k_bins; i++) {
		logk = logkmin + (double)i * logdk;
		k_cur = bin_to_x(&rescaling_k_bin_info, i);
		//k_cur = exp(logk);
		Pk_cur = volume*splint_Pk(&f1, k_cur);
		Pk_cur2 = volume*splint_Pk(&f2, k_cur);
		Pk_cur3 = volume*splint_Pk(&f3, k_cur);
		Pk_cur4 = volume*splint_Pk(&f4, k_cur);
		Pk_cur5 = volume*splint_Pk(&f5, k_cur);
		fprintf(ftmp, "%le \t %le \t %le \t %le \t %le \t %le \n", k_cur, Pk_cur, Pk_cur2, Pk_cur3, Pk_cur4, Pk_cur5);
		
	
	}
	
	fclose(ftmp);
	
	exit(0);
	
	printf("splined Pk\n");	
	
	rescaling_R_bin_info = prep_bins(0.001, 30.0, 500, true); //logbin
	prep_variance_splines(&Pk_spline_current, &variance_spline_current, NULL, NULL, false, z_test);
	reverse_spline(&variance_spline_current, &variance_spline_current_reversed);

	

	
	//for 2-halo term with ZA, this sets R_nl to the value at sigma_nonlin = 1.0
	//set_ZA_params();
	R_nl = splint_generic(&variance_spline_current_reversed, 1.0);	
	printf("Rnl is: %le\n", R_nl);	
	
	sigma_exp_smoothed = sqrt(expected_variance_smoothed(R_nl, &Pk_spline_current));
	printf("expected variance in smoothed displacement: %lf \n", sigma_exp_smoothed);

	//double halo_mass = 5E+13;
	
	char out_model[100];
	sprintf(out_model, "%s%s", home_directory, "HOD_rands/model_pk.dat");
	haloModel_out(out_model, 0.025, 0.6, 100);		
	
	//stuff for one-halo term NFW etc.
	//set_halo_params(halo_mass);
	//nfw cumsums for this halo mass, can be done outside of here
	//char nfw_cumsums_out[100];
	//sprintf(nfw_cumsums_out, "/home/jonas/Testing_GR/code/Jonas/output/nfw_cumsums/cumsums_mass_%1.1e.dat", halo_mass);
	//prep_nfw_cumsums(nfw_cumsums_out);
	//nfw_spline_reverse = input_spline(nfw_cumsums_out, "%le \t %le \n", true); //reversed
	//remove(nfw_cumsums_out);
	char out_this[100];
	for (int i = 0; i < 50; i++) {
		printf("\n run %d \n", i);
		char out_temp[] = "/home/jonas/Testing_GR/code/Jonas/output/HOD_rands/";
		randoms_halo_model(out_temp, i);
		/*
		sprintf(out_this, "%s%d%s", out_temp, i, ".dat");
		
		populateParticles();	
	
		fftw_allocateMemory_ZA();	
		restore_disp_ZA();	
		fftw_free(displacements);
		fftw_free(displacements_fourier);

		//allocate FFT memory
		fftw_allocateMemory();

		populateGridCIC();
		gridIntoOverdensity();
		overdensity_pk(false);	
		pk_to_file_logbin(out_this);	
		
		clearPk();		
		clearGrid();
		clearParticles();
	*/
	}
	
	
	
	/*
	char root_path_rands[100];
	sprintf(root_path_rands, "%s%s", home_directory, "HOD_rands/"); 
	for (int i = 0; i < 5; i++) {
		printf("run %d\n", i);		
		randoms_halo_model(root_path_rands, i);
	}	
	*/
	
	
	/* GET DATA FROM THE MOCKS HERE
	char in[100];
	char out[100];
	char in_final[100];	
	char out_final[100];
	char out_filename[100];		
	
	char storage[] = "/home/jonas/Testing_GR/code/Jonas/output/plots";
	char type[] = "F6";
	bool redshift_space = true;

	if (redshift_space) {
		sprintf(out_filename, "%s_combined_redshift_%s.dat", type, "%d");
	} else {
		sprintf(out_filename, "%s_combined_%s.dat", type, "%d");
	}
		
	sprintf(in, "%sAllHalo_%s_1500Mpc_a0.7_%s.txt", mock_directory, type, "%d");
	sprintf(out, "%s/%s/%s", storage, type, out_filename);
	
	fftw_allocateMemory();
	for (int i = 1; i < 7; i++) {
		printf("run %d\n", i);
		sprintf(out_final, out, i);
		sprintf(in_final, in, i);
		read_full_xyz(in_final, 1, mocks_format);
		if (redshift_space) toRedshift();		
		halos_run(out_temp_pk, 1.0);		
		halos_run(out_temp_folded_pk, 4.0);
		combine_folded(out_temp_pk, out_temp_folded_pk, out_final, 0.0, 1.2, 0.45, true);	
		remove(out_temp_pk);
		remove(out_temp_folded_pk);
		clearMemory();	
		fftw_clearMemory();	
	}
	*/
	
	
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
