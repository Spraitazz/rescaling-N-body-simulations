

/*
int halo_mass_func() {
	double dm, this_mass;
	int bins, index;
	bins = 50;
	char temp_out[100];
	//sprintf(temp_out, "%s%s", tempdir, "hmf_temp.dat");
	const char format[] = "%le \t %le \n";
	double* n_m = calloc((size_t)bins, sizeof(*n_m));
	double* hmf_m = calloc((size_t)(bins-1), sizeof(*hmf_m));
	
	dm = (Mmax - Mmin)/(double)bins;
	for (int i = 0; i < particle_no; i++) {
		this_mass = particles[i][6];
		index = (int) floor((this_mass - Mmin) / dm);
		if (index == bins) index = bins - 1;		
		n_m[index] += 1.0;
	}
	
	//n = N/V
	for (int i = 0; i < bins; i++) {
		n_m[i] /= volume;
	}
	
	FILE *out = fopen(temp_out, "w");
	for (int i = 0; i < bins-1; i++) {
		hmf_m[i] = ((n_m[i+1] - n_m[i])/dm);
		//take midpoint of bin
		this_mass = Mmin + ((double)i + 0.5)*dm;
		fprintf(out, format, this_mass, hmf_m[i]); 
	}
	
	fclose(out);
	//HMF_spline = input_spline_file(temp_out, format, false);	
	//remove(temp_out);
	return 0;
}

*/

//GR -> F5
//scale box -> ZA -> HOD?
int rescale_testRun() {	

	double s_test = 0.9;
	double z_test = 0.3;
	double z_tmp;
	
	char Pk_current_path[100];
	sprintf(Pk_current_path, "%s/data/linear_matter_pk_sig8_0.593_z_0.75.dat", home_directory);
	
	char Pk_target_path[100];
	sprintf(Pk_target_path, "%s/data/linear_matter_pk_0.76.dat", home_directory);	
		
	//Pk spline for target, model
	Pk_current = input_spline_file(Pk_current_path, Pk_model_format, NORMAL);
	//Pk_target = input_spline_file(Pk_target_path, Pk_model_format, false);	
	Pk_target = input_spline_file(Pk_current_path, Pk_model_format, NORMAL);		
	
	z_current = 0.75;
	z_target = 0.75;	
	
	//find optimal s, z		
	R1_primed = 0.05;
	R2_primed = 5.0;
	vary_z_current = true;
	//global variables
	k_bins = 300;
	z_bins = 200;
	Dplus_bins = 300;
	R_bins = 300;
	k_min = 0.01;
	k_max = 1.0;
	z_min = 0.0;
	z_max = 3.0;
	current_gravity = GR;
	target_gravity = GR;
	
	//binning info
	rescaling_k_bin_info = prep_bins(k_min, k_max, k_bins, LOG_BIN);	
	rescaling_z_bin_info = prep_bins(z_min, z_max, z_bins, NORMAL_BIN);
	rescaling_R_bin_info = prep_bins(0.6*R1_primed, 1.5*R2_primed, R_bins, LOG_BIN); 
	
	//integral workspaces
	z_to_t_workspace = gsl_integration_workspace_alloc(z_to_t_workspace_size);
	OLV_workspace = gsl_integration_workspace_alloc(OLV_workspace_size);
	
	//redshift-time relation reversed:	
	prep_redshift_spline(rescaling_z_bin_info);
	
	//for target simulation:
	Pk_target = prep_Pk_constz(target_gravity, z_target, z_test, Pk_target, rescaling_k_bin_info);
	Pk_target_extended = extend_spline_model(Pk_target);	
	variance_spline_target = prep_variances_test(rescaling_R_bin_info, Pk_target_extended, s_test);
	
	//for current simulation:		
	Pk_splines_zBins = malloc((size_t)z_bins * sizeof(*Pk_splines_zBins));
	Pk_splines_zBins_extended = malloc((size_t)z_bins * sizeof(*Pk_splines_zBins_extended));
	variance_splines_zBins = malloc((size_t)z_bins * sizeof(*variance_splines_zBins));
	
	double min_sq_diff = DBL_MAX;
	double sq_diff;

	for (int i = 0; i < z_bins; i++) {
		z_tmp = bin_to_x(rescaling_z_bin_info, i);
		printf("Preparing OLV at different redshifts. zbin %d, z: %lf \n", i, z_tmp);
		
		Pk_splines_zBins[i] = prep_Pk_constz(current_gravity, z_current, z_tmp, Pk_current, rescaling_k_bin_info);
		
		sq_diff = pow(splint_generic(Pk_splines_zBins[i], 0.1) - splint_generic(Pk_target, 0.1), 2.0);
		if (sq_diff < min_sq_diff) {
			min_sq_diff = sq_diff;
			z_guess = z_tmp;
		}
		
		Pk_splines_zBins_extended[i] = extend_spline_model(Pk_splines_zBins[i]);

		variance_splines_zBins[i] = prep_variances(rescaling_R_bin_info, Pk_splines_zBins_extended[i]);
	}
	
	printf("z guess: %lf \n", z_guess);
	
		
	
	dsq_multimin(vary_z_current, z_guess, variance_spline_target, variance_splines_zBins);

	free(Pk_splines_zBins);
	free(Pk_splines_zBins_extended);
	free(variance_splines_zBins);
	gsl_integration_workspace_free(z_to_t_workspace);	
	gsl_integration_workspace_free(OLV_workspace);	
	return 0;		
}

int rescale_model() {	
	
	char Pk_current_path[100];
	sprintf(Pk_current_path, "%s/data/linear_matter_pk_sig8_0.593_z_0.75.dat", home_directory);
	
	char Pk_target_path[100];
	sprintf(Pk_target_path, "%s/data/linear_matter_pk_0.76.dat", home_directory);	
	
	//Pk spline for target, model
	Pk_current = input_spline_file(Pk_current_path, Pk_model_format, NORMAL);
	Pk_target = input_spline_file(Pk_target_path, Pk_model_format, NORMAL);	
	
	z_current = 0.75;
	z_target = 0.76;	
	
	//find optimal s, z		
	R1_primed = m_to_R(1e12, z_target);
	R2_primed = m_to_R(5e15, z_target);
	
	printf("R1': %le, R2': %le \n", R1_primed, R2_primed);
	
	//varying z in current sim, box scaling always for current sim
	vary_z_current = true;

	k_bins = 300;
	z_bins = 200;
	Dplus_bins = 300;
	R_bins = 300;
	k_min = 0.01;
	k_max = 1.0;
	z_min = 0.0;
	z_max = 5.0;
	current_gravity = GR;
	target_gravity = GR;	
	
	//binning info
	rescaling_k_bin_info = prep_bins(k_min, k_max, k_bins, LOG_BIN);	
	rescaling_z_bin_info = prep_bins(z_min, z_max, z_bins, NORMAL_BIN);
	rescaling_R_bin_info = prep_bins(0.6*R1_primed, 1.5*R2_primed, R_bins, LOG_BIN); 
	
	//integral workspaces
	z_to_t_workspace = gsl_integration_workspace_alloc(z_to_t_workspace_size);
	OLV_workspace = gsl_integration_workspace_alloc(OLV_workspace_size);
	
	//redshift-time relation reversed:	
	prep_redshift_spline(rescaling_z_bin_info);
	
	//for target simulation:	
	Pk_target_extended = extend_spline_model(Pk_target);	
	variance_spline_target = prep_variances(rescaling_R_bin_info, Pk_target_extended);
	
	//for current simulation:		
	Pk_splines_zBins = malloc((size_t)z_bins * sizeof(*Pk_splines_zBins));
	Pk_splines_zBins_extended = malloc((size_t)z_bins * sizeof(*Pk_splines_zBins_extended));
	variance_splines_zBins = malloc((size_t)z_bins * sizeof(*variance_splines_zBins));
	
	
	
	double min_sq_diff = DBL_MAX;
	double sq_diff;

	double z_cur;
	for (int i = 0; i < z_bins; i++) {
		z_cur = bin_to_x(rescaling_z_bin_info, i);
		printf("zbin %d, z: %lf \n", i, z_cur);

		Pk_splines_zBins[i] = prep_Pk_constz(current_gravity, z_current, z_cur, Pk_current, rescaling_k_bin_info);	
		
		sq_diff = pow(splint_generic(Pk_splines_zBins[i], 0.1) - splint_generic(Pk_target, 0.1), 2.0);
		if (sq_diff < min_sq_diff) {
			min_sq_diff = sq_diff;
			z_guess = z_cur;
		}

		//Pk_splines_zBins_extended[i] = extend_spline_model(Pk_splines_zBins[i]);
		
		//variance_splines_zBins[i] = prep_variances(rescaling_R_bin_info, Pk_splines_zBins_extended[i]);
	}
	
	printf("z guess: %lf \n", z_guess);
	//R1_primed = m_to_R(1e12, z_target);
	//R2_primed = m_to_R(5e15, z_target);
		
	
	for (int i = 0; i < z_bins; i++) {
		z_cur = bin_to_x(rescaling_z_bin_info, i);
		printf("zbin %d, z: %lf \n", i, z_cur);		

		Pk_splines_zBins_extended[i] = extend_spline_model(Pk_splines_zBins[i]);
		
		variance_splines_zBins[i] = prep_variances(rescaling_R_bin_info, Pk_splines_zBins_extended[i]);
	}
	

	//minimization to find s, z_rescaled (global vars)
	dsq_multimin(vary_z_current, z_guess, variance_spline_target, variance_splines_zBins);
	/*
	double dsqmin = dsq_minimized_value;
	for (int i = 0; i < 10; i++) {
		printf("\n RETRYING MINIMIZATION \n");
		dsq_multimin(vary_z_current, variance_spline_target, variance_splines_zBins);
		if (dsq_minimized_value < dsqmin) {
			dsqmin = dsq_minimized_value;
		}
	}
	*/

	int z_ind_final = x_to_bin(rescaling_z_bin_info, z_rescaled);
	
	Pk_current_extended = extend_spline_model(Pk_current);
	
	//current 
	char* current_outPath = prep_path("/data/rescaling/delsq_current.dat");
	print_delsq(Pk_current_extended, current_outPath, 0.01, 0.1, 200, false);
	free(current_outPath);
	
	//target
	char* target_outPath = prep_path("/data/rescaling/delsq_target.dat");
	print_delsq(Pk_target_extended, target_outPath, 0.01, 0.1, 200, false);
	free(target_outPath);
	
	//current with redshift scale
	char* current_outPath_z = prep_path("/data/rescaling/delsq_current_z.dat");
	print_delsq(Pk_splines_zBins_extended[z_ind_final], current_outPath_z, 0.01, 0.1, 200, false);
	free(current_outPath_z);	
	
	//current after redshift AND box scale
	char current_outPath_zs[100];
	sprintf(current_outPath_zs, "%s/data/rescaling/delsq_current_zs.dat", home_directory);
	print_delsq(Pk_splines_zBins_extended[z_ind_final], current_outPath_zs, 0.01, 0.1, 200, true); 

	
	generate_sigma_plots();	
	
	free(Pk_splines_zBins);
	free(Pk_splines_zBins_extended);
	free(variance_splines_zBins);
	free(rescaling_k_bin_info);
	free(rescaling_R_bin_info);
	free(rescaling_z_bin_info);
	gsl_integration_workspace_free(z_to_t_workspace);	
	gsl_integration_workspace_free(OLV_workspace);	
	return 0;		
}




//CATALOGUE TO MODEL

int rescale_catalogue_to_model() {		

	char Pk_target_path[100];
	sprintf(Pk_target_path, "%s/data/linear_matter_pk_sig8_0.593_z_0.75.dat", home_directory);		

	char sim_current[100];	
	sprintf(sim_current, "%s/%s", mock_directory, "AllHalo_GR_1500Mpc_a0.7_1.txt");	
	
	char outFile_scaled[100];
	sprintf(outFile_scaled, "%s/data/rescaling/measured_pk_scaled.dat", home_directory);
	char outFile[100];
	sprintf(outFile, "%s/data/rescaling/measured_pk.dat", home_directory);
	
	z_current = a_to_z(0.7);
	z_target = 0.75;
	
	Particle_Catalogue *current_catalogue;	
	
	//current simulation data input
	current_catalogue = input_catalogue_file(sim_current, 1, mocks_format);	
	
	//min and max masses, here for current catalogue, but actually must be for target catalogue
	mass_min_max(current_catalogue);
	printf("mmin: %le,  mmax: %le\n", Mmin,  Mmax);
	R1_primed = m_to_R(Mmin, z_current);
	R2_primed = m_to_R(Mmax, z_current);
	printf("R1': %le, R2': %le \n", R1_primed, R2_primed);		
	
	overdensity_allocateMemory();
	haloes_measure_Pk(current_catalogue, outFile, 1.0, CIC);			
	overdensity_freeMemory();
	
	Pk_current = input_spline_file(outFile, k_Pkmono_format, NORMAL);
	/*
	for (int i = 0; i < 100; i++) {
		double kcur = Pk_current->xmin + (double)i * ((Pk_current->xmax - Pk_current->xmin) / 99.0);
		printf("k: %lf, Pk: %lf \n", kcur, splint_generic(Pk_current, kcur));
	}
	*/
	//remove(outFile);	
	Pk_target = input_spline_file(Pk_target_path, Pk_model_format, NORMAL);	
	
	printf("current kmin: %lf, kmax: %lf \n", Pk_current->xmin, Pk_current->xmax);
	
	//find optimal s, z		
	vary_z_current = true;
	
	k_bins = 300;
	z_bins = 200;
	Dplus_bins = 300;
	R_bins = 300;
	k_min = 0.02;
	k_max = 0.2;
	z_min = 0.0;
	z_max = 4.0;	
	current_gravity = GR;
	target_gravity = GR;
	
	/*
	if (current_gravity != GR || target_gravity != GR) {
		prep_Dplus_kz();	
	}
	*/
	
	rescaling_k_bin_info = prep_bins(k_min, k_max, k_bins, LOG_BIN);
	rescaling_z_bin_info = prep_bins(z_min, z_max, z_bins, NORMAL_BIN);
	rescaling_R_bin_info = prep_bins(0.6*R1_primed, 1.5*R2_primed, R_bins, LOG_BIN); 
	//HOD_mass_bin_info = prep_bins(Mmin, Mmax, 100, NORMAL_BIN);
	
	//integral workspaces
	z_to_t_workspace = gsl_integration_workspace_alloc(z_to_t_workspace_size);
	OLV_workspace = gsl_integration_workspace_alloc(OLV_workspace_size);
	
	//redshift-time relation reversed:	
	prep_redshift_spline(rescaling_z_bin_info);
	
	//for current simulation:		
	Pk_splines_zBins = malloc((size_t)z_bins * sizeof(*Pk_splines_zBins));
	Pk_splines_zBins_extended = malloc((size_t)z_bins * sizeof(*Pk_splines_zBins_extended));
	variance_splines_zBins = malloc((size_t)z_bins * sizeof(*variance_splines_zBins));
	Pk_current_extended = extend_spline(Pk_current);	
	
	//for target simulation:	
	Pk_target_extended = extend_spline_model(Pk_target);	
	variance_spline_target = prep_variances(rescaling_R_bin_info, Pk_target_extended);

	//finding a redshift to start minimising from
	double z_cur, z_guess, sq_diff, min_sq_diff;
	min_sq_diff = DBL_MAX;
	for (int i = 0; i < z_bins; i++) {	
		z_cur = bin_to_x(rescaling_z_bin_info, i);
		Pk_splines_zBins[i] = prep_Pk_constz(current_gravity, z_current, z_cur, Pk_current, rescaling_k_bin_info);	
		
		//check more points
		sq_diff = pow(splint_generic(Pk_splines_zBins[i], 0.1) - splint_generic(Pk_target, 0.1), 2.0);
		if (sq_diff < min_sq_diff) {
			min_sq_diff = sq_diff;
			z_guess = z_cur;
		}			
	}	
	
	printf("z guess: %lf \n", z_guess);			
	
	//prepare OLV bins at different redshifts
	for (int i = 0; i < z_bins; i++) {
		z_cur = bin_to_x(rescaling_z_bin_info, i);
		printf("zbin %d, z: %lf \n", i, z_cur);		

		Pk_splines_zBins_extended[i] = extend_spline(Pk_splines_zBins[i]);		
		variance_splines_zBins[i] = prep_variances(rescaling_R_bin_info, Pk_splines_zBins_extended[i]);
	}
	
	//current 
	char current_outPath[100];
	sprintf(current_outPath, "%s/data/rescaling/cat_to_model/delsq_current.dat", home_directory);
	print_delsq(Pk_current_extended, current_outPath, k_min, k_max, 200, false);
	
	//target
	char target_outPath[100];
	sprintf(target_outPath, "%s/data/rescaling/cat_to_model/delsq_target.dat", home_directory);
	print_delsq(Pk_target_extended, target_outPath, k_min, k_max, 200, false);	
	
	//minimization
	dsq_multimin(vary_z_current, z_guess, variance_spline_target, variance_splines_zBins);
	int z_ind_final = x_to_bin(rescaling_z_bin_info, z_rescaled);
	
	//current s scaled
	char current_outPath_s[100];
	sprintf(current_outPath_s, "%s/data/rescaling/cat_to_model/delsq_current_s.dat", home_directory);
	print_delsq(Pk_current_extended, current_outPath_s, k_min, k_max, 200, true);
	
	//scale box
	volume_limits[0] *= s;
	volume_limits[1] *= s;
	volume_limits[2] *= s;
	volume = volume_limits[0]*volume_limits[1]*volume_limits[2];
	
	//ensure PBC
	for (int i = 0; i < current_catalogue->particle_no; i++) {
		PBC(&(current_catalogue->particles[i].x), &(current_catalogue->particles[i].y), &(current_catalogue->particles[i].z));
	}	
	
	//scale masses, OMEGA M DEPENDS ON Z!!!!!!!!!!!
	for (int i = 0; i < particle_no; i++) {		
		current_catalogue->particles[i].mass *= pow(s, 3.0) * (omega_m_primed/omega_m_z(z_rescaled));
	}
	
	mass_min_max(current_catalogue);
	printf("AFTER RESCALING. M min: %le, M max: %le\n", Mmin, Mmax);
	
	overdensity_allocateMemory();
	haloes_measure_Pk(current_catalogue, outFile_scaled, 1.0, CIC);			
	overdensity_freeMemory();	
	
	SplineInfo* Pk_current_rescaled = input_spline_file(outFile_scaled, k_Pkmono_format, NORMAL);
	
	SplineInfo_Extended* Pk_current_rescaled_extended = extend_spline(Pk_current_rescaled);
	for (int i = 0; i < 200; i++) {
		double kcur = 0.02 + (double)i * ((0.2 - 0.02) / 199.0);
		printf("k: %lf, Pk: %lf \n", kcur, volume*splint_Pk(Pk_current_rescaled_extended, kcur));
	}
	//remove(outFile);
	
	//current s scaled, measured
	char current_outPath_s_measured[100];
	sprintf(current_outPath_s_measured, "%s/data/rescaling/cat_to_model/delsq_current_s_measured.dat", home_directory);
	print_delsq(Pk_current_rescaled_extended, current_outPath_s_measured, k_min, k_max, 200, false);
	
	/*
	//current with redshift scale
	char current_outPath_z[100];
	sprintf(current_outPath_z, "%s/data/rescaling/cat_to_model/delsq_current_z.dat", home_directory);
	print_delsq(Pk_splines_zBins_extended[z_ind_final], current_outPath_z, k_min, k_max, 200, false);
	
	//deviation from wanted pk
	double someNumber = 0.0;
	double k_cur;
	for (int i = 0; i < k_bins; i++) {
		k_cur = bin_to_x(rescaling_k_bin_info, i);
		someNumber += pow(splint_Pk(Pk_splines_zBins_extended[z_ind_final], k_cur) - splint_Pk(Pk_target_extended, k_cur), 2.0);
	}

	printf("some number after z scale: %le \n", someNumber);
	someNumber = 0.0;
	
	//current after redshift AND box scale
	char current_outPath_zs[100];
	sprintf(current_outPath_zs, "%s/data/rescaling/cat_to_model/delsq_current_zs.dat", home_directory);
	print_delsq(Pk_splines_zBins_extended[z_ind_final], current_outPath_zs, k_min, k_max, 200, true); 
	
	for (int i = 0; i < k_bins; i++) {
		k_cur = bin_to_x(rescaling_k_bin_info, i);
		someNumber += pow(splint_Pk(Pk_splines_zBins_extended[z_ind_final], k_cur) - splint_Pk(Pk_target_extended, k_cur/s), 2.0);
	}
	
	printf("some number after zs scale: %le \n", someNumber);
	*/
	
	
	exit(0);	
	
	
	
	
	
	
	
	/*	
	//scale velocities
	for (int i = 0; i < particle_no; i++) {
		particles[i][3] = particles[i][3]*s*(H_primed*a_primed*fg_primed)/(H*a*fg);
		particles[i][4] = particles[i][4]*s*(H_primed*a_primed*fg_primed)/(H*a*fg);
		particles[i][5] = particles[i][5]*s*(H_primed*a_primed*fg_primed)/(H*a*fg);	
	}
	*/
	

	//calc_b_eff(z_rescaled);
	//bias_plot();
	
	
	
	
	/*
	populateGridCIC();
	gridIntoOverdensity();
	overdensity_pk(false);
	pk_to_file_logbin(out_final);
*/
	
	/*
	for (int i = 0; i < particle_no; i++) {
		PBC(&particles[i][0], &particles[i][1], &particles[i][2]);
	}	
	*/
	
	//scale masses, velocities with s
	//does this come before or after HOD??
	//NEED SOME PARAMETERS HERE!!
	//scale_masses();
	//scale_velocities();
	
	
	//Pk_out(out_final);
	
	/*
	
	//move around haloes
	//displacements and velocities! Don't necessarily need to free this yet		
	fftw_allocateMemory_ZA();
	restore_disp_ZA();
	fftw_free(displacements);
	fftw_free(displacements_fourier);
	*/
	
	free(Pk_splines_zBins);
	free(Pk_splines_zBins_extended);
	free(variance_splines_zBins);
	free(rescaling_k_bin_info);
	free(rescaling_R_bin_info);
	free(rescaling_z_bin_info);
	gsl_integration_workspace_free(z_to_t_workspace);	
	gsl_integration_workspace_free(OLV_workspace);	
	return 0;		
}



//CATALOGUE TO CATALOGUE

int rescale_catalogue_to_catalogue() {		

	char Pk_target_path[100];
	sprintf(Pk_target_path, "%s/data/linear_matter_pk_sig8_0.593_z_0.75.dat", home_directory);		

	char sim_current[100];	
	sprintf(sim_current, "%s%s", mock_directory, "AllHalo_GR_1500Mpc_a0.7_1.txt");	
	
	//char out_final[] = "/home/jonas/Testing_GR/code/Jonas/output/rescaling/GR_rescaled_test.dat";
	
	char tmpFile[100];
	sprintf(tmpFile, "%s/data/rescaling/tmp_unfolded.dat", home_directory);
	char tmpFile2[100];
	sprintf(tmpFile2, "%s/data/rescaling/tmp_folded.dat", home_directory);
	char outFile[100];
	sprintf(outFile, "%s/data/rescaling/measured_pk.dat", home_directory);
	
	z_current = a_to_z(0.7);
	z_target = 0.75;	
	
	//current simulation data input
	//read_full_xyz(sim_current, 1, mocks_format);
	//min and max masses, MUST COME BEFORE CUMSUMS
	//mass_min_max(&particles, particle_no);
	printf("mmin: %le,  mmax: %le\n", Mmin,  Mmax);
	R1_primed = m_to_R(Mmin, z_current);
	R2_primed = m_to_R(Mmax, z_current);
	printf("R1': %le, R2': %le \n", R1_primed, R2_primed);		
	
	overdensity_allocateMemory();
	//haloes_measure_Pk(&particles, particle_no, tmpFile, 1.0, CIC);	
	//haloes_measure_Pk(tmpFile2, 4.0, false);
	//combine_folded(tmpFile, tmpFile2, pk_format_folding, outFile, 0.0, 1.1, 0.3, true);	//MAKE SURE THESE ARE MANUALLY SET FOR BEST P(k)	
	overdensity_freeMemory();
	remove(tmpFile);
	remove(tmpFile2);
	
	Pk_current = input_spline_file(outFile, "%le \t %le \t %*e \n", false);
	Pk_target = input_spline_file(Pk_target_path, Pk_model_format, false);
	remove(outFile);	
	
	//find optimal s, z		
	vary_z_current = true;
	
	k_bins = 200;
	z_bins = 300;
	Dplus_bins = 300;
	R_bins = 300;
	k_min = 0.001;
	k_max = 10.0;
	z_min = 0.0;
	z_max = 5.0;
	
	current_gravity = GR;
	target_gravity = GR;
	
	rescaling_k_bin_info = prep_bins(k_min, k_max, k_bins, LOG_BIN);
	rescaling_z_bin_info = prep_bins(z_min, z_max, z_bins, NORMAL_BIN);
	rescaling_R_bin_info = prep_bins(0.6*R1_primed, 1.5*R2_primed, R_bins, LOG_BIN); 
	HOD_mass_bin_info = prep_bins(Mmin, Mmax, 100, NORMAL_BIN);
	
	//for current simulation:		
	Pk_splines_zBins = malloc((size_t)z_bins * sizeof(*Pk_splines_zBins));
	Pk_splines_zBins_extended = malloc((size_t)z_bins * sizeof(*Pk_splines_zBins_extended));
	variance_splines_zBins = malloc((size_t)z_bins * sizeof(*variance_splines_zBins));	

	double z_tmp;
	for (int i = 0; i < z_bins; i++) {
		z_tmp = bin_to_x(rescaling_z_bin_info, i);
		printf("zbin %d, z: %lf \n", i, z_tmp);
		
		Pk_splines_zBins[i] = Pk_current;
		prep_Pk_constz(current_gravity, z_current, z_tmp, Pk_splines_zBins[i], rescaling_k_bin_info);	
		Pk_splines_zBins_extended[i] = extend_spline_model(Pk_splines_zBins[i]);	
			 
		variance_splines_zBins[i] = prep_variances(rescaling_R_bin_info, Pk_splines_zBins_extended[i]);
	}
	
	//for target simulation:	
	Pk_target_extended = extend_spline_model(Pk_target);
	
	//variance_spline_target = prep_variances(&rescaling_R_bin_info, &Pk_target_extended);	

	//variance_spline_const = variance_spline_target;		
	
	//dsq_multimin(vary_z_current, &variance_spline_const, &variance_splines_zBins);
	
	//got s, z_rescaled
	int z_ind_final = x_to_bin(rescaling_z_bin_info, z_rescaled);
	variance_spline_current = variance_splines_zBins[z_ind_final];
	variance_spline_current_reversed = reverse_spline(variance_spline_current); //IMPORTANT!!
	
	/*
	//scale box
	volume_limits[0] *= s;
	volume_limits[1] *= s;
	volume_limits[2] *= s;
	volume = volume_limits[0]*volume_limits[1]*volume_limits[2];
	
	for (int i = 0; i < particle_no; i++) {
		PBC(&particles[i][0], &particles[i][1], &particles[i][2]);
	}
		
	//scale masses
	for (int i = 0; i < particle_no; i++) {
		//OMEGA M DEPENDS ON Z!!!!!!!!!!!
		particles[i][6] = pow(s, 3.0)*omega_m_primed*particles[i][6]/omega_m_z(z_rescaled);
	}
	
	//scale velocities
	for (int i = 0; i < particle_no; i++) {
		particles[i][3] = particles[i][3]*s*(H_primed*a_primed*fg_primed)/(H*a*fg);
		particles[i][4] = particles[i][4]*s*(H_primed*a_primed*fg_primed)/(H*a*fg);
		particles[i][5] = particles[i][5]*s*(H_primed*a_primed*fg_primed)/(H*a*fg);	
	}
	*/
	

	//calc_b_eff(z_rescaled);
	//bias_plot();
	
	
	
	
	/*
	populateGridCIC();
	gridIntoOverdensity();
	overdensity_pk(false);
	pk_to_file_logbin(out_final);
*/
	
	/*
	for (int i = 0; i < particle_no; i++) {
		PBC(&particles[i][0], &particles[i][1], &particles[i][2]);
	}	
	*/
	
	//scale masses, velocities with s
	//does this come before or after HOD??
	//NEED SOME PARAMETERS HERE!!
	//scale_masses();
	//scale_velocities();
	
	
	//Pk_out(out_final);
	
	/*
	
	//move around haloes
	//displacements and velocities! Don't necessarily need to free this yet		
	fftw_allocateMemory_ZA();
	restore_disp_ZA();
	fftw_free(displacements);
	fftw_free(displacements_fourier);
	*/
	
	remove(tmpFile);
	return 0;		
}



/*
int HOD_catalogue() {
	double cos_theta_sat, theta_sat, phi_sat, r_sat, unif_rnd, mean_r_sat;
	double erf_central, mean_N_satellite, mean_N_central, halo_mass;
	double x_sat, y_sat, z_sat, dm_nfw;
	//velocity assignment for satellites variables
	double result_I, error, f_cdm, sigma_sq_sat, sigma_sq_sat_tcdm, vx, vy, vz, sigma;
	int has_central, satellite_number, total_sats, total_centrals, nfw_bins, nfw_index;
	
	//prepare NFW profile cumulative sums
	nfw_bins = 100;
	nfw_cumsums_reversed = malloc((size_t)nfw_bins * sizeof(*nfw_cumsums_reversed));
	nfw_cumsums_catalogue(100, "/home/jonas/Testing_GR/code/Jonas/output/nfw_cumsums/");
	dm_nfw = (Mmax - Mmin)/(double)nfw_bins;
	
	//find out how many centrals and allocate them
	total_centrals = 0;
	for (int i = 0; i < particle_no; i++) {
		halo_mass = particles[i][6];
		erf_central = gsl_sf_erf((log10(halo_mass) - log_Mmin)/sigma_logm);
		mean_N_central = 0.5 * (1.0 + erf_central);		
		has_central = (int)gsl_ran_bernoulli(gsl_ran_r, (unsigned int)mean_N_central);
		if (has_central == 1) {
			particles[total_centrals][0] = particles[i][0];
			particles[total_centrals][1] = particles[i][1];
			particles[total_centrals][2] = particles[i][2];
			particles[total_centrals][3] = particles[i][3];
			particles[total_centrals][4] = particles[i][4];
			particles[total_centrals][5] = particles[i][5];
			particles[total_centrals][6] = particles[i][6];			
			total_centrals += 1;
		}
	}
	
	//either this, or manually zero the rest of memory
	central_no = total_centrals;
	particles = realloc(particles, (size_t)central_no * sizeof(*particles));
	
	//find out how many satellites to allocate for each HALO OR CENTRAL?????
	total_sats = 0;
	int* no_satellites = calloc((size_t)central_no, sizeof(*no_satellites));	
	for (int i = 0; i < central_no; i++) { //<----- OR LOOP THROUGH HALOES
		halo_mass = particles[i][6];
		mean_N_satellite = pow((halo_mass - M0)/M1, alpha);
		satellite_number = (int)gsl_ran_poisson(gsl_ran_r, (unsigned int)mean_N_satellite);			
		no_satellites[i] = satellite_number;
		total_sats += satellite_number;		
	}
	
	//realloc to add memory for satellites
	satellite_no = total_sats;	
	particles = realloc(particles, (size_t)(central_no + satellite_no) * sizeof(*particles));
	//calloc zeros out the values of satellites
	for (int i = central_no; i < central_no + satellite_no; i++) {
		particles[i] = (double*) calloc(7, sizeof(**particles));
	}		
	
	//prepare integral for velocity assignment of satellites		
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);	
	gsl_function F;
  	F.function = &I_integral_HOD;
  	F.params = 0;   
  	//end velocity assignment stuff			

	//add satellites, this counter is to place satellites in the now expanded particles array
	int satellite_index = central_no;
	for (int i = 0; i < central_no; i++) {
		halo_mass = particles[i][6];
		nfw_index = (int) floor((halo_mass - Mmin)/dm_nfw);
		//each central i will have no_satellites[i] as found before
		for (int j = 0; j < no_satellites[i]; j++) {
			//cos theta is uniform random in [-1,1]
			cos_theta_sat = randomDouble(-1.0, 1.0);
			theta_sat = acos(cos_theta_sat);
			
			//phi is uniform random in [0,2pi]
			phi_sat = randomDouble(0.0, 2.0*M_PI);
			
			//selects a cumulative sum in range [cumsum_min, cumsum_max]
			//splint (inverse) then relates that cumsum to some r value between 0 and r_virial
			unif_rnd = randomDouble(nfw_cumsums_reversed[nfw_index].xmin, nfw_cumsums_reversed[nfw_index].xmax);

			r_sat = splint_generic(nfw_cumsums_reversed[nfw_index], unif_rnd);
			
			//track the mean radial displacement
			mean_r_sat += r_sat;			
		
			//actual position is with respect to the i-th central
			x_sat = r_sat * sin(theta_sat) * cos(phi_sat) + particles[i][0];
			y_sat = r_sat * sin(theta_sat) * sin(phi_sat) + particles[i][1];
			z_sat = r_sat * cos(theta_sat) + particles[i][2];			

			PBC(&x_sat, &y_sat, &z_sat);
		
			//place the satellite in particle array		
			particles[satellite_index][0] = x_sat;
			particles[satellite_index][1] = y_sat;
			particles[satellite_index][2] = z_sat;	
			
			//finally, assign velocities			
			gsl_integration_qags(&F, r_sat/rs_global, 1000.0, 0, 1e-5, 1000, w, &result_I, &error); 
	
			f_cdm = log(1.0 + cdm_global) - cdm_global/(1.0 + cdm_global);	
			//sigma squared from the paper with the rest of HOD recipe
			sigma_sq_sat = (G*halo_mass*cdm_global*cdm_global*r_sat/(R_virial_global*R_virial_global*f_cdm))*pow(1.0 + cdm_global*r_sat/R_virial_global, 2.0)*result_I;
			
			//sigma squared from tcdm paper, DIFFERENT
			sigma_sq_sat_tcdm = (G*halo_mass/(3.0*R_virial_global))*cdm_global*(1.0 - 1.0/(1.0+pow(cdm_global,2.0)) - 2.0*log(1.0 + cdm_global)/(1.0 + cdm_global))/(2.0*pow(log(1.0 + cdm_global) - cdm_global/(1.0 + cdm_global),2.0));
			
			double ratio = sigma_sq_sat_tcdm/sigma_sq_sat;
			if (ratio > 1.1 || ratio < 0.9) {
				//printf("ratio: %lf\n", ratio);
			}
	
			//this is the relative velocity to the halo, so add together with the host halo velocity
			sigma = sqrt(sigma_sq_sat);
			vx = gsl_ran_gaussian(gsl_ran_r, sigma_sq_sat);
			vy = gsl_ran_gaussian(gsl_ran_r, sigma_sq_sat);
			vz = gsl_ran_gaussian(gsl_ran_r, sigma_sq_sat);	
			
			particles[satellite_index][3] = particles[i][3] + vx;
			particles[satellite_index][4] = particles[i][4] + vy;
			particles[satellite_index][5] = particles[i][5] + vz;		
			
			satellite_index += 1;
		}		
	}
	
	return 0;
}
*/
