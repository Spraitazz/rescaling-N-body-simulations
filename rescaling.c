

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
/*
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
*/

/*
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
	//R1_primed = m_to_R(1e12, z_target);
	//R2_primed = m_to_R(5e15, z_target);
	
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
	
	double dsqmin = dsq_minimized_value;
	for (int i = 0; i < 10; i++) {
		printf("\n RETRYING MINIMIZATION \n");
		dsq_multimin(vary_z_current, variance_spline_target, variance_splines_zBins);
		if (dsq_minimized_value < dsqmin) {
			dsqmin = dsq_minimized_value;
		}
	}
	

	int z_ind_final = x_to_bin(rescaling_z_bin_info, z_rescaled);
	
	Pk_current_extended = extend_spline_model(Pk_current);
	
	//current 
	char* current_outPath = prep_path("/data/rescaling/delsq_current.dat");
	//print_delsq(Pk_current_extended, current_outPath, 0.01, 0.1, 200, false);
	free(current_outPath);
	
	//target
	char* target_outPath = prep_path("/data/rescaling/delsq_target.dat");
	//print_delsq(Pk_target_extended, target_outPath, 0.01, 0.1, 200, false);
	free(target_outPath);
	
	//current with redshift scale
	char* current_outPath_z = prep_path("/data/rescaling/delsq_current_z.dat");
	//print_delsq(Pk_splines_zBins_extended[z_ind_final], current_outPath_z, 0.01, 0.1, 200, false);
	free(current_outPath_z);	
	
	//current after redshift AND box scale
	char current_outPath_zs[100];
	sprintf(current_outPath_zs, "%s/data/rescaling/delsq_current_zs.dat", home_directory);
	//print_delsq(Pk_splines_zBins_extended[z_ind_final], current_outPath_zs, 0.01, 0.1, 200, true); 

	
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
*/



int rescale_randoms_model_to_model() {	

	char Pk_current_path[100];
	sprintf(Pk_current_path, "%s/data/models/linear_matter_pk_sig8_0.593_z_0.75.dat", home_directory);
	Spline *Pk_current = input_spline_file(Pk_current_path, Pk_model_format, NORMAL, MODEL);	
	//this will go from redshift 0.75 to redshift 3.5 later	

	Parameters *cur_params = malloc(sizeof(*cur_params));		
	cur_params->H0 = 100.0*0.678;
	cur_params->omega_m_0 = 0.24;
	cur_params->omega_v_0 = 0.76;
	cur_params->omega_r_0 = 0.0;
	cur_params->omega = 1.0;	
	cur_params->gamma = 0.545;
	cur_params->z = 3.5;	
	
	char Pk_target_path[100];
	sprintf(Pk_target_path, "%s/data/models/matter_pk_om_m_015_z_10_sigma8_069.dat", home_directory);
	Spline *Pk_target = input_spline_file(Pk_target_path, Pk_model_format, NORMAL, MODEL);	
	
	
	Parameters *targ_params = malloc(sizeof(*targ_params));
	targ_params->H0 = 100.0*0.678;
	targ_params->omega_m_0 = 0.15;
	targ_params->omega_v_0 = 0.85;
	targ_params->omega_r_0 = 0.0;
	targ_params->omega = 1.0;	
	targ_params->gamma = 0.545;
	targ_params->z = 1.0;
	
	Mmin = 2e12;
	Mmax = 2e15;
	printf("mmin: %le,  mmax: %le\n", Mmin,  Mmax);
	R1_primed = m_to_R(Mmin, targ_params);
	R2_primed = m_to_R(Mmax, targ_params);
	printf("R1': %le, R2': %le \n", R1_primed, R2_primed);	
	
	//find optimal s, z		
	bool vary_z_current = true;	
	int k_bins = 300;
	int z_bins = 50;
	Dplus_bins = 300;
	int R_bins = 300;
	double k_min = 0.06; //kfund
	double k_max = 0.4; //knyq/4
	double z_min = 0.0;
	double z_max = 5.0;	
	int current_gravity = GR;
	int target_gravity = GR;	
	
	BinInfo *k_binInfo = prep_bins(k_min, k_max, k_bins, LOG_BIN);
	BinInfo *z_binInfo = prep_bins(z_min, z_max, z_bins, NORMAL_BIN);
	BinInfo *R_binInfo = prep_bins(0.6*R1_primed, 1.5*R2_primed, R_bins, LOG_BIN);
	
	Pk_current = prep_Pk_constz(current_gravity, 0.75, cur_params->z, Pk_current, k_binInfo, z_binInfo, cur_params);
	
	Spline *variance_current = prep_variances(R_binInfo, Pk_current);
	Spline *variance_target = prep_variances(R_binInfo, Pk_target);
	Spline *variance_spline_current_reversed = reverse_spline(variance_current);	
	R_nl = splint_generic(variance_spline_current_reversed, 1.0);		
	printf("Rnl is: %le \n", R_nl);


	//particles, ZA disp field storage
	particle_no = 128*128*128;	
	double** disps_init;
	
	char Pk_current_model_out[100];
	sprintf(Pk_current_model_out, "%s/data/rescaling/model_to_model/Pk_current_model.dat", home_directory);
	//BinInfo* Pk1_bins = prep_bins(0.03, 0.5, 200, LOG_BIN);
	haloModel_out(Pk_current_model_out, k_binInfo, 0.0, 1.0, REDSHIFT_SPACE, Pk_current, cur_params);
	
	Particle_Catalogue *current_catalogue;
	
	char Pk_current_out[100];	
	sprintf(Pk_current_out, "%s/data/rescaling/model_to_model/Pk_current.dat", home_directory);

	current_catalogue = populateParticles_crystalLattice(particle_no);

	//ZA + vels
	calloc2D_double(&disps_init, particle_no, 3);
	ZA_allocateMemory();
	ZA_displacements(current_catalogue, &disps_init, NULL, Pk_current, false);
	apply_ZA_displacements(current_catalogue, &disps_init, REAL_SPACE);
	ZA_velocities(current_catalogue, &disps_init, cur_params);
	ZA_freeMemory();

	//measure Pk
	overdensity_allocateMemory();
	toRedshift(current_catalogue, cur_params);
	haloes_measure_Pk(current_catalogue, Pk_current_out, 1.0, CIC);
	overdensity_freeMemory();

	free2D_double(&disps_init, particle_no);	
	//free(current_catalogue->particles);
	//free(current_catalogue);
	

	//Proceeding with the measured Pk	
	Spline *Pk_current_ZA = input_spline_file(Pk_current_out, k_Pkmono_format, NORMAL, MEASURED);
	
	//printf("current kmin: %lf, kmax: %lf \n", Pk_current->splineInfo->xmin, Pk_current->splineInfo->xmax);		
	
	//for current simulation:		
	Spline **Pk_splines_zBins = malloc((size_t)z_bins * sizeof(*Pk_splines_zBins));
	Spline **variance_splines_zBins = malloc((size_t)z_bins * sizeof(*variance_splines_zBins));	

	//finding a redshift to start minimising from
	double z_to, z_guess, sq_diff, min_sq_diff;
	min_sq_diff = DBL_MAX;
	for (int i = 0; i < z_bins; i++) {	
		z_to = bin_to_x(z_binInfo, i);
		printf("preparing Pk at z: %lf \n", z_to);
		
		Pk_splines_zBins[i] = prep_Pk_constz(current_gravity, cur_params->z, z_to, Pk_current_ZA, k_binInfo, z_binInfo, cur_params);	
		variance_splines_zBins[i] = prep_variances(R_binInfo, Pk_splines_zBins[i]);
		
		//check more points
		sq_diff = pow(splint_Pk(Pk_splines_zBins[i], 0.1) - splint_Pk(Pk_target, 0.1), 2.0);
		if (sq_diff < min_sq_diff) {
			min_sq_diff = sq_diff;
			z_guess = z_to;
		}			
	}	
	printf("z guess: %lf \n", z_guess);
	
	//minimization
	dsq_multimin(vary_z_current, z_guess, variance_target, variance_splines_zBins, z_binInfo);	
	
	//print after z scaling	
	Spline *Pk_current_ZA_z_rescaled = prep_Pk_constz(current_gravity, cur_params->z, z_rescaled, Pk_current_ZA, k_binInfo, z_binInfo, cur_params);
	char Pk_current_z_scaled_out[100];
	sprintf(Pk_current_z_scaled_out, "%s/data/rescaling/model_to_model/Pk_current_z_scaled.dat", home_directory);
	print_spline_file(Pk_current_ZA_z_rescaled, k_binInfo, Pk_current_z_scaled_out);
	
	//model after z scaling
	char Pk_current_model_z_scaled_out[100];
	sprintf(Pk_current_model_z_scaled_out, "%s/data/rescaling/model_to_model/Pk_current_model_z_scaled.dat", home_directory);
	
	
	//scale box
	volume_limits[0] *= s;
	volume_limits[1] *= s;
	volume_limits[2] *= s;
	volume = volume_limits[0]*volume_limits[1]*volume_limits[2];
	
	//move particles, ensure PBC
	for (int i = 0; i < current_catalogue->particle_no; i++) {
		current_catalogue->particles[i].x *= s;
		current_catalogue->particles[i].y *= s;
		current_catalogue->particles[i].z *= s;
		//PBC(&(current_catalogue->particles[i].x), &(current_catalogue->particles[i].y), &(current_catalogue->particles[i].z));
	}	
	
	//scale masses, OMEGA M DEPENDS ON Z!!!!!!!!!!!
	for (int i = 0; i < particle_no; i++) {		
		current_catalogue->particles[i].mass *= pow(s, 3.0) * (omega_m_z(targ_params->z, targ_params)/omega_m_z(cur_params->z, cur_params));
	}
	
	//mass_min_max(current_catalogue);
	//printf("AFTER RESCALING. M min: %le, M max: %le\n", Mmin, Mmax);
	
	//scale velocities, redshift not yet changed
	scale_velocities(s, current_catalogue, cur_params, targ_params);	
	
	//after rescaling, redshift changed
	Parameters *cur_params_rescaled = malloc(sizeof(*cur_params));		
	cur_params_rescaled->H0 = 100.0*0.678;
	cur_params_rescaled->omega_m_0 = 0.15;
	cur_params_rescaled->omega_v_0 = 0.85;
	cur_params_rescaled->omega_r_0 = 0.0;
	cur_params_rescaled->omega = 1.0;	
	cur_params_rescaled->gamma = 0.545;
	cur_params_rescaled->z = z_rescaled;
	
	char Pk_current_s_scaled_out[100];
	sprintf(Pk_current_s_scaled_out, "%s/data/rescaling/model_to_model/Pk_current_ZA_s_scaled.dat", home_directory);
	overdensity_allocateMemory();
	haloes_measure_Pk(current_catalogue, outFile_scaled, 1.0, CIC);			
	overdensity_freeMemory();	
	
	Spline* Pk_current_ZA_s_rescaled = input_spline_file(outFile_scaled, k_Pkmono_format, NORMAL, MEASURED);	
	
	//current s scaled, measured
	//char current_outPath_s_measured[100];
	//sprintf(current_outPath_s_measured, "%s/data/rescaling/model_to_model/delsq_current_s_measured.dat", home_directory);
	//print_delsq(Pk_current_rescaled, current_outPath_s_measured, rescaling_k_bin_info, false);
	
	
	exit(0);
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
	
	//RESCALED PARAMS HERE
	//scale velocities. Current z is z after rescaling? Parameters are now those of the target cosmology ? (H specifically)
	//scale_velocities(s, current_catalogue, cur_params, targ_params, z_current, z_rescaled);
	toRedshift(current_catalogue, targ_params);
	
	double** individual_displacements;
	calloc2D_double(&individual_displacements, particle_no, 3);
	//dont forget to change to fractional displacements for ZA
	ZA_displacements(current_catalogue, &individual_displacements, Pk_current_rescaled, Pk_target, true);
	//mass_bias_displacements(current_catalogue, &individual_displacements);
	apply_ZA_displacements(current_catalogue, &individual_displacements, REAL_SPACE);
	ZA_velocities(current_catalogue, &individual_displacements, targ_params);
	
	
	free(Pk_splines_zBins);
	free(variance_splines_zBins);
	free(k_binInfo);
	free(R_binInfo);
	free(z_binInfo);
	//gsl_integration_workspace_free(OLV_workspace);	
	return 0;		
}




//CATALOGUE TO MODEL
/*
int rescale_catalogue_to_model() {	

	Parameters cur_params = {
		h = 0.678, //dimensionless hubble parameter
		H0 = 100.0*0.678, //[kms^-1 Mpc^-1]
		omega_m_0 = 0.24,
		omega_v_0 = 0.76,
		omega_r_0 = 0.0,
		omega = 1.0,	
		fg = pow(0.24, 0.545)	
	};	
	
	Parameters targ_params = {
		h = 0.678, //dimensionless hubble parameter
		H0 = 100.0*0.678, //[kms^-1 Mpc^-1]
		omega_m_0 = 0.24,
		omega_v_0 = 0.76,
		omega_r_0 = 0.0,
		omega = 1.0,	
		fg = pow(0.24, 0.4)	
	};	

	char Pk_target_path[100];
	sprintf(Pk_target_path, "%s/data/linear_matter_pk_sig8_0.593_z_0.75.dat", home_directory);

	char sim_current[100];	
	sprintf(sim_current, "%s/%s", mock_directory, "AllHalo_GR_1500Mpc_a0.7_1.txt");	
	volume_limits[0] = 1500.0;  // h^-1 Mpc 
	volume_limits[1] = 1500.0;
	volume_limits[2] = 1500.0;
	volume = volume_limits[0]*volume_limits[1]*volume_limits[2];
	
	char outFile_scaled[100];
	sprintf(outFile_scaled, "%s/data/rescaling/cat_to_model/measured_pk_scaled.dat", home_directory);
	char outFile[100];
	sprintf(outFile, "%s/data/rescaling/cat_to_model/measured_pk.dat", home_directory);
	
	z_current = a_to_z(0.7);
	z_target = 0.75;
	
	//current simulation data input
	Particle_Catalogue *current_catalogue = input_catalogue_file(sim_current, 1, mocks_format);	
	
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
	print_delsq(Pk_current_extended, current_outPath, rescaling_k_bin_info, false);
	
	//target
	char target_outPath[100];
	sprintf(target_outPath, "%s/data/rescaling/cat_to_model/delsq_target.dat", home_directory);
	print_delsq(Pk_target_extended, target_outPath, rescaling_k_bin_info, false);	
	
	//minimization
	dsq_multimin(vary_z_current, z_guess, variance_spline_target, variance_splines_zBins);
	int z_ind_final = x_to_bin(rescaling_z_bin_info, z_rescaled);
	BinInfo* rescaled_k_bins = prep_bins(k_min/s, k_max/s, k_bins, LOG_BIN);
	
	//current s scaled
	char current_outPath_s[100];
	sprintf(current_outPath_s, "%s/data/rescaling/cat_to_model/delsq_current_s.dat", home_directory);
	print_delsq(Pk_current_extended, current_outPath_s, rescaling_k_bin_info, true);
	
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
	//print_spline(Pk_current_rescaled, rescaling_k_bin_info);
	//exit(0);
	
	SplineInfo_Extended* Pk_current_rescaled_extended = extend_spline(Pk_current_rescaled);
	
	
	//current s scaled, measured
	char current_outPath_s_measured[100];
	sprintf(current_outPath_s_measured, "%s/data/rescaling/cat_to_model/delsq_current_s_measured.dat", home_directory);
	print_delsq(Pk_current_rescaled_extended, current_outPath_s_measured, rescaling_k_bin_info, false);
	
	
	
	
	exit(0);		
	
	
	//scale velocities. Current z is z after rescaling? Parameters are now those of the target cosmology ? (H specifically)
	scale_velocities(s, current_catalogue, &cur_params, &targ_params, z_current, z_rescaled);
	toRedshift(current_catalogue, &targ_params, z_rescaled);
	
	double** individual_displacements;
	calloc2D_double(&individual_displacements, particle_no, 3);
	//dont forget to change to fractional displacements for ZA
	//ZA_displacements(current_catalogue, &individual_displacements, Pk_target_extended, NGP);
	mass_bias_displacements(current_catalogue, &individual_displacements);
	apply_ZA_displacements(current_catalogue, &individual_displacements, REAL_SPACE);
	ZA_velocities(current_catalogue, &individual_displacements, &targ_params, z_rescaled);
	

	//calc_b_eff(z_rescaled);
	//bias_plot();
	
	

	
	//scale masses, velocities with s
	//does this come before or after HOD??
	//NEED SOME PARAMETERS HERE!!
	//scale_masses();
	//scale_velocities();
	
	
	//Pk_out(out_final);

	
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
*/


//CATALOGUE TO CATALOGUE
/*
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
	//R1_primed = m_to_R(Mmin, z_current);
	//R2_primed = m_to_R(Mmax, z_current);
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
	
	
	remove(tmpFile);
	return 0;		
}
*/


