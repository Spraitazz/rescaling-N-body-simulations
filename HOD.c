double rho_box, halo_mass;


int prep_bin_info() {
	//binning
	k_bins = 500;	
	Dplus_bins = 500;
	//k_min = Pk_current->xmin;
	//k_max = Pk_current->xmax;	
	k_min = 0.01;
	k_max = 1.5;
	z_min = 0.0;
	z_max = 3.5;
	z_bins = 300;
	R_bins = 300;
	R1_primed = 0.05;
	R2_primed = 5.0;
	rescaling_k_bin_info = prep_bins(k_min, k_max, k_bins, LOG_BIN); 
	rescaling_R_bin_info = prep_bins(R1_primed, R2_primed, R_bins, LOG_BIN);
	rescaling_z_bin_info = prep_bins(z_min, z_max, z_bins, NORMAL_BIN);
	return 0;
}

int HOD_main() {
	
	// model spline for Pk
	char Pk_in_model[100];
	sprintf(Pk_in_model, "%s/data/linear_matter_pk_sig8_0.593_z_0.75.dat", home_directory);
	Pk_current = input_spline_file(Pk_in_model, Pk_model_format, NORMAL);	
	
	char tmpFile[100];
	char tmpFile2[100];
	char tmpFile3[100];

	char Pk1[100];
	char Pk2[100];
	char Pk3[100];
	
	//PARAMETERS
	//redshift_space = false;
	z_current_global = 0.75;

	prep_bin_info();	
	
	//redshift 0.75 -> 3.5
	//z_to_t_workspace = gsl_integration_workspace_alloc(z_to_t_workspace_size);
	//prep_redshift_spline(&rescaling_z_bin_info);
	//prep_Pk_constz(GR, 0.75, z_current_global, &Pk_current, &rescaling_k_bin_info);	
	Pk_current_extended = extend_spline_model(Pk_current);
	//gsl_integration_workspace_free(z_to_t_workspace);
	
	OLV_workspace = gsl_integration_workspace_alloc(OLV_workspace_size);
	variance_spline_current = prep_variances(rescaling_R_bin_info, Pk_current_extended);
	variance_spline_current_reversed = reverse_spline(variance_spline_current);
	gsl_integration_workspace_free(OLV_workspace);
	
	R_nl = splint_generic(variance_spline_current_reversed, 1.0);
	if (R_nl < 1e-10) {
		printf("0 R_nl?? Rmin: %lf, Rmax: %lf, sigma^2 min: %le, sigma^2 max: %le \n", splint_generic(variance_spline_current_reversed, variance_spline_current_reversed->xmin), splint_generic(variance_spline_current_reversed, variance_spline_current_reversed->xmax), variance_spline_current_reversed->xmin, variance_spline_current_reversed->xmax);
		exit(0);
	} else {
		printf("\n Rnl is: %le \n", R_nl);	
	}
	
	sigma_exp_smoothed = sqrt(expected_variance_smoothed(R_nl, Pk_current_extended));
	printf("\n expected standard deviation of smoothed displacement: %lf \n", sigma_exp_smoothed);
	
	halo_no = 25000;
	particle_no = halo_no;
	prep_oneHalo();
	
	for (int i = 0; i < 100; i++) {	

		sprintf(Pk1, "%s/data/HOD_rands/onehalo/measured_pk_folded1_%d.dat", home_directory, i);
		sprintf(Pk2, "%s/data/HOD_rands/onehalo/measured_pk_folded4_%d.dat", home_directory, i);
		sprintf(Pk3, "%s/data/HOD_rands/onehalo/measured_pk_folded16_%d.dat", home_directory, i);

		printf("\n ------------------------------------------------------- \n RUN %d \n ------------------------------------------------------- \n", i);			

		//need to reset this before every run, as changes in Pk measurement 
		particle_no = halo_no;		
		initParticles();	
		
		populateHaloes(halo_mass);	

		add_centrals();			
		add_satellites();
		//increases particle_no	-> central_no + satellite_no
		prep_HOD_Pk(true);
		
		printf("\n Measuring Pk. particles: %d, centrals: %d, satellites: %d. Output to %s \n", particle_no, central_no, satellite_no, Pk1);	
		overdensity_allocateMemory();
		haloes_measure_Pk(Pk1, 1.0, NGP);
		haloes_measure_Pk(Pk2, 4.0, NGP);
		haloes_measure_Pk(Pk3, 16.0, NGP);		
		overdensity_freeMemory();
	
		//!!
		free2D_double(&particles, particle_no);
		
	}
	
	weighted_shotnoise_sats /= 100.0;		
	//weighted_shotnoise_sats = weighted_shotnoise();	
	
	char out_model[100];
	sprintf(out_model, "%s/data/HOD_rands/model_pk.dat", home_directory);
	haloModel_out(out_model, rescaling_k_bin_info, weighted_shotnoise_sats, REAL_SPACE, Pk_current_extended);
	//HOD_model(out_model, 0.02, 1.15, 200);	
	
	
	return 0;
}

//sets some numbers
int prep_oneHalo() {

	rho_box = rho_bar_z(z_current_global);
	halo_mass = rho_box * volume / halo_no;
	//halo_mass = 5e14;
	printf("\n halo mass: %le \n", halo_mass);	

	mean_N_sats = pow((halo_mass - M0)/M1, alpha);	
	
	//getting Mstar for cdm, using the overdensity linear variance obtained from the model power spectrum
	R_mstar = splint_generic(variance_spline_current_reversed, pow(delta_c, 2.0));
	
	//m* is the nonlinear mass scale at z = 0 defined such that sigma^2 (m*) = delta_c ^ 2
	//inverse spline returns R when sigma^2 has this value, then get m*
	Mstar = R_to_m(R_mstar, z_current_global);
	
	nfw_spline_reverse = reverse_spline(prep_nfw_cumsums(halo_mass, 1000));
	
	printf("\n m*: %le, cdm: %le, Rvir: %le, rs: %le \n", Mstar, cdm_global, R_virial_global, rs_global);	
 
	return 0;
}

//adds centrals according to HOD paper
int add_centrals() {
	double mean_N_central, this_halo_mass;
	int has_central, total_centrals;

	//find out how many centrals and allocate them
	total_centrals = 0;
	int* have_centrals = malloc((size_t)halo_no * sizeof(*centrals));
	
	for (int i = 0; i < halo_no; i++) {
		this_halo_mass = particles[i][6];
		mean_N_central = 0.5 * (1.0 + gsl_sf_erf((log10(this_halo_mass) - log_Mmin)/sigma_logm));		
		has_central = (int) gsl_ran_bernoulli(gsl_ran_r, mean_N_central);
		has_central = 1;
		
		// boolean list over haloes. 
		have_centrals[i] = has_central;

		if (has_central == 1) {			
			total_centrals += 1;
		}
	}
	
	central_no = total_centrals;	
	total_centrals = 0;
	calloc2D_double(&centrals, central_no, 7);
	
	for (int i = 0; i < halo_no; i++) {
		if (have_centrals[i] == 1) {
			centrals[total_centrals][0] = particles[i][0];
			centrals[total_centrals][1] = particles[i][1]; 
			centrals[total_centrals][2] = particles[i][2]; 
			centrals[total_centrals][3] = particles[i][3]; 
			centrals[total_centrals][4] = particles[i][4]; 
			centrals[total_centrals][5] = particles[i][5]; 
			centrals[total_centrals][6] = particles[i][6];
			total_centrals += 1;  
		}
	}
	
	free(have_centrals);
	return 0;
}

//each central Pk is weighted by the number of satellites assigned to its halo
double weighted_shotnoise() {

	double sum_weights = 0.0;
	double sum_weights_sq = 0.0;
	int no_tmp;
	for (int i = 0; i < halo_no; i++) {
		no_tmp = gsl_ran_poisson(gsl_ran_r, mean_N_sats); 
		sum_weights += (double) no_tmp;
		sum_weights_sq += (double) (no_tmp * no_tmp);		
	}
	return volume * (sum_weights_sq / pow(sum_weights, 2.0));
}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      

int add_satellites() {

	double cos_theta_sat, theta_sat, phi_sat, r_sat, unif_rnd, mean_r_sat;
	double x_sat, y_sat, z_sat;
	//velocity assignment for satellites variables
	double result_I, error, rho_nfw, f_cdm, sigma_sq_sat, sigma_sq_sat_tcdm, vx, vy, vz, sigma;
	
	//how many satellites for each halo? HOD paper
	int* no_satellites = calloc(halo_no, sizeof(*no_satellites));	 
	satellite_no = 0;
	double sum_weights = 0.0;
	double sum_weights_sq = 0.0;	
	for (int i = 0; i < halo_no; i++) {
		// e.g. 5 satellites per central. doesn't have to be stochastic. 
		no_satellites[i] = gsl_ran_poisson(gsl_ran_r, mean_N_sats); // take this out for the time being.
		sum_weights += (double) no_satellites[i];
		sum_weights_sq += (double) (no_satellites[i] * no_satellites[i]);
		//no_satellites[i] = 4;
		satellite_no += no_satellites[i];
	}
	
	double SN = sum_weights_sq / pow(sum_weights, 2.0);
	//printf("\n shotnoise here: %lf \n", SN*volume);
	weighted_shotnoise_sats += SN*volume;
	
	if (satellite_no > 0) {	
		calloc2D_double(&satellites, satellite_no, 7);	
		printf("\n adding %d satellites\n", satellite_no);	
	
		//prepare integral for velocity assignment of satellites
		//gsl_function F;
	  	//F.function = &I_integral_HOD;
	  	//F.params = 0;   
	  	//end velocity assignment stuff					

		int satellite_index = 0;
		mean_r_sat = 0.0;
		for (int i = 0; i < halo_no; i++) {
			//each central i will have no_satellites[i] as found before
		
		
			for (int j = 0; j < no_satellites[i]; j++) {
				//picking values uniformly randomly distributed on a sphere
				cos_theta_sat = randomDouble(-1.0, 1.0);
				theta_sat = acos(cos_theta_sat);
				phi_sat = randomDouble(0.0, 2.0*pi);	
			
				r_sat = splint_generic(nfw_spline_reverse, randomDouble(0.0, 1.0));
				//this is in units of rs
				r_sat *= rs_global;			
				//track the mean radial displacement
				mean_r_sat += r_sat;				
		
				//actual position is with respect to the i-th central
				x_sat = r_sat * sin(theta_sat) * cos(phi_sat) + particles[i][0];
				y_sat = r_sat * sin(theta_sat) * sin(phi_sat) + particles[i][1];
				z_sat = r_sat * cos(theta_sat) + particles[i][2];		
		
				//PBC
				PBC(&x_sat, &y_sat, &z_sat);
		
				//place the satellite in particle array		
				satellites[satellite_index][0] = x_sat;
				satellites[satellite_index][1] = y_sat;
				satellites[satellite_index][2] = z_sat;
			
				//finally, assign velocities	
			
				//gsl_integration_qags(&F, r_sat/rs_global, 1000.0, 0, 1e-5, 1000, onehalo_workspace, &result_I, &error); 
	
				//f_cdm = log(1.0 + cdm_global) - cdm_global/(1.0 + cdm_global);	
				//sigma_sq_sat = (G*halo_mass*cdm_global*cdm_global*r_sat/(R_virial_global*R_virial_global*f_cdm))*pow(1.0 + cdm_global*r_sat/R_virial_global, 2.0)*result_I;
			
				// nothing to do with tcdm. mead et al sigma^2(M).
				sigma_sq_sat_tcdm = (G*halo_mass/(3.0*R_virial_global))*cdm_global*(1.0 - 1.0/(1.0+pow(cdm_global,2.0)) - 2.0*log(1.0 + cdm_global)/(1.0 + cdm_global))/(2.0*pow(log(1.0 + cdm_global) - cdm_global/(1.0 + cdm_global),2.0));
			
				//double ratio = sigma_sq_sat_tcdm/sigma_sq_sat;
				//if (ratio > 1.1 || ratio < 0.9) {
					//printf("ratio: %lf\n", ratio);
			//	}
	
				//this is the relative velocity to the halo, so add together with the host halo velocity
				sigma = sqrt(sigma_sq_sat_tcdm);
				vx = gsl_ran_gaussian(gsl_ran_r, sigma_sq_sat);
				vy = gsl_ran_gaussian(gsl_ran_r, sigma_sq_sat);
				vz = gsl_ran_gaussian(gsl_ran_r, sigma_sq_sat);	
			
				satellites[satellite_index][3] = particles[i][3] + vx;
				satellites[satellite_index][4] = particles[i][4] + vy;
				satellites[satellite_index][5] = particles[i][5] + vz;
				satellites[satellite_index][6] = particles[i][6];		
			
				satellite_index += 1;
			}		
		}
	
		//sanity check?
		mean_r_sat /= (double) satellite_no;
		printf("\n mean r of sats: %lf \n", mean_r_sat); 
	
		free(no_satellites);
		
	} else {
		printf("no satellites to be added for this mass halo \n");
	}
	return 0;
}

int prep_HOD_Pk(bool sats) {

	//for measuring either centrals only, or centrals + sats, allocate the right amount of memory
	particle_no = central_no;
	if (sats) particle_no += satellite_no;	
	particles = realloc(particles, (size_t)particle_no * sizeof(*particles));
	
	for (int i = 0; i < particle_no; i++) {
		particles[i] = malloc(7 * sizeof(*particles[i]));
	}
	
	for (int i = 0; i < central_no; i++) {
		memcpy(particles[i], centrals[i], 7 * sizeof(*particles[i]));		
	}
	//after measuring satellites + centrals Pk, everything is freed before next loop of assignments
	if (sats) free2D_double(&centrals, central_no);
	
	if (sats) {
		for (int i = 0; i < satellite_no; i++) {
			memcpy(particles[i + central_no], satellites[i], 7 * sizeof(*particles[i + central_no]));	
		}
		free2D_double(&satellites, satellite_no);	
	}
	
	for (int i = 0; i < particle_no; i++) {
		if (particles[i][0] == 0.0 || particles[i][1] == 0.0 || particles[i][2] == 0.0) {
			printf("zero. ind: %d, x: %lf, y: %lf, z: %lf \n", i, particles[i][0], particles[i][1], particles[i][2]);
		}
	}
	
	return 0;
}
