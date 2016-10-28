//selects halos with mass above m_min, outputs to new file.
//NOT CURRENTLY USED
/*
int select_halos_mmin(char outFile[]) {
	double log_Mmin, Mmin;
	log_Mmin = 14.06;
	Mmin = pow(10., log_Mmin);	
	//catalogue should be read in
	FILE *f = fopen(outFile, "w");	
	for (int i = 0; i < particle_no; i++) {
		if (particles[i][6] > Mmin) {
			fprintf(f, mocks_format, particles[i][0], particles[i][1], particles[i][2],
			particles[i][3], particles[i][4], particles[i][5], particles[i][6]);
			fprintf(f, "\n");
		}
	}	
	fclose(f);
	return 0;
}
*/


//sets the global variables cdm_global and R_virial_global for use in all functions given a halo mass
//VARIANCE SHOULD BE SPLINED, STORED
int set_halo_params(double halo_mass, Spline *variance_spline, Parameters *params) {	

	Spline *variance_reversed = reverse_spline(variance_spline);
	double R_mstar = splint_generic(variance_reversed, pow(delta_c, 2.0));
	double Mstar = R_to_m(R_mstar, params);
		
	//del_nl_z() ????
	cdm_global = (c0 / (1.0 + params->z)) * pow(halo_mass / Mstar, beta); // eqn. B6
	R_virial_global = pow(3.0*halo_mass/(4.0*M_PI*rho_bar_z(params)*del_nl), 1.0/3.0); // eqn. B5
	rs_global = R_virial_global/cdm_global;	

	free(variance_reversed);
	return 0;
}
/*
int print_cdm(char outFile[]) {

	Parameters *cur_params_z0 = malloc(sizeof(*cur_params_z0));
	cur_params_z0->H0 = 67.8;
	cur_params_z0->omega_m_0 = 0.25;
	cur_params_z0->omega_v_0 = 0.75;
	cur_params_z0->omega_r_0 = 0.0;
	cur_params_z0->omega = 0.0;
	cur_params_z0->gamma = 0.545;
	cur_params_z0->z = 0.0;	
	
	Parameters *cur_params_z06 = malloc(sizeof(*cur_params_z06));
	cur_params_z06->H0 = 67.8;
	cur_params_z06->omega_m_0 = 0.25;
	cur_params_z06->omega_v_0 = 0.75;
	cur_params_z06->omega_r_0 = 0.0;
	cur_params_z06->omega = 0.0;
	cur_params_z06->gamma = 0.545;
	cur_params_z06->z = 0.6;

	char model1[100];
	sprintf(model1, "%s/data/models/matter_pk_om_m_025_z_0_sigma8_06.dat", home_directory);
	SplineInfo* Pk_model1_z06 = input_spline_file(model1, Pk_model_format, NORMAL);
	
	char model2[100];
	sprintf(model2, "%s/data/models/matter_pk_om_m_015_z_0_sigma8_1.dat", home_directory);
	SplineInfo* Pk_model2_z06 = input_spline_file(model2, Pk_model_format, NORMAL);
	
	Parameters *targ_params_z0 = malloc(sizeof(*targ_params_z0));
	targ_params_z0->H0 = 67.8;
	targ_params_z0->omega_m_0 = 0.15;
	targ_params_z0->omega_v_0 = 0.85;
	targ_params_z0->omega_r_0 = 0.0;
	targ_params_z0->omega = 0.0;
	targ_params_z0->gamma = 0.545;
	targ_params_z0->z = 0.0;	
	
	Parameters *targ_params_z06 = malloc(sizeof(*targ_params_z0));
	targ_params_z06->H0 = 67.8;
	targ_params_z06->omega_m_0 = 0.15;
	targ_params_z06->omega_v_0 = 0.85;
	targ_params_z06->omega_r_0 = 0.0;
	targ_params_z06->omega = 0.0;
	targ_params_z06->gamma = 0.545;
	targ_params_z06->z = 0.6;	
	
	
	//binning
	k_bins = 500;	
	Dplus_bins = 500;
	z_bins = 500;
	k_min = 0.01;
	k_max = 1.0;	
	z_min = 0.0;
	z_max = 3.0;
	double R_min = 0.01;
	double R_max = 5.0;
	rescaling_k_bin_info = prep_bins(k_min, k_max, k_bins, LOG_BIN);
	rescaling_R_bin_info = prep_bins(R_min, R_max, 500, LOG_BIN);
	rescaling_z_bin_info = prep_bins(z_min, z_max, z_bins, NORMAL_BIN);

	z_to_t_workspace = gsl_integration_workspace_alloc(z_to_t_workspace_size);
	OLV_workspace = gsl_integration_workspace_alloc(OLV_workspace_size);

	prep_redshift_spline(rescaling_z_bin_info);

	//at z = 1.2
	SplineInfo *Pk_model1_z12 = prep_Pk_constz(GR, 0.0, 0.6, Pk_model1_z06, rescaling_k_bin_info);	

	SplineInfo_Extended *Pk_model1_z06_ext = extend_spline_model(Pk_model1_z06);
	SplineInfo_Extended *Pk_model1_z12_ext = extend_spline_model(Pk_model1_z12);
	
	SplineInfo *Pk_model2_z12 = prep_Pk_constz(GR, 0.0, 0.6, Pk_model2_z06, rescaling_k_bin_info);	
	SplineInfo_Extended *Pk_model2_z06_ext = extend_spline_model(Pk_model2_z06);
	SplineInfo_Extended *Pk_model2_z12_ext = extend_spline_model(Pk_model2_z12);	

	
	SplineInfo *variance_model1_z06 = prep_variances(rescaling_R_bin_info, Pk_model1_z06_ext);
	SplineInfo *variance_model1_z12 = prep_variances(rescaling_R_bin_info, Pk_model1_z12_ext);
	
	SplineInfo *variance_model2_z06 = prep_variances(rescaling_R_bin_info, Pk_model2_z06_ext);
	SplineInfo *variance_model2_z12 = prep_variances(rescaling_R_bin_info, Pk_model2_z12_ext);
	//print_spline(variance_model2_z12, rescaling_R_bin_info);	
		
	Mmin = 5e12;
	Mmax = 1e15;
	int bins = 200;
	FILE *f = fopen(outFile, "w");
	double dm = (Mmax - Mmin)/(double)(bins-1);
	//model1 0.6, model2 0.6, model1 1.2, model2 1.2
	double m, cdm1, cdm2, cdm3, cdm4;
	for (int i = 0; i < bins; i++) {
		m = Mmin + (double)i * dm;

		set_halo_params(m, variance_model1_z06, cur_params_z0);
		cdm1 = cdm_global;		

		set_halo_params(m, variance_model2_z06, targ_params_z0);
		cdm2 = cdm_global;
		
		set_halo_params(m, variance_model1_z12, cur_params_z06);
		cdm3 = cdm_global;		

		set_halo_params(m, variance_model2_z12, targ_params_z06);
		cdm4 = cdm_global;
		
		fprintf(f, "%le \t %le \t %le \t %le \t %le \n", m, cdm1, cdm2, cdm3, cdm4);
	
	}
	fclose(f);
	return 0;

}
*/

//cumulative sum in rho NFW up to r, given a halo of mass m, r here in units of rs
double cumsum_nfw(double r) {	
	return (log(1.0 + r) - r*pow(1.0 + r, -1.0))*pow(log(1.0 + cdm_global) - cdm_global/(1.0 + cdm_global), -1.0);		
}

//outputs integral of rho_nfw * dV from 0 to r for multiple r going from 0 to r_virial
Spline* prep_nfw_cumsums(double halo_mass, int bins) {

	//set_halo_params(halo_mass);	
	
	double rstart, rstop, dr, r, cumsum;	
	
	double* Rs = malloc((size_t)bins * sizeof(*Rs));
	double* cumsums = malloc((size_t)bins * sizeof(*cumsums));
	
	rstart = 0.0;
	rstop = R_virial_global;
	dr = (rstop - rstart) / (double) (bins-1);
	
	//r will be in units of rs
	for (int i = 0; i < bins; i++) {
		r = rstart + (double)i * dr;
		r /= rs_global;
		cumsum = cumsum_nfw(r);
		Rs[i] = r;
		cumsums[i] = cumsum;
		if (cumsum < 0.0 || cumsum > 1.0) {
			printf("cumsum: %lf \n", cumsum);
			exit(0);
		}
	}
	
	return input_spline_values(bins, Rs, cumsums, MEASURED);
}

//prepares nfw cumulative sums for a given catalogue, places into array of SplineInfo objects
//MUST KNOW Mmin and Mmax beforehand, stored globally
int nfw_cumsums_catalogue(int bins, const char root_dir[]) {
	double dm, this_mass;
	char out_file[100];	
	dm = (Mmax - Mmin)/(double)(bins-1);
	for (int i = 0; i < bins; i++) {
		this_mass = Mmin + (double)i * dm;
		//printf("this mass is: %le\n", this_mass);
		sprintf(out_file, "%scumsums_mass_%d.dat", root_dir, i);
		//set the r_vir and the cdm parameters!!!
		//set_halo_params(this_mass);		
		//prep_nfw_cumsums(out_file);	
		//nfw_cumsums_reversed[i] = input_spline(out_file, "%le \t %le \n", true);
		remove(out_file);
	}
	return 0;
}

//for model power spectrum, COPIED FROM MIKE
double ukm_nfw_profile(double k){
    // For u(k|m) = (1/M)\int_vol \rho(x|m) exp^{-ikx}
    // Spherically symmetric profile, truncated at the virial radius. 
    // u(k|m) = \int_0^{r_{vir}} dr 4pi r^2 [sin(kr)/kr] rho(r|m)/m
    
    double Interim = sin(k*rs_global)*(gsl_sf_Si((1.0 + cdm_global)*k*rs_global) - gsl_sf_Si(k*rs_global)) - sin(cdm_global*k*rs_global)*pow((1.0 + cdm_global)*k*rs_global, -1.0) + cos(k*rs_global)*(gsl_sf_Ci((1.0 + cdm_global)*k*rs_global) - gsl_sf_Ci(k*rs_global));
    
    // for the NFW profile, defined by rhos, rs and c, mass is not independent. 
    Interim /= log(1.0 + cdm_global) - cdm_global/(1.0 + cdm_global);
    return Interim;
}

double halo_model_Pk(double k, double shot_noise, double twohalo_prefactor, int redshift_space, int order, Spline *Pk_spline, double fg) {
   
    double interim = shot_noise;  
  
	//one-halo term                                      				
	//interim *= pow(ukm_nfw_profile(k), 2.0);	
	
	//double onehalo_term = shot_noise * pow(ukm_nfw_profile(k), 2.0);
	double twohalo_term = splint_Pk(Pk_spline, k) * volume * exp(-k*k*R_nl*R_nl) * twohalo_prefactor;	 
	
	//2halo + redshift space
	if (redshift_space == REDSHIFT_SPACE) {
		if (order == 0) {
			
			//beta = fg because no galaxy bias assumed
			//onehalo_term += pow(ukm_nfw_profile(k), 2.0)*kaiserGauss_Monofactor(0.0, 0.0);				
			//twohalo_term *= kaiserGauss_Monofactor(0.0, fg); //k dependence only with velocity dispersion	
			twohalo_term *= KaiserLorentz_Monofactor(0.0, fg);
		
		} else if (order == 2) {	
		
			//onehalo_term += pow(ukm_nfw_profile(k), 2.0)*kaiserGauss_Quadfactor(0.0, 0.0);
			//twohalo_term *= kaiserGauss_Quadfactor(0.0, fg); //k dependece for lorentz/gaussian 
			twohalo_term *= KaiserLorentz_Quadfactor(0.0, fg);	
		
		} else if (order == 4) {
			// Lorentzian -> Gauss.  eval. for k=0 -> Kaiser factor.   
			//twohalo_term *= kaiserGauss_Hexfactor(0.0, fg); 
			twohalo_term *= KaiserLorentz_Hexfactor(0.0, fg);
		
		} else {
			printf("only monopole (order 0), quadrupole (order 2) and hexadecapole (order 4) available \n");
			exit(0);
		}
	} else if (redshift_space != REAL_SPACE) {
		printf("can only have REDSHIFT_SPACE or REAL_SPACE in halo model\n");
	}
	
	interim += twohalo_term;	
	//shotnoise * onehalo + twohalo	
	
	return interim; 	
}

//just outputs the halo model power spectrum to a given file
int haloModel_out(char out[], BinInfo *k_bins, double shot_noise, double twohalo_prefactor, int redshift_space, Spline *Pk_spline, Parameters *par) {
	
	double k, model_mono, model_quad, model_hex, fg;
	fg = pow(par->omega_m_0, par->gamma);	
	FILE *f = fopen(out, "w");
	if (out == NULL) {
		printf("wrong file in haloModel_out: %s \n", out);
		exit(0);
	}	

	for (int i = 0; i < k_bins->bins; i++) {
		k = bin_to_x(k_bins, i);		
		model_mono = halo_model_Pk(k, shot_noise, twohalo_prefactor, redshift_space, 0, Pk_spline, fg);
		model_quad = halo_model_Pk(k, shot_noise, twohalo_prefactor, redshift_space, 2, Pk_spline, fg);
		model_hex = halo_model_Pk(k, shot_noise, twohalo_prefactor, redshift_space, 4, Pk_spline, fg);
		//fprintf(f, "%le \t %le \t %le \n", k, model_mono, model_quad);
		fprintf(f, "%le \t %le \t %le \t %le \n", k, model_mono, model_quad, model_hex);
	}
	
	fclose(f);	
	return 0;
}

double I_integral_HOD(double t, void *params) {
	params = params;
	double f_t = log(1.0 + t) - t/(1.0 + t);
	return f_t/(pow(t,3.0)*pow(1.0 + t,2.0));
}


