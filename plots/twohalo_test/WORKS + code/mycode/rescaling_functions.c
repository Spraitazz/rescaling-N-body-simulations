/*
double Dplus_k_z(double k, double z, bool z_const) {
	int Dplus_cur_ind;
	double Dplus_cur, logkmin, logkmax, logdk, logk;
	
	logkmin = log(Dplus_k_min);
	logkmax = log(Dplus_k_max);	
	logdk = (logkmax - logkmin) / (double)Dplus_k_bins;	
	
	if (z_const) {
		if (k > Dplus_const_z_spline.xmax) {
			Dplus_cur = splint_generic(Dplus_const_z_spline, Dplus_const_z_spline.xmax);		
		} else if (k < Dplus_const_z_spline.xmin)  {
			Dplus_cur = splint_generic(Dplus_const_z_spline, Dplus_const_z_spline.xmin);		
		} else {
			Dplus_cur = splint_generic(Dplus_const_z_spline, k);
		}
	} else {
		//printf("trying to get D(k,z)\n\n");
		//find correct k-bin
		//COULD BE USEFUL TO USE INFO FROM 2 NEARST K BINS TOO
		logk = log(k);
		Dplus_cur_ind = (int) floor((logk - logkmin) / logdk);		
		if (Dplus_cur_ind >= Dplus_k_bins) {
			//this definitely happens
			Dplus_cur_ind = Dplus_k_bins - 1;
			//printf("Wtf k bins at dplus(k, z), logk %lf, logkmin %lf, logdk %lf\n", logk, logkmin, logdk);
		} else if (Dplus_cur_ind < 0) {
			//printf("same wtf %d, k %lf, logk %lf, logkmin %lf, logdk %lf\n", Dplus_cur_ind, k, logk, logkmin, logdk);
			Dplus_cur_ind = 0;
		}
	
		//select the Dplus value at the z value given
		if (z < Dplus_z_min) {
			//Dplus_cur = splint_generic(Dplus_splines[Dplus_cur_ind], Dplus_splines[Dplus_cur_ind].xmin);
		} else if (z > Dplus_z_max) {
			//Dplus_cur = splint_generic(Dplus_splines[Dplus_cur_ind], Dplus_splines[Dplus_cur_ind].xmax);
		} else {
			//Dplus_cur = splint_generic(Dplus_splines[Dplus_cur_ind], z);
		}		
	}	
	return Dplus_cur;
}
*/

//volume term comes from VP(k) used for plotting
double delsq_lin(SplineInfo_Extended *pk_spline, double k) {
	double Pk = splint_Pk(pk_spline, k);
	return volume * Pk * pow(k, 3.0)/(2.0*pi*pi);
}

double OLV(double k, void *params) {
	double margin, T_sq, Dplus_cur, R, z, delsq;	
	SplineInfo_Extended *pk_spline;	
	
	OLV_parameters *parameters = (OLV_parameters*)params;	
	R = parameters->R;
	pk_spline = parameters->spline;		
	margin = pow(10.0,-6.0);	

	delsq = delsq_lin(pk_spline, k);
		
	// if small, explicit taylor
	if (k*R < margin) {
		T_sq = 1.0 - pow(k*R,2.0)/5.0 + 3.0*pow(k*R, 4.0)/175.0; //+ O(kR^6)
	} else {
		T_sq = pow((3.0 / pow(k*R, 3.0))*(sin(k*R) - k*R*cos(k*R)), 2.0);
	}
	//additional 1/k factor from d(lnk) = 1/k dk
	return T_sq * delsq / k;	
}


//integrates the overdensity linear variance function at a given r from k = 0.0001 to k = 12.0
double integrate_OLV(double R, SplineInfo_Extended *pk_spline) {
	double error, result, margin, mm, cc;	
	OLV_parameters params = {pk_spline, R}; 	
    gsl_function F;	
    F.function = &OLV;  
    F.params = &params;    
    gsl_integration_qagiu(&F, 0.0, 0, 1e-5, 2000, OLV_workspace, &result, &error);    
	return result;
}

int prep_variances(BinInfo *R_bin_info, char OLV_out[], SplineInfo_Extended *pk_spline) {	
	FILE* outFile = fopen(OLV_out, "w");
	double result, R;	
	OLV_workspace = gsl_integration_workspace_alloc(2000);	
	for (int i = 0; i < R_bin_info->bins; i++) {		
		R = bin_to_x(R_bin_info, i);		
		result = integrate_OLV(R, pk_spline);			
		fprintf(outFile, "%lf \t %lf \n", R, result);		
	}
	gsl_integration_workspace_free(OLV_workspace);
	fclose(outFile);
	return 0;	
}

double OLV_smoothed(double k, void *params) {
	double R_nl, Pk, delsq;
	SplineInfo_Extended *pk_spline;	
		
	OLV_parameters *parameters = (OLV_parameters*)params;	
	R_nl = parameters->R;
	pk_spline = parameters->spline;
	
	delsq = delsq_lin(pk_spline, k);
	
	return delsq * exp(-k*k*R_nl*R_nl) / pow(k, 3.0);
}

double expected_variance_smoothed(double R_nl, SplineInfo_Extended *pk_spline) {
	double error, result;
	OLV_workspace = gsl_integration_workspace_alloc(2000);
	OLV_parameters params = {pk_spline, R_nl}; 
	gsl_function F;	
    F.function = &OLV_smoothed;  
    F.params = &params;    
    gsl_integration_qagiu(&F, 2.0*pi/volume_limits[0], 0, 1e-5, 2000, OLV_workspace, &result, &error); 	
    //gsl_integration_qags(&F, 2.0*pi/volume_limits[0], 1.0, 0, 1e-5, 2000, OLV_workspace, &result, &error); 	
    
	gsl_integration_workspace_free(OLV_workspace);
	return result;
}

double variance_R_z(double R, double z, bool current) {
	double result, z_ind, dz_this;
	int index;

	if ((current && vary_z_current) || ((!current) && (!vary_z_current))) {
		//using the data with varied z
		if (z > z_max) {
			index = z_bins-1;
		} else if (z < z_min) {
			index = 0;
		} else {
			index = (int)floor((z - z_min)/dz);		
			z_ind = (double)index * dz + z_min;
			dz_this = fabs((z - z_ind)/dz);						
			result = (1.0 - dz_this)*splint_generic(&variance_splines[index], R);
			if (z > z_ind && index < z_bins - 1) {			
				result += dz_this*splint_generic(&variance_splines[index+1], R);			
			} else if (z < z_ind && index > 0) {
				result += dz_this*splint_generic(&variance_splines[index-1], R);
			} 
		}
		if (result == 0.0) {
			//printf("sigma 0. index: %d, variance at ind: %le, R: %lf, z_ind: %lf, z: %lf, dz: %lf\n", index, splint_generic(variance_splines[index], R), R, z_ind, z, dz);
		}
	} else {
		//z is constant
		if (current) {
			result = splint_generic(&variance_spline_current, R);
		} else {
			result = splint_generic(&variance_spline_target, R);
		}
	}
	return result;
}

//CHANGE TO HAVE OLV SPLINES FOR DIFFERENT Z, Z PRIME
//z is set in main.c!!!
double delta_sq_int_func(double R, void *params) {
	double sigma, sigma_primed, s, z_const, z_var, vary_param;
	int index;	
	double *parameters = (double *)params;
	s = parameters[0];
	z_var = parameters[1];
	z_const = parameters[2];
	vary_param = parameters[3];	

	if (vary_param > 0.0) {
		//varying z in current simulation
		sigma = sqrt(variance_R_z(R/s, z_var, true));	
		sigma_primed = sqrt(variance_R_z(R, z_const, false));
	} else {
		//varying z in target simulation
		sigma = sqrt(variance_R_z(R/s, z_const, true));//at z, using current Pk for variance
		sigma_primed = sqrt(variance_R_z(R, z_var, false)); //at z_primed, using target Pk
	}

	if (sigma_primed == 0.0) {
		//printf("sigma primed is 0. vary param: %lf, z_var: %lf, z_const: %lf\n", vary_param, z_var, z_const);
		//is this the correct behaviour?
		return 1.0/R;
	} else {
		return pow(1.0 - (sigma / sigma_primed), 2.0) / R;
	}
}

double delta_sq_rms(const gsl_vector *inputs, void *params) {
	//dont want negative z
	if (gsl_vector_get(inputs, 1) < z_min) {
		return 1.0;
	} else {	
		double R1_primed, R2_primed, z_const, vary_param, error, result;		
		double *parameters = (double*)params;
		R1_primed = parameters[0];
		R2_primed = parameters[1];
		z_const = parameters[2];
		vary_param = parameters[3];

		
		// QAGP integration??
		/*
		double pts[OLV_z_bins + 2];
		pts[0] = Dplus_z_min;
		for (int i = 0; i < Dplus_z_bins; i++) {
			pts[i+1] = Dplus_z_min + Dplus_dz * (double)(i+1);
		}
		pts[Dplus_z_bins+1] = Dplus_z_max;
		*/

		//[0] - s, [1] - z (varied), [2] - z (constant) [3] - sets whether the varied z is in target or current sim
		double integral_params[4] = {gsl_vector_get(inputs, 0), gsl_vector_get(inputs, 1), z_const, vary_param};

		gsl_function F;	
		F.function = &delta_sq_int_func;
		F.params = integral_params;  	
		//gsl_integration_qags(&F, R1_primed, R2_primed, 0, 1e-4, 1000, dsq_workspace, &result, &error);
		gsl_integration_cquad(&F, R1_primed, R2_primed, 0, 1e-6, dsq_workspace_cquad, &result, &error, NULL);		
		result /= log(R2_primed / R1_primed);	
		//printf("Result of integral: %lf\n", result);
		
		return result;
	}
}

int dsq_multimin() {		
	size_t iter = 0;
	int status;
	double size, vary_param, z_const, z_var_init;
	dsq_workspace_cquad = gsl_integration_cquad_workspace_alloc(200);
	
	if (vary_z_current) {
		//providing z_primed, varying z in current sim
		vary_param = 1.0;
		z_const = z_target;
		z_var_init = z_current;
	} else {
		//providing z, varying z_primed (target z)
		vary_param = -1.0;
		z_const = z_current;
		z_var_init = z_target;
	}

	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2; 
	gsl_multimin_fminimizer *dsq_minimizer = gsl_multimin_fminimizer_alloc(T, 2); //using 2 variables (dimensions)
	
	//z here is the known one, NOT VARIED
	double dsq_parameters[4] = {R1_primed, R2_primed, z_const, vary_param};
	
	//inputs are the parameters along which the function is minimized
	gsl_vector *inputs;
	gsl_vector *step_sizes;	
	inputs = gsl_vector_alloc(2);
	//param 1 - s, param 2 - z, initial guesses are 1 for scaling parameter, and z = z'
	gsl_vector_set(inputs, 0, 1.0);
	gsl_vector_set(inputs, 1, z_var_init);	
	//param 1 - ds, param 2 - dz
	step_sizes = gsl_vector_alloc(2);
	gsl_vector_set(step_sizes, 0, 0.01);
	gsl_vector_set(step_sizes, 1, dz);	

	gsl_multimin_function dsq_func;
	dsq_func.n = 2;
	dsq_func.f = delta_sq_rms;
	dsq_func.params = dsq_parameters;		
	gsl_multimin_fminimizer_set(dsq_minimizer, &dsq_func, inputs, step_sizes);

	do {
		iter++;
		status = gsl_multimin_fminimizer_iterate(dsq_minimizer);	
		if (status) break;	
		size = gsl_multimin_fminimizer_size(dsq_minimizer);
		status = gsl_multimin_test_size(size, 1e-4);
	
		printf("MINIMISATION parameters, s: %lf, z: %lf, function: %lf\n", gsl_vector_get(dsq_minimizer->x, 0), gsl_vector_get(dsq_minimizer->x, 1), dsq_minimizer->fval);
	} while (status == GSL_CONTINUE && iter < 100);	

	s = gsl_vector_get(dsq_minimizer->x, 0);
	z_rescaled = gsl_vector_get(dsq_minimizer->x, 1);
	dsq_minimized_value = dsq_minimizer->fval;
	printf("s: %lf, z_rescaled: %lf, dsq at s, z: %le\n", s, z_rescaled, dsq_minimized_value);
	gsl_integration_cquad_workspace_free(dsq_workspace_cquad);	
	return 0;
}

int rescale_box() {
	return 0;
}

int rescale_redshift() {
	return 0;
}

int prep_Pk_constz(int gravity, double z_fixed, SplineInfo *Pk_spline_in, BinInfo *k_bin_info) {
	if (z_fixed < 0.0) {
		printf("z < 0 \n");
		return 1;
	}	
	double k_cur, Pk_cur, Dplus_cur;
	
	//redshift-time relation reversed:
	z_to_t_workspace = gsl_integration_workspace_alloc(1000);
	prep_redshift_spline(z_min-0.1, z_max+0.1, 1000);	
	
	//z is constant
	int bins = k_bin_info->bins;
	double* ks = malloc((size_t)bins * sizeof(*ks)); //THESE CANNOT BE ARRAYS, OTHERWISE C HAPPENS
	double* Pks = malloc((size_t)bins * sizeof(*Pks));
	if (gravity == GR) {
		//Dplus wont depend on k, spline gives z vs Dplus
		SplineInfo tmp_GR_Dplus = Dplus_spline(GR, 0.0);
		Dplus_cur = splint_generic(&tmp_GR_Dplus, z_fixed);
		if (Dplus_cur <= 0.0) {
			printf("\n Dplus = %lf, z = %lf. TEST DPLUS at z = 0.5: %lf \n", Dplus_cur, z_fixed, splint_generic(&tmp_GR_Dplus, 0.5));
			return 1;
		}
		
		//multiply by P(k), bin splines into z-bins			
		for (int j = 0; j < bins; j++) {
			k_cur = bin_to_x(k_bin_info, j);
			Pk_cur = pow(Dplus_cur, 2.0) * splint_generic(Pk_spline_in, k_cur);
			ks[j] = k_cur;
			Pks[j] = Pk_cur;
		}			

		*Pk_spline_in = input_spline_values(bins, ks, Pks);		
	
	} else {
		//Dplus depends on k		
		for (int i = 0; i < bins; i++) {
			k_cur = bin_to_x(k_bin_info, i);
			SplineInfo tmp = Dplus_spline(gravity, k_cur);
			Dplus_cur = splint_generic(&tmp, z_fixed);
			Pk_cur = splint_generic(Pk_spline_in, k_cur);
			ks[i] = k_cur;
			Pks[i] = pow(Dplus_cur, 2.0) * Pk_cur;
		}
			
		*Pk_spline_in = input_spline_values(bins, ks, Pks);			
				
	}	
	gsl_integration_workspace_free(z_to_t_workspace);
	return 0;
}


int prep_Pk_varz(int gravity, SplineInfo *Pk_spline_in, SplineInfo **Pk_splines_zBins, BinInfo *k_bin_info, BinInfo *z_bin_info) {
	double k_cur, z_cur, Pk_cur, Dplus_cur;
	
	//redshift-time relation reversed:
	z_to_t_workspace = gsl_integration_workspace_alloc(1000);
	prep_redshift_spline(z_min, z_max + 0.1, 1000);	

	if (Pk_splines_zBins == NULL) {
		printf("if varying z, provide the address for an array of SplineInfo (&SplineInfo*) \n");
		return 1;
	}

	int z_bins = z_bin_info->bins;
	if (gravity == GR) {
		SplineInfo tmp_GR_Dplus = Dplus_spline(GR, 0.0);						
		for (int i = 0; i < z_bins; i++) {
			double ks[z_bins];
			double Pks[z_bins];
			z_cur = bin_to_x(z_bin_info, i);
			for (int j = 0; j < k_bin_info->bins; j++) {
				k_cur = bin_to_x(k_bin_info, j);
				Pk_cur = pow(splint_generic(&tmp_GR_Dplus, z_cur), 2.0) * splint_generic(Pk_spline_in, k_cur);
				//have P(k) vs k in this z-bin
				ks[j] = k_cur;
				Pks[j] = Pk_cur;					
			}					
			(*Pk_splines_zBins)[i] = input_spline_values(k_bin_info->bins, ks, Pks);			
		}	
	} else {		
		//first bin Dpluses by k
		SplineInfo* Dplus_splines_kBins = malloc((size_t)k_bin_info->bins * sizeof(*Dplus_splines_kBins));
		for (int j = 0; j < k_bin_info->bins; j++) {
			k_cur = bin_to_x(k_bin_info, j);
			Dplus_splines_kBins[j] = Dplus_spline(gravity, k_cur);					
		}
		
		//now bin P(k) by z
		for (int i = 0; i < z_bins; i++) {
			double ks[z_bins];
			double Pks[z_bins];
			z_cur = bin_to_x(z_bin_info, i);
			for (int j = 0; j < k_bin_info->bins; j++) {
				k_cur = bin_to_x(k_bin_info, j);					
				Dplus_cur = splint_generic(&Dplus_splines_kBins[j], z_cur);
				if (Dplus_cur <= 0.0) {
					printf("Dplus is %le in prep_Pk_varz() \n", Dplus_cur);
				}
				Pk_cur = pow(Dplus_cur, 2.0) * splint_generic(Pk_spline_in, k_cur);
				//have P(k) vs k in this z-bin
				ks[j] = k_cur;
				Pks[j] = Pk_cur;					
			}					
			(*Pk_splines_zBins)[i] = input_spline_values(k_bin_info->bins, ks, Pks);			
		}	
	}
	gsl_integration_workspace_free(z_to_t_workspace);
	return 0;
}


//prepare [bins] splines, each at different z, for the varied-z simulation (can be either target or current)
//prepare one spline for the other simulation at the given z
/*
int prep_variance_splines(SplineInfo_Extended *Pk_spline_in, SplineInfo *variance_spline_out, SplineInfo_Extended* *Pk_splines_zBins, SplineInfo* *variance_splines_zBins, bool vary_z, double z_fixed) {
	double dz, z_cur;
	char temp_out[100];
	sprintf(temp_out, "%s%s", tempdir, "variance_spline_temp.dat");	
	
	if (vary_z) {
		if (Pk_splines_zBins == NULL || variance_splines_zBins == NULL) {
			printf("z varied but spline arrays for either input Pk or output variances not provided\n");
			return 1;
		}		
		for (int i = 0; i < rescaling_z_bin_info.bins; i++) {
			z_cur = bin_to_x(rescaling_z_bin_info, i);
			prep_variances(rescaling_R_bin_info, temp_out, Pk_splines_zBins[i]);
			(*variance_splines_zBins)[i] = input_spline_file(temp_out, "%le \t %le \n", false);
		}	
	} else {
		if (z_fixed == 0) {
			printf("if z is fixed, provide a value \n");
			return 1;
		}
		prep_variances(rescaling_R_bin_info, temp_out, Pk_spline_in);
		*variance_spline_out = input_spline_file(temp_out, "%le \t %le \n", false);	
	}

	remove(temp_out);
	return 0;
}
*/

int generate_sigma_plots() {
	int bins = 100;
	double sigma, sigma_resc_size, sigma_resc_redshift, sigma_primed, R;
	double dR = (R2_primed - R1_primed)/(double)(bins-1);
	char plots_out[] = "/home/jonas/Testing_GR/code/Jonas/output/rescaling/sigmas_TEST_NEW.dat";
	FILE *f = fopen(plots_out, "w");
	char smth[100];
	if (vary_z_current) {
		sprintf(smth, "%s", "sigma CURRENT redshift rescaled");
	} else {
		sprintf(smth, "%s", "sigma TARGET redshift rescaled");
	}
	fprintf(f, "R \t sigma current \t sigma target \t sigma current size rescaled \t %s \n", smth); 
	for (int i = 0; i < bins; i++) {
		R = R1_primed + (double)i * dR;		
		sigma = sqrt(variance_R_z(R, z_current, true));
		sigma_primed = sqrt(variance_R_z(R, z_target, false));
		//rescaling size specifically for current box
		sigma_resc_size = sqrt(variance_R_z(R/s, z_current, true));
		if (vary_z_current) {
			//also after size redshift apply rescaled z to current	
			sigma_resc_redshift = sqrt(variance_R_z(R/s, z_rescaled, true));
		} else {		
			//size rescale on current, redshift rescale on target	
			sigma_resc_redshift = sqrt(variance_R_z(R, z_rescaled, false));	
		}

		//sigma_resc = sigma_primed = 1.0;
		fprintf(f, "%le \t %le \t %le \t %le \t %le\n", R, sigma, sigma_primed, sigma_resc_size, sigma_resc_redshift);
	}
	fclose(f);
	return 0;
}

