SplineInfo_Extended* extend_spline_model(SplineInfo *spline) {
	
	SplineInfo_Extended *model_spline = malloc(sizeof(*model_spline));
	double pk_loA, pk_hiA, pk_lon, pk_hin, xmin, xmax, xend_min, xstart_max;
		
	xmin = spline->xmin;
	xmax = spline->xmax;
	xend_min = exp(log(xmin) + 0.2 * (log(xmax) - log(xmin))); //10% of range, can vary this
	xstart_max = xmax - 0.1 * (xmax - xmin);	
	
	powerlaw_regression(spline->lines, xmin, xend_min, 1.0, spline->x_vals, spline->y_vals, &pk_loA, &pk_lon);    
    powerlaw_regression(spline->lines, xstart_max, xmax, 1.0, spline->x_vals, spline->y_vals, &pk_hiA, &pk_hin);
    
    model_spline->spline = spline;
	model_spline->pk_hiA = pk_hiA;
	model_spline->pk_loA = pk_loA;
	model_spline->pk_hin = pk_hin;
	model_spline->pk_lon = pk_lon;
	model_spline->model = true;
	return model_spline;
}

SplineInfo_Extended* extend_spline(SplineInfo *spline) {
	SplineInfo_Extended* extended_spline = malloc(sizeof(*extended_spline));
	extended_spline->spline = spline;
	extended_spline->pk_hiA = 0.0;
	extended_spline->pk_loA = 0.0;
	extended_spline->pk_hin = 0.0;
	extended_spline->pk_lon = 0.0;
	extended_spline->model = false;
	return extended_spline;
}

//uses spline from what's available + power laws at low and high k. (all here return VP(k))
double splint_Pk_model(SplineInfo_Extended *pk_model, double k) {
	double Pk;
	if (pk_model->model == false) {
		printf("asking to splint a model, but giving a non-model extended SplineInfo\n");
		exit(0);
	} else {
		if (k < pk_model->spline->xmin) {
			Pk = pk_model->pk_loA*pow(k, 3.0 + pk_model->pk_lon);		
		} else if (k > pk_model->spline->xmax) {
			Pk = pk_model->pk_hiA*pow(k, 3.0 + pk_model->pk_hin);		
		} else {
			Pk = splint_generic(pk_model->spline, k);
		}
	}
	return Pk/volume;
}

double splint_Pk(SplineInfo_Extended *pk_spline, double k) {
	if (pk_spline->model) {
		return splint_Pk_model(pk_spline, k);
	} else {
		return splint_generic(pk_spline->spline, k)/volume;
	}
}

//volume term comes from VP(k) used for plotting
double delsq_lin(SplineInfo_Extended *pk_spline, double k) {
	double Pk = splint_Pk(pk_spline, k);	
	return volume * Pk * pow(k, 3.0)/(2.0*pi*pi);
}

double OLV(double k, void *params) {
	double margin, T_sq, R, delsq, result;
	
	OLV_parameters *parameters = (OLV_parameters*)params;	
	R = parameters->R;
	margin = 1e-6;	

	delsq = delsq_lin(parameters->spline, k);
	if (delsq < 1e-10) {
		//printf("delsq %lf at k %lf \n", delsq, k);
	}
		
	// if small, explicit taylor
	if (k*R < margin) {
		T_sq = 1.0 - pow(k*R,2.0)/5.0 + 3.0*pow(k*R, 4.0)/175.0; //+ O(kR^6)
	} else {
		T_sq = pow((3.0 / pow(k*R, 3.0))*(sin(k*R) - k*R*cos(k*R)), 2.0);
	}
	
	if (T_sq <= 0.0) {
		printf("T squared %lf at k %lf \n", T_sq, k);
		exit(0);
	}
	
	//additional 1/k factor from d(lnk) = 1/k dk	
	result = T_sq * delsq / k;	
	return result;	
}


//integrates the overdensity linear variance function at a given r from k = 0 to infinity, or arbitrary chosen upper limit
double integrate_OLV(double R, SplineInfo_Extended *pk_spline) {
	double error, result, margin;
	int errcode;	
	OLV_parameters params = {pk_spline, R};	
		
    gsl_function F;	
    F.function = &OLV;  
    F.params = &params; 
    margin = 1e-7;    

    gsl_error_handler_t* old_handler;
    old_handler = gsl_set_error_handler_off(); 
    
    if (OLV_workspace == NULL) {
    	printf("initialize OLV workspace \n");
    	exit(0);
    } else {
		errcode = gsl_integration_qagiu(&F, 0.0, 0, margin, OLV_workspace_size, OLV_workspace, &result, &error); 
	}
    
    switch (errcode) {
    	case GSL_EMAXITER:
    	
    		printf("the maximum number of subdivisions was exceeded. in OLV integral\n");
    		exit(0);
    		
    	case GSL_EROUND:   
    	 		
    		while (errcode == GSL_EROUND && margin < 1e-3) {
    			margin *= 2.0;
    			errcode = gsl_integration_qagiu(&F, 0.0, 0, margin, OLV_workspace_size, OLV_workspace, &result, &error); 		
    		}
    		
    		if (errcode == GSL_EROUND) {
    			printf("cannot reach tolerance because of roundoff error, or roundoff error was detected in the extrapolation table. Refuse to go below 1e-3 relative error, something is wrong. OLV integral\n");
    			exit(0);
    		}
    		
    		break;  		
    	case GSL_ESING:
    		return 0.0;
    		printf("a non-integrable singularity or other bad integrand behavior was found in the integration interval. OLV integral, R: %lf, result: %lf \n", R, result);
    		exit(0);
    		
    	case GSL_EDIVERGE:
    		return 0.0;
    		printf("the integral is divergent, or too slowly convergent to be integrated numerically. OLV integral, R: %lf, result: %lf \n", R, result);
    		exit(0);    
    		      
    }    
	return result;
}

//prepares OLV using the R binning info (upper limit, lower limit, bin no), and the power spectrum spline (functions.c)
SplineInfo* prep_variances(BinInfo *R_bin_info, SplineInfo_Extended *pk_spline) {

	double result, R;		

	double* Rs = malloc((size_t)R_bin_info->bins * sizeof(*Rs));
	double* OLVs = malloc((size_t)R_bin_info->bins * sizeof(*OLVs));
	int actual_data = 0;
	
	for (int i = 0; i < R_bin_info->bins; i++) {		
		R = bin_to_x(R_bin_info, i);				
		result = integrate_OLV(R, pk_spline);
		if (result > 1e-10) {
			Rs[i] = R;
			OLVs[i] = result;
			actual_data += 1;	
		} else {
			printf("OLV: %le at R: %lf \n", result, R);
		}		
	}	

	SplineInfo* toReturn = input_spline_values(actual_data, Rs, OLVs);
	return toReturn;	
}

//DELETE ME 
SplineInfo* prep_variances_test(BinInfo *R_bin_info, SplineInfo_Extended *pk_spline, double s_test) {	

	double result, R;

	double* Rs = malloc((size_t)R_bin_info->bins * sizeof(*Rs));
	double* OLVs = malloc((size_t)R_bin_info->bins * sizeof(*OLVs));
	int actual_data = 0;
	
	for (int i = 0; i < R_bin_info->bins; i++) {			
		R = bin_to_x(R_bin_info, i);	
		result = integrate_OLV(R/s_test, pk_spline);
		if (result > 1e-10) {
			Rs[i] = R;
			OLVs[i] = result;	
			actual_data += 1;
		}
	}
	
	SplineInfo* toReturn = input_spline_values(actual_data, Rs, OLVs);
	return toReturn;	
}
//DELETE ME

double OLV_smoothed(double k, void *params) {
	double R_nl_this, delsq;
		
	OLV_parameters *parameters = (OLV_parameters*)params;	
	R_nl_this = parameters->R;	
	
	delsq = delsq_lin(parameters->spline, k);
	
	return delsq * exp(-k*k*R_nl_this*R_nl_this) / pow(k, 3.0);
}

double expected_variance_smoothed(double R_nl_this, SplineInfo_Extended *pk_spline) {
	double error, result;
	OLV_parameters params = {pk_spline, R_nl_this}; 
	gsl_function F;	
    F.function = &OLV_smoothed;  
    F.params = &params;    
   	//double cellSize = volume_limits[0] / (double) cells[0];    
    //double k_nyquist = pi / cellSize;
    //gsl_integration_qags(&F, 2.0*pi/volume_limits[0], k_nyquist, 0, 1e-6, OLV_workspace_size, OLV_workspace, &result, &error); 
    gsl_integration_qagiu(&F, 2.0*pi/volume_limits[0], 0, 1e-6, OLV_workspace_size, OLV_workspace, &result, &error); 
	return result;
}

double delta_sq_int_func(double R, void *params) {
	double sigma, sigma_primed, s_cur, z_var;
	int index;	
	bool vary_z_cur;	
	
	Dsq_Params *parameters = (Dsq_Params*)params;
	s_cur = parameters->s;
	z_var = parameters->z_var;
	vary_z_cur = parameters->vary_z_current;	
	
	index = x_to_bin(rescaling_z_bin_info, z_var);

	if (vary_z_cur) {
		sigma = sqrt(splint_generic(parameters->variances_varz[index], R/s_cur));	
		sigma_primed = sqrt(splint_generic(parameters->variance_const, R));
	} else {
		sigma = sqrt(splint_generic(parameters->variance_const, R/s_cur));//at z, using current Pk for variance
		sigma_primed = sqrt(splint_generic(parameters->variances_varz[index], R)); //at z_primed, using target Pk
	}

	if (sigma_primed < 1e-10) {
		//is this the correct behaviour?
		return 1.0/R;
	} else {
		return pow(1.0 - sigma/sigma_primed, 2.0) / R;
	}
}

double delta_sq_rms(const gsl_vector *inputs, void *params) {

	double z_var, s_cur, error, result, R1_pr, R2_pr;	
	bool vary_z_cur;
		
	Multimin_Params *parameters = (Multimin_Params*) params;
	R1_pr = parameters->R1_primed;
	R2_pr = parameters->R2_primed;
	vary_z_cur = parameters->vary_z_current;
	
	s_cur = gsl_vector_get(inputs, 0);
	z_var = gsl_vector_get(inputs, 1);		
		
	//dont want negative z.outside of redshift range -> not good, return maximum possible value
	if (z_var < z_min || z_var > z_max) {
		return 1.0;
	} else {
		
		Dsq_Params* integral_params = malloc(sizeof(*integral_params));
		integral_params->s = s_cur; 
		integral_params->z_var = z_var;
		integral_params->vary_z_current = vary_z_cur;
		integral_params->variance_const = parameters->variance_const;
		integral_params->variances_varz = parameters->variances_varz;

		gsl_function F;	
		F.function = &delta_sq_int_func;
		F.params = integral_params;  	
		
		//cquad works much better here
		gsl_integration_cquad(&F, R1_pr, R2_pr, 0, 1e-6, dsq_workspace_cquad, &result, &error, NULL);		
		result /= log(R2_pr / R1_pr);	
		free(integral_params);		
		return result;
	}
}

int dsq_multimin(bool vary_z_cur, double z_init, SplineInfo* variance_const, SplineInfo** variances_varz) {
		
	size_t iter = 0;
	int status;
	int variables = 2;
	double size, nm_simplex_stop_size;
	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2; 
	gsl_multimin_fminimizer *dsq_minimizer = NULL;
	gsl_vector *inputs, *step_sizes;
	gsl_multimin_function dsq_func;
	
	//the size of the n (=2) dimensional simplex around a minimum for when the minimization has found "the minimum" (not necessarily the global one
	nm_simplex_stop_size = 1e-7;
	
	dsq_workspace_cquad = gsl_integration_cquad_workspace_alloc(500);
	
	//the parameters to send to the delta squared integral
	Multimin_Params* dsq_parameters = malloc(sizeof(*dsq_parameters));
	dsq_parameters->vary_z_current = vary_z_cur;
	dsq_parameters->R1_primed = R1_primed;
	dsq_parameters->R2_primed = R2_primed;
	dsq_parameters->variance_const = variance_const;
	dsq_parameters->variances_varz = variances_varz;
	
	//inputs are the parameters along which the function is minimized
	//param 1 - s, param 2 - z, initial guesses are 1 for scaling parameter, and z = z'	
	inputs = gsl_vector_alloc(variables);
	gsl_vector_set(inputs, 0, 1.0);
	gsl_vector_set(inputs, 1, z_init);	
	
	//param 1 - ds, param 2 - dz
	step_sizes = gsl_vector_alloc(variables);
	gsl_vector_set(step_sizes, 0, 0.02);
	gsl_vector_set(step_sizes, 1, rescaling_z_bin_info->dx);	

	//the function to minimize	
	dsq_func.n = variables;
	dsq_func.f = delta_sq_rms;
	dsq_func.params = dsq_parameters;	
	
	//allocate memory for minimizer
	dsq_minimizer = gsl_multimin_fminimizer_alloc(T, variables);

	//set minimizer variables
	gsl_multimin_fminimizer_set(dsq_minimizer, &dsq_func, inputs, step_sizes);

	//minimization loop
	do {
	
		iter++;
		
		//tries new parameters, evaluates function given
		status = gsl_multimin_fminimizer_iterate(dsq_minimizer);	
		if (status) break;	
		
		//check if minimum found - criteria for stopping is simplex size
		size = gsl_multimin_fminimizer_size(dsq_minimizer);
		status = gsl_multimin_test_size(size, nm_simplex_stop_size);
		
		if (iter % 20 == 0){
			printf("MINIMISATION. loop %ld. parameters, s: %lf, z: %lf, RMS diff in OLVs: %lf\n", iter, gsl_vector_get(dsq_minimizer->x, 0), gsl_vector_get(dsq_minimizer->x, 1), dsq_minimizer->fval);
			if (iter % 100 == 0) {
				nm_simplex_stop_size *= 10.0;
			}
		}
		
	} while (status == GSL_CONTINUE && iter < 200);	

	//final values of s, z, dsq after minimization
	s = gsl_vector_get(dsq_minimizer->x, 0);
	z_rescaled = gsl_vector_get(dsq_minimizer->x, 1);
	dsq_minimized_value = dsq_minimizer->fval;
	printf("s: %lf, z_rescaled: %lf, dsq at s, z: %le\n", s, z_rescaled, dsq_minimized_value);
	
	//free memory used
	gsl_multimin_fminimizer_free(dsq_minimizer);
	gsl_integration_cquad_workspace_free(dsq_workspace_cquad);		
	free(dsq_parameters);
	
	//can be used to check if everything went ok
	return status;
}

/*
//a is declared in main.c
int scale_velocities(double s, double H, double H_primed, double a_primed) {
	double fg, fg_primed;
	//fg = d lng/ d lna
	for (int i = 0; i < particle_no; i++) {
		particles[i][3] = particles[i][3]*s*(H_primed*a_primed*fg_primed)/(H*a*fg);
		particles[i][4] = particles[i][4]*s*(H_primed*a_primed*fg_primed)/(H*a*fg);
		particles[i][5] = particles[i][5]*s*(H_primed*a_primed*fg_primed)/(H*a*fg);	
	}
	return 0;
}
*/

/*
int prep_Dplus_kz(int gravity, BinInfo* zBins, BinInfo* kBins) {

	double k_cur;
	double zs[zBins->bins];
	double ks[kBins->bins];
	double Dpluses[Dplus_bins];
	
	for (int i = 0; i < kBins->bins; i++) {
		k_cur = bin_to_x(kBins, i);
		Dplus_calc(gravity, k_cur, &zs, &Dpluses);
	}
	
	
	const gsl_interp2d_type *T = gsl_interp2d_bilinear;
	gsl_spline2d *spline = gsl_spline2d_alloc(T, nx, ny);
	gsl_spline2d_init(spline, 
	
	int res = gsl_interp2d_eval_e()
	if (res == GSL_EDOM) {
	
	}
	
	return 0;
}
*/

SplineInfo* prep_Pk_constz(int gravity, double z_from, double z_to, SplineInfo *Pk_spline_in, BinInfo *k_bin_info) {
	if (z_from < 0.0 || z_to < 0.0) {
		printf("z < 0 in prep_Pk_constz. z from: %lf, z to: %lf \n", z_from, z_to);
		exit(0);
	}	
	double k_cur, Pk_cur, Dplus_cur, Dplus_after;
	SplineInfo* Dplus_tmp;		
	
	double* ks = malloc((size_t)k_bin_info->bins * sizeof(*ks)); //THESE CANNOT BE ARRAYS, OTHERWISE C HAPPENS
	double* Pks = malloc((size_t)k_bin_info->bins * sizeof(*Pks));

	if (gravity == GR) {
		//Dplus wont depend on k, spline gives z vs Dplus
		Dplus_tmp = Dplus_spline(GR, 0.0);
		Dplus_cur = splint_generic(Dplus_tmp, z_from);
		Dplus_after = splint_generic(Dplus_tmp, z_to);
		
		if (Dplus_cur <= 0.0 || Dplus_after <= 0.0) {
			printf("\nDplus <= 0. Dplus_cur = %lf, Dplus_after = %lf z = %lf. TEST DPLUS at z = 0.5: %lf \n", Dplus_cur, Dplus_after, z_from, splint_generic(Dplus_tmp, 0.5));
			exit(0);
		}
		
		//multiply by P(k), bin splines into z-bins			
		for (int j = 0; j < k_bin_info->bins; j++) {
			k_cur = bin_to_x(k_bin_info, j);
			Pk_cur = splint_generic(Pk_spline_in, k_cur);
			Pk_cur /= pow(Dplus_cur, 2.0);
			Pk_cur *= pow(Dplus_after, 2.0);
			ks[j] = k_cur;
			Pks[j] = Pk_cur;
		}
		
		//free(Dplus_tmp);			
	
	} else if (gravity == F4 || gravity == F5 || gravity == F6) {
		//Dplus depends on k		
		for (int i = 0; i < k_bin_info->bins; i++) {
			k_cur = bin_to_x(k_bin_info, i);
			
			Dplus_tmp = Dplus_spline(gravity, k_cur);
			Dplus_cur = splint_generic(Dplus_tmp, z_from);			
			Dplus_after = splint_generic(Dplus_tmp, z_to);
			
			if (Dplus_cur <= 0.0 || Dplus_after <= 0.0) {
				printf("\nDplus <= 0. Dplus_cur = %lf, Dplus_after = %lf z = %lf. TEST DPLUS at z = 0.5: %lf \n", Dplus_cur, Dplus_after, z_from, splint_generic(Dplus_tmp, 0.5));
				exit(0);
			}
			
			Pk_cur = splint_generic(Pk_spline_in, k_cur);
			Pk_cur /= pow(Dplus_cur, 2.0);
			Pk_cur *= pow(Dplus_after, 2.0);
			
			ks[i] = k_cur;
			Pks[i] = Pk_cur;
			
			free(Dplus_tmp);
		}	
	} else {
		printf("wrong gravity input in prep_pk_constz(), gravity: %d \n", gravity);
		exit(0);
	}	

	SplineInfo *toReturn = input_spline_values(k_bin_info->bins, ks, Pks);
	return toReturn;
}

int generate_sigma_plots() {
	int bins = 100;
	int z_ind_init = x_to_bin(rescaling_z_bin_info, z_current);
	int z_ind_scaled = x_to_bin(rescaling_z_bin_info, z_rescaled);
	printf("z ind init: %d, z ind scaled: %d \n", z_ind_init, z_ind_scaled);
	double sigma, sigma_resc_z, sigma_resc_both, sigma_primed, R;
	char plots_out[100];
	char smth[100];
	
	double dR = (R2_primed - R1_primed)/(double)(bins-1);
	
	sprintf(plots_out, "%s/data/rescaling/sigmas_TEST_NEW.dat", home_directory);
	FILE *f = fopen(plots_out, "w");

	if (vary_z_current) {
		sprintf(smth, "%s", "sigma CURRENT redshift rescaled, current length scaled");
	} else {
		sprintf(smth, "%s", "sigma TARGET redshift rescaled, current length scaled");
	}
	
	fprintf(f, "R \t sigma current \t sigma target \t sigma current z rescaled \t %s \n", smth); 
	for (int i = 0; i < bins; i++) {
		R = R1_primed + (double)i * dR;		
		sigma = sqrt(splint_generic(variance_splines_zBins[z_ind_init], R));
		sigma_primed = sqrt(splint_generic(variance_spline_target, R));
		//rescaling size specifically for current box
		sigma_resc_z = sqrt(splint_generic(variance_splines_zBins[z_ind_scaled], R));
		if (vary_z_current) {
			//also after size redshift apply rescaled z to current	
			sigma_resc_both = sqrt(splint_generic(variance_splines_zBins[z_ind_scaled], R/s));
		} else {		
			//size rescale on current, redshift rescale on target	
			//sigma_resc_redshift = sqrt(variance_R_z(R, z_rescaled, false));	
		}

		//sigma_resc = sigma_primed = 1.0;
		fprintf(f, "%le \t %le \t %le \t %le \t %le\n", R, sigma, sigma_primed, sigma_resc_z, sigma_resc_both);
	}
	fclose(f);
	return 0;
}

