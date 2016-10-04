double z_to_a(double z) {
	return 1.0/(z + 1.0); 
}

double a_to_z(double a) {
	return 1.0/a - 1.0;
}

double H_a(double a) {
	return H0*sqrt(omega_v + omega_m*pow(a,-3.0) + omega_r*pow(a,-4.0) + (1.0 - omega)*pow(a,-2.0));	
}

double z_to_t_int(double z, void * params) {
	params = params;	
	return 1.0/((1.0 + z)*H_a(z_to_a(z)));
}

double z_to_t(double z) {
	double error, result;
	gsl_function F;	
	F.function = &z_to_t_int;
	F.params = 0;  	
	gsl_integration_qagiu(&F, z, 0, 1e-7, z_to_t_workspace_size, z_to_t_workspace, &result, &error);	
	return result;	
}

double t_to_z(double t) {
	return splint_generic(redshift_spline_reversed, t);
}

double a_to_t(double a) {
	return z_to_t(a_to_z(a));
}

int prep_redshift_spline(BinInfo *z_bin_info) {
	redshift_spline_reversed = prep_spline_generic(z_bin_info, REVERSED, &z_to_t);
	
	if (splint_generic(redshift_spline_reversed, redshift_spline_reversed->xmin) != z_bin_info->xmax) {
		printf("redshift spline didn't go well \n");
		printf("%lf \t %lf \n", splint_generic(redshift_spline_reversed, redshift_spline_reversed->xmin), z_bin_info->xmax);
		exit(0);
	}
	
	if (splint_generic(redshift_spline_reversed, redshift_spline_reversed->xmax) != z_bin_info->xmin) {
		printf("redshift spline didn't go well \n");
		printf("%lf \t %lf \n", splint_generic(redshift_spline_reversed, redshift_spline_reversed->xmax), z_bin_info->xmin);
		exit(0);
	}
	
	return 0;
}

/*
double t_to_a_analytic(double t) {
	double t_lambda = 2.0/(3.0*H0*sqrt(omega_v));
	return pow(omega_m/omega_v, 1.0/3.0)*pow(sinh(t/t_lambda),2.0/3.0);
}
*/

//http://arxiv.org/pdf/1412.5195v3.pdf (19)
//f[0] - Dplus dot, f[1] - Dplus double dot
//y[0] - Dplus evaluated at some value of a, y[1] - Dplus dot, evaluated at some a
int growth_function(double t, const double Dplus_vals[], double dDplus_dt[], void *params) {
	double lambda_zero, lambda, k, omega_m_z, bracket_factor, a_cur, z_cur, t0;
	double R_bar, R0_bar;
	double* parameters = (double*) params;
	lambda_zero = parameters[0];
	k = parameters[1];
	t0 = parameters[2];
	
	//generic
	z_cur = t_to_z(t);	
	a_cur = z_to_a(z_cur);
	
	//omega_m + omega_v
	//a_cur = t_to_a_analytic(t);
	//z_cur = a_to_z(a_cur);
	
	//EdS
	//a_cur = pow(t/t0, 2.0/3.0);
	
	H = H_a(a_cur);
	omega_m_z = omega_m*pow(a_cur, -3.0)*pow(H0/H, 2.0);	
	
	//(7), (17)
	R0_bar = 3.0*H0*H0*(omega_m + 4.0*omega_v); 
	R_bar = 3.0*H0*H0*(omega_m*pow(a_cur, -3.0) + 4.0*omega_v);
	lambda = sqrt(pow(lambda_zero,2.0)*pow(R0_bar/R_bar,3.0));	
	bracket_factor = lambda*lambda*k*k/(a_cur*a_cur);
	
	dDplus_dt[0] = Dplus_vals[1]; //the value of the first derivative 
	dDplus_dt[1] = (3.0/2.0)*H*H*omega_m_z*(1.0 + (1.0/3.0)*(bracket_factor)/(1.0 + bracket_factor))*Dplus_vals[0] - 2.0*H*Dplus_vals[1]; //the value of the second derivative
	return GSL_SUCCESS;
}



//int Dplus_calc(int gravity, double k, double* zs, double* Dpluses)
int Dplus_calc(int gravity, double k, double** zs, double** Dpluses) {
	double t_start, t_end, Dplus_a_start, Dplusdot_a_start, dt, t_curr, t_next;
	double z_start, z_end, t0, lambda_zero, z_cur;
	int status;
	
	z_start = z_max;
	z_end = z_min;
	t_start = z_to_t(z_start);	
	t_end = z_to_t(z_end);
	t0 = t_start;	
	dt = (t_end - t_start)/(double)(Dplus_bins - 1);	
	
	//starting at a redshift where EdS applies, Dplus = (t/t0)^(2/3)
	Dplus_a_start = 1.0;
	Dplusdot_a_start = 2.0/(3.0*t0);
	
	switch (gravity) {
		case F4:
			lambda_zero = lambda_F4;
			break;
		case F5:
			lambda_zero = lambda_F5;
			break;
		case F6:
			lambda_zero = lambda_F6;
			break;
		case GR:
			lambda_zero = lambda_GR;
			break;
		default:
			printf("Gravity can only be one of the 4 (integers): F4, F5, F6, GR\n");
			return 1;	
	}		
	
	double current_Dplus_values[2] = {Dplus_a_start, Dplusdot_a_start};
	double parameters[] = {lambda_zero, k, t0};	
		
	gsl_odeiv2_system sys = {growth_function, NULL, 2, &parameters}; //rk8pd
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, dt, 1e-7, 0.0);		

	t_curr = t_start;
	(*zs)[0] = z_start;
	(*Dpluses)[0] = current_Dplus_values[0];			
	for (int i = 1; i < Dplus_bins; i++) {	
		//the next iteration goes from t_curr to t_next in some number of steps chosen by the driver
		t_next = t_curr + dt;			
		status = gsl_odeiv2_driver_apply(d, &t_curr, t_next, current_Dplus_values);
		z_cur = t_to_z(t_curr);	
		(*zs)[i] = z_cur;
		(*Dpluses)[i] = current_Dplus_values[0];
		
		if (status != GSL_SUCCESS) {
			printf ("error, return value=%d\n", status);
			break;
		}	
	}	
	gsl_odeiv2_driver_free(d);
	return 0;
}

SplineInfo* Dplus_spline(int gravity, double k) {

	//freeing these arrays (pointers) gives bad behaviour afterwards. Initialize elsewhere?
	double* zs = malloc((size_t)Dplus_bins * sizeof(*zs));
	double* Dpluses = malloc((size_t)Dplus_bins * sizeof(*Dpluses));

	Dplus_calc(gravity, k, &zs, &Dpluses);
	//starts from high z, goes to z = 0
	for (int i = 0; i < Dplus_bins; i++) {
	
		if (Dpluses[i] <= 0.0) {
			printf("Dplus %lf at index %d \n", Dpluses[i], i);
			exit(0);
		}
		
		//normalise to 1 at z = 0
		Dpluses[i] /= Dpluses[Dplus_bins - 1];
		
		if (zs[i] == 0.0 && Dpluses[i] != 1.0) {
			printf("Dplus not normalized correctly. ind: %d, z: %lf, Dplus: %lf \n", i, zs[i], Dpluses[i]);
			exit(0);
		}
	}
	SplineInfo* toReturn = input_spline_values(Dplus_bins, zs, Dpluses);	
	return toReturn;
}

/*
int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params) {
	double lambda_zero, k, omega_m_z, bracket_factor, a_cur, z_cur;
	double R_bar, R0_bar, lambda, fR0;
	double* parameters = (double*) params;
	lambda_zero = parameters[0];
	k = parameters[1];	
	
	//z_cur = splint_generic(redshift_spline_reversed, t);
	//a_cur = z_to_a(z_cur);
	a_cur = t_to_a_analytic(t);
	z_cur = a_to_z(a_cur);
	
	//(7), (17)
	R0_bar = 3.0*H0*H0*(omega_m + 4.0*omega_v);
	R_bar = 3.0*H0*H0*(omega_m*pow(a_cur, -3.0) + 4.0*omega_v);
	lambda = sqrt(pow(lambda_zero,2.0)*pow(R0_bar/R_bar,3.0));
	
	omega_m_z = omega_m*pow(a_cur, -3.0)*pow(H0/H, 2.0);
	omega_m_z = (1.0 - pow(H/H0, 2.0))/(1.0 - pow(a_cur, -3.0));
	H = H_z(z_cur);
	bracket_factor = lambda*lambda*k*k/(a_cur*a_cur);
	//printf("bracket: %le, H: %le, omega_m_Z: %le\n", bracket_factor, H, omega_m_z);
	
	gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 2, 2);
	gsl_matrix * m = &dfdy_mat.matrix; 
	gsl_matrix_set(m, 0, 0, 0.0);
	gsl_matrix_set(m, 0, 1, 1.0);
	gsl_matrix_set(m, 1, 0, (3.0/2.0)*H*H*omega_m_z*(1.0 + (1.0/3.0)*(bracket_factor)/(1.0 + bracket_factor)));
	gsl_matrix_set(m, 1, 1, -2.0*H);
	dfdt[0] = 0.0;
    dfdt[1] = 0.0;
	return GSL_SUCCESS;
}
*/




