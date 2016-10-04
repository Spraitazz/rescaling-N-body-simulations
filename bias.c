//http://arxiv.org/pdf/1308.5183v3.pdf
double f_v(double v) {
	double A, q, p;
	A = 0.216;
	q = 0.707;
	p = 0.3;
	return A*(1.0 + 1.0/pow(q*v*v, p))*exp(-q*v*v/2.0);
}

double ln_fv_deriv(double v) {
	double q = 0.707;
	double p = 0.3;
	return (1.0/(1.0 + 1.0/pow(q*v*v, p)))*(1.0/pow(q,p))*(-2.0*p)*pow(v, -2.0*p - 1.0) - q*v;
}

double m_to_v(double mass, double z) {
	double R = m_to_R(mass, z);
	double sigma_R = sqrt(splint_generic(variance_spline_current, R));
	double v = delta_c/sigma_R;
	return v;
}

double v_to_m(double v, double z) {
	double sigma_R = delta_c / v;
	double R = splint_generic(variance_spline_current_reversed, sigma_R*sigma_R);
	double mass = R_to_m(R, z);
	return mass;
}

//http://arxiv.org/pdf/1308.5183v3.pdf
double b_v(double v) {
	return 1.0 - 1.0/delta_c - (v/delta_c)*ln_fv_deriv(v);
}

double bf_over_m(double v, void * params) {
	double* z = (double*) params;	
	return b_v(v)*f_v(v)/v_to_m(v, z[0]);
}

double f_over_m(double v, void * params) {
	double* z = (double*) params;	
	return f_v(v)/v_to_m(v, z[0]);
}

int calc_b_eff(double z) {
	double vmin, vmax;
	vmin = m_to_v(Mmin, z);
	vmax = m_to_v(Mmax, z);
	
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
	double error_top, result_top, error_bot, result_bot;
	double params[1] = {z};		

	gsl_function F;	
	F.function = &bf_over_m;
	F.params = params;        

    gsl_integration_qags(&F, vmin, vmax, 0, 1e-7, 1000, w, &result_top, &error_top);		
	
	F.function = &f_over_m;
	
	gsl_integration_qags(&F, vmin, vmax, 0, 1e-7, 1000, w, &result_bot, &error_bot);
	
	b_eff = result_top / result_bot;
	printf("b_eff (effective mass bias for real catalogue ZA): %lf\n", b_eff);		
	
	gsl_integration_workspace_free(w);
	
	return 0;
}



int bias_plot(double z) {
	int bins = 100;
	double dm = (Mmax - Mmin) / (double)bins;
	double vmin;
	char out[] = "/home/jonas/Testing_GR/code/Jonas/output/rescaling/bias_plot.dat";
	FILE *f = fopen(out, "w");
	for (int i = 0; i < bins; i++) {
		Mmin += i*dm;
		vmin = m_to_v(Mmin, z);
		calc_b_eff(z);
		fprintf(f, "%le \t %le \n", vmin, b_eff);
	}
	fclose(f);
	return 0;
}

