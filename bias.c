//http://arxiv.org/pdf/1308.5183v3.pdf
double f_v(double v) {
	double A, q, p;
	A = 0.216;
	q = 0.707;
	p = 0.3;
	return A*(1.0 + 1.0/pow(q*v*v, p))*exp(-q*v*v/2.0);
}

double f_v_new(double v) {
	return 0.0;
}

double ln_fv_deriv(double v) {
	double q = 0.707;
	double p = 0.3;
	return (1.0/(1.0 + 1.0/pow(q*v*v, p)))*(1.0/pow(q,p))*(-2.0*p)*pow(v, -2.0*p - 1.0) - q*v;
}

double m_to_v(double mass, Parameters *params, Spline* variance_spline) {
	double R = m_to_R(mass, params);
	double sigma_R = sqrt(splint_generic(variance_spline, R));
	double v = delta_c/sigma_R;
	return v;
}

double v_to_m(double v, Bias_Params *bias_params) {
	double sigma_R = delta_c / v;
	Spline* variance_reversed = reverse_spline(bias_params->variance);
	double R = splint_generic(variance_reversed, sigma_R*sigma_R);
	double mass = R_to_m(R, bias_params->parameters);
	return mass;
}

//http://arxiv.org/pdf/1308.5183v3.pdf
double b_v(double v) {
	return 1.0 - 1.0/delta_c - (v/delta_c)*ln_fv_deriv(v);
}

double b_m(double m, Parameters *params, Spline* variance_spline) {
	return b_v(m_to_v(m, params, variance_spline));
}

double bf_over_m(double v, void * params) {
	Bias_Params *parameters = (Bias_Params*) params;
	return b_v(v)*f_v(v)/v_to_m(v, parameters);
}

double f_over_m(double v, void * params) {
	Bias_Params *parameters = (Bias_Params*) params;
	return f_v(v)/v_to_m(v, parameters);
}

double calc_b_eff(Parameters *params, Spline* variance_spline, double Mmin, double Mmax) {
	double vmin, vmax;
	vmin = m_to_v(Mmin, params, variance_spline);
	vmax = m_to_v(Mmax, params, variance_spline);
	
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
	double error_top, result_top, error_bot, result_bot;
	
	Bias_Params *send = malloc(sizeof(*send));
	send->parameters = params;
	send->variance = variance_spline;

	gsl_function F;	
	F.function = &bf_over_m;
	F.params = send;        

    gsl_integration_qags(&F, vmin, vmax, 0, 1e-7, 1000, w, &result_top, &error_top);		
	
	F.function = &f_over_m;
	
	gsl_integration_qags(&F, vmin, vmax, 0, 1e-7, 1000, w, &result_bot, &error_bot);
	

	double b_eff = result_top / result_bot;
	//printf("b_eff (effective mass bias for real catalogue ZA): %lf\n", b_eff);		
	
	free(send);
	gsl_integration_workspace_free(w);
	
	return b_eff;
}



int bias_plot(double z) {
	int bins = 100;
	double dm = (Mmax - Mmin) / (double)bins;
	double vmin;
	char out[] = "/home/jonas/Testing_GR/code/Jonas/output/rescaling/bias_plot.dat";
	FILE *f = fopen(out, "w");
	for (int i = 0; i < bins; i++) {
		Mmin += i*dm;
		//vmin = m_to_v(Mmin, z);
		//calc_b_eff(z, SplineInfo* variance_spline);
		fprintf(f, "%le \t %le \n", vmin, b_eff);
	}
	fclose(f);
	return 0;
}

