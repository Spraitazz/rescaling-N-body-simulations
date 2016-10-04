/*
Mostly a file for playing around with gsl minimisation, ignore this.



*/

double noise_func_multimin(const gsl_vector *v, void *params) {
	double A, p, offsets, k, shotnoise, nbar, nbar_sats, nbar_full, mean_N_satellite, model_pk;
	double *parameters = (double *)params;
	nbar = (particle_no - satellite_no) / volume;
	satellite_no = 1000000;
	nbar_full = particle_no / volume;
	nbar_sats = satellite_no / volume;
	mean_N_satellite = pow((halo_mass - M0)/M1, alpha);

	//x could be halo mass, changes monopole, need to output list for diff masses + spline/splint
	A = gsl_vector_get(v, 0);
	p = gsl_vector_get(v, 1);
	
	//high-k values more important, penalise for getting them wrong??
	
	for (int i = 0; i < spectrum_size_final; i++) {				
		if (power_spectrum_monopole[i] != 0.0) {
			k = bin_k_average[i];
			shotnoise = A*(1.0/nbar_sats)*pow(ukm_nfw_profile(k), p);
			model_pk = haloModel_pk(k, nbar, beta, 0);

			offsets += pow(power_spectrum_monopole[i] - shotnoise - model_pk, 2.0);
		} 
	}
	//printf("offsets: %lf\n", offsets);
	return offsets;
}

//TRY WITH AVERAGES OF 20 SIMS
double noise_func(double shotnoise_param, void *params) {
	(void)(params); /* avoid unused parameter warning */
	double offsets, k, shotnoise;
	double nbar = (particle_no - satellite_no) / volume;
	double nbar_full = particle_no / volume;
	for (int i = 0; i < spectrum_size_final; i++) {				
		if (power_spectrum_monopole[i] != 0.0) {
			k = bin_k_average[i];
			//NEED TO INCOROPORATE UKM nfw profile
			offsets += pow((power_spectrum_monopole[i] - (shotnoise_param)/nbar_full) - haloModel_pk(k, nbar, beta, 0), 2.0);
		} 
	}
	return offsets;
}

//http://arxiv.org/pdf/1308.5183v3.pdf eq. 15 for f(v), then https://arxiv.org/pdf/astro-ph/0001493v1.pdf eq. 17
int onehalo_shotnoise() {
	double A, q, p, fv, v;
	A = 0.216;
	q = 0.707;
	p = 0.3;
	
	fv = A*(1.0 + 1.0/pow(q*v*v, p))*exp(-q*v*v/2.0);

	return 0;
}

typedef struct {
	double A_min;
	double p_min;
	double f_min;
} MinimiserOutput;

MinimiserOutput minimise_noise_multidim() {
	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2; //nmsimplex2 may be faster??
	gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc(T, 2);
	
	gsl_vector *x;
	gsl_vector *step_size;
	
	gsl_multimin_function my_func;	
	
	x = gsl_vector_alloc(2);
	//param 1 - A, param 2 - p
	gsl_vector_set(x, 0, randomDouble(0.0,1.0));
	gsl_vector_set(x, 1, randomDouble(-2.0,2.0));
	
	step_size = gsl_vector_alloc(2);
	gsl_vector_set(step_size, 0, 0.01);
	gsl_vector_set(step_size, 1, 0.1);
	
	//replace with actual parameters?
	double params[2] = {0.0, 0.0};
	
	my_func.n = 2;
	my_func.f = noise_func_multimin;
	my_func.params = params;
	
	gsl_multimin_fminimizer_set(s, &my_func, x, step_size);
	
	
	
	size_t iter = 0;
	int status;
	double size;
	
	do {
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);
		
		if (status) break;
		
		size = gsl_multimin_fminimizer_size(s);
		status = gsl_multimin_test_size(size, 1e-3);
		
		if (status == GSL_SUCCESS) {
			printf("found min\n");
		}
		
		printf ("%5lu %.5f %.5f %10.5f\n", iter,
              gsl_vector_get(s->x, 0), 
              gsl_vector_get(s->x, 1), 
              gsl_multimin_fminimizer_minimum(s));
	
	} while (status == GSL_CONTINUE && iter < 100);
	
	MinimiserOutput out = {
		.A_min = gsl_vector_get(s->x, 0),
		.p_min = gsl_vector_get(s->x, 1),
		.f_min =  gsl_multimin_fminimizer_minimum(s)	
	};
	
	//printf("A: %lf, p: %lf, minvalue: %lf\n", gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1), gsl_multimin_fminimizer_minimum(s));

	return out;
}

int test_multimin() {
	power_spectrum_monopole = malloc(26 * sizeof(*power_spectrum_monopole));
	bin_k_average = malloc(26 * sizeof(*bin_k_average));
	FILE *pk_avgs = fopen("./plots/toaverage/randoms_average.dat", "r");
	for (int i = 0; i < 26; i++) {
		ph = fscanf(pk_avgs, "%lf \t %lf \n", &bin_k_average[i], &power_spectrum_monopole[i]);
	}
	fclose(pk_avgs);

	spectrum_size_final = 26;
	MinimiserOutput global_min, this_min;
	global_min = minimise_noise_multidim();
	
	for (int i = 0; i < 2000; i++) {
		this_min = minimise_noise_multidim();
		if (this_min.f_min < global_min.f_min) {
			global_min = this_min;
		}
	}
	
	printf("A min: %lf, p min: %lf, f min: %lf\n", global_min.A_min, global_min.p_min, global_min.f_min);
	
	return 0;

}

int minimise_noise() {
	int status;
	double param_lower, param_upper, param_min;
	const gsl_min_fminimizer_type *T = gsl_min_fminimizer_brent;
	gsl_min_fminimizer *s = gsl_min_fminimizer_alloc(T);
	
	gsl_function F;
	F.function = &noise_func;
	F.params = 0;
	
	param_lower = 0.0;
	param_upper = 1.5;
	param_min = 1.1;
	
	gsl_min_fminimizer_set(s, &F, param_min, param_lower, param_upper);
	
	size_t iter = 0;
	
	do {
		iter++;
		status = gsl_min_fminimizer_iterate(s);
		printf("%d\n", status);
		param_min = gsl_min_fminimizer_x_minimum(s);
		param_lower = gsl_min_fminimizer_x_lower(s);
		param_upper = gsl_min_fminimizer_x_upper(s);
		printf("between [%lf, %lf] param min: %lf\n", param_lower, param_upper, param_min);
		
		status = gsl_min_test_interval(param_lower, param_upper, 1000000.0, 0.01);
		
		if (status == GSL_SUCCESS) {
			printf("converged to parameter minimum: %lf\n", param_min);
		}
		
	} while (status == GSL_CONTINUE && iter < 1000);
	
	return 0;
}


int minimisation_main() {
    test_multimin();
	
	FILE *f = fopen("/home/jonas/Testing_GR/code/Jonas/output/mintest/testpk.dat", "w");
	double A = 0.633;
	double p = -0.455;
	double k;
	//double satellite_no = 1000000;
	double nbar_sats = satellite_no / volume;
	//double mean_N_satellite = pow((halo_mass - M0)/M1, alpha);
	double shotnoise;
	for (int i = 0; i < 26; i++) {
		k =  bin_k_average[i];
		shotnoise = A*(1.0/nbar_sats)*pow(ukm_nfw_profile(k), p);
		fprintf(f, "%lf \t %lf \n", k, power_spectrum_monopole[i] - shotnoise);
	
	}
	fclose(f);
	
	return 0;
}
