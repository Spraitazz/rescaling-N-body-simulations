int initParticles() {
	particles = realloc(particles, (size_t)particle_no*sizeof(*particles));
	for(int i = 0; i < particle_no; i++) {
		particles[i] = (double *) calloc(7, sizeof(**particles));
	}	
	return 0;
}

int init_rng() {
	gsl_rng_env_setup();
	gsl_ran_T = gsl_rng_default;
	gsl_ran_r = gsl_rng_alloc(gsl_ran_T);
	return 0;
}

int clearParticles() {
	for (int i = 0; i < particle_no; i++) {
		particles[i][0] = 0.0;
		particles[i][1] = 0.0;
		particles[i][2] = 0.0;
		particles[i][3] = 0.0;
		particles[i][4] = 0.0;
		particles[i][5] = 0.0;
		particles[i][6] = 0.0;
	}
	return 0;
}

void homeDirChk() {
	char tmp[200];
	sprintf(tmp, "%s/t", home_directory);
	FILE *f = fopen(tmp, "w");
	if (f == NULL) {
		printf("check home directory setting, currently (%s) \n", home_directory);
		exit(0);
	} else {
		fclose(f);
		remove(tmp);
	}
}

int initGrid() {
	grid = calloc((size_t)(cells[0]*cells[1]*cells[2]), sizeof(*grid));	
	return 0;
}

int clearGrid() {
	for (int i = 0; i < cells[0]*cells[1]*cells[2]; i++) {
		grid[i] = 0.0;
	}
	return 0;
}

int overdensity_allocateMemory() {
	//allocate
	overdensity = (fftw_complex*) fftw_malloc((size_t)(cells[0]*cells[1]*cells[2])*sizeof(fftw_complex));
	overdensity_fourier = (fftw_complex*) fftw_malloc((size_t)(cells[0]*cells[1]*cells[2])*sizeof(fftw_complex));
	//zero
	overdensity_clearMemory();
	//set plan
	overdensity_plan = fftw_plan_dft_3d(cells[0], cells[1], cells[2], overdensity, overdensity_fourier, FFTW_FORWARD, FFTW_ESTIMATE);
	return 0;
}

int overdensity_clearMemory() {
	for (int i = 0; i < cells[0]*cells[1]*cells[2]; i++) {
		overdensity_fourier[i][0] = 0.0;
		overdensity_fourier[i][1] = 0.0;
		overdensity[i][0] = 0.0;
		overdensity[i][1] = 0.0;
	}
	return 0;
}

int overdensity_freeMemory() {
	fftw_free(overdensity);
	fftw_free(overdensity_fourier);
	return 0;
}

int ZA_allocateMemory() {
	displacements_fourier = (fftw_complex*)fftw_malloc((size_t)(cells_displ[0]*cells_displ[1]*cells_displ[2])*sizeof(fftw_complex));
	displacements = (fftw_complex*) fftw_malloc((size_t)(cells_displ[0]*cells_displ[1]*cells_displ[2])*sizeof(fftw_complex));
	ZA_clearMemory();
	ZA_displacements_plan = fftw_plan_dft_3d(cells_displ[0], cells_displ[1], cells_displ[2], displacements_fourier, displacements, FFTW_BACKWARD, FFTW_ESTIMATE);
	return 0;
}

int ZA_clearMemory() {
	for (int i = 0; i < cells_displ[0]*cells_displ[1]*cells_displ[2]; i++) {
		displacements_fourier[i][0] = 0.0;
		displacements_fourier[i][1] = 0.0;
		displacements[i][0] = 0.0;
		displacements[i][1] = 0.0;
	}
	return 0;
}

int ZA_freeMemory() {
	fftw_free(displacements_fourier);
	fftw_free(displacements);
	return 0;
}

int clearPk() {
	for (int i = 0; i < spectrum_size; i++) {	
		power_spectrum_monopole[i] = 0.0;
		power_spectrum_quadrupole[i] = 0.0;
		fk_abs_squares[i] = 0.0;		
		k_number[i] = 0;
		bin_k_average[i] = 0.0;
		binsums_L2[i][0] = 0.0;
		binsums_L2[i][1] = 0.0;
		binsums_L2[i][2] = 0.0;
		bin_k_min[i] = 1e10;
		bin_k_max[i] = 0.0;
		toRemove_bins[i] = false;
	}
	return 0;
}

//allocates memory given how many bins are wanted, returns step size dk
double init_spectrum_storage_len(int bins, double kmin, double kmax) {
	double dk = (kmax - kmin) / (double) (bins-1);
	fk_abs_squares = calloc((size_t)bins, sizeof(*fk_abs_squares));
	power_spectrum_monopole = calloc((size_t)bins, sizeof(*power_spectrum_monopole));
	power_spectrum_quadrupole = calloc((size_t)bins, sizeof(*power_spectrum_quadrupole));	
	k_number = calloc((size_t)bins, sizeof(*k_number));	
	bin_k_average = calloc((size_t)bins, sizeof(*bin_k_average));
	bin_k_min = malloc((size_t)bins * sizeof(*bin_k_min));
	//start with something BIG
	for (int i = 0; i < bins; i++) {
		bin_k_min[i] = 1e10;
	}
	bin_k_max = calloc((size_t)bins, sizeof(*bin_k_max));
	toRemove_bins = calloc((size_t)bins, sizeof(*toRemove_bins));
	calloc2D_double(&binsums_L2, bins, 3);		
	return dk;
}

int free_spectrum_storage() {
	free(fk_abs_squares);
	free(power_spectrum_monopole);
	free(power_spectrum_quadrupole);
	free(k_number);
	free(bin_k_average);
	free(bin_k_min);
	free(bin_k_max);
	free(toRemove_bins);
	free2D_double(&binsums_L2, 3);
	return 0;
}

void calloc2D_double(double*** arr, int dim1, int dim2) {
	*arr = malloc((size_t) dim1 * sizeof(**arr));
	for (int i = 0; i < dim1; i++) {
		(*arr)[i] = calloc((size_t)dim2, sizeof(***arr));
	}	
}

void free2D_double(double*** arr, int dim2) {
	for (int i = 0; i < dim2; i++) {
		free((*arr)[i]);
	}
	free(*arr);
	arr = NULL;
}

int clearMemory() {	
	clearPk();		
	clearGrid();
	clearParticles();
	return 0;
}

int memCheck() {
	bool pass = true;
	printf("particles check: \n");
	for (int i = 0; i < particle_no; i++) {
		if (particles[i][0] != 0.0 || particles[i][1] != 0.0 || particles[i][2] != 0.0 || particles[i][3] != 0.0 || particles[i][4] != 0.0) {
			printf("particles not empty\n");
			pass = false;
			break;
		}
	}
	if (particle_no > 0) {
		printf("particles: %d\n", particle_no);
	}	
	if (central_no > 0) {
		printf("centrals: %d\n", central_no);
	}
	if (satellite_no > 0) {
		printf("satellites: %d\n", satellite_no);
	}
	printf("grid check: \n");
	for (int i = 0; i < cells[0] * cells[1] * cells[2]; i++) {
		if (grid[i] != 0.0) {
			printf("grid not empty\n");
			pass = false;
			break;
		}	
	}
	printf("overdensity check: \n");	
	for (int i = 0; i < cells[0] * cells[1] * cells[2]; i++) {
		if (overdensity[i][0] != 0.0 || overdensity_fourier[i][0] != 0.0) {
			printf("overdensity arrays not empty\n");
			pass = false;
			break;
		}	
	}
	printf("ZA disp arrays check: \n");
	for (int i = 0; i < cells_displ[0]*cells_displ[1]*cells_displ[2]; i++) {
		if (displacements[i][0] != 0.0 || displacements_fourier[i][0] != 0.0) {
			printf("displacement arrays not empty\n");
			pass = false;
			break;
		}	
	}
	if (Jenkins_foldfactor != 1.0) {
		printf("fold factor not 1\n");
		pass = false;
	}
	printf("Checking Pk calc data: \n");
	for (int i = 0; i < spectrum_size; i++) {
		if (fk_abs_squares[i] != 0.0 || power_spectrum_monopole[i] != 0.0 || power_spectrum_quadrupole[i] != 0.0 || k_number[i] != 0 || toRemove_bins[i] != false || binsums_L2[i][0] != 0.0 || binsums_L2[i][1] != 0.0 || binsums_L2[i][2] != 0.0) {
			printf("Pk calculation data not empty\n");
			pass = false;
			break;
		} 
	}
	if (pass) {
		printf("memory empty\n");
		return 0;
	} else {
		return 0;
	}
}

