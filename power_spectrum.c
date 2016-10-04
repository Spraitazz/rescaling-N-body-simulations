int populateParticles() {	
	for (int i = 0; i < particle_no; i++) {	
		//positions		
		particles[i][0] = randomDouble(0.0, volume_limits[0]);	
		particles[i][1] = randomDouble(0.0, volume_limits[1]);
		particles[i][2] = randomDouble(0.0, volume_limits[2]);
		
		// velocities
		particles[i][3] = 0.0;	
		particles[i][4] = 0.0;
		particles[i][5] = 0.0;	
		
		// mass
		particles[i][6] = 0.0;
	}				
	return 0;
}

int populateHaloes(double halo_mass) {	
	for (int i = 0; i < halo_no; i++) {	
		//positions		
		particles[i][0] = randomDouble(0.0, volume_limits[0]);	
		particles[i][1] = randomDouble(0.0, volume_limits[1]);
		particles[i][2] = randomDouble(0.0, volume_limits[2]);
		
		// velocities
		particles[i][3] = 0.0;	
		particles[i][4] = 0.0;
		particles[i][5] = 0.0;	
		
		// mass
		particles[i][6] = halo_mass;
	}				
	return 0;
}

int populateParticles_lattice() {
	if (particle_no != cells[0]*cells[1]*cells[2]) {
		printf("for lattice particle number should equal cell number\n");
		exit(0);
	}	
	double cellSizes[3];
	cellSizes[0] = volume_limits[0] / (double) cells[0];
	cellSizes[1] = volume_limits[1] / (double) cells[1];
	cellSizes[2] = volume_limits[2] / (double) cells[2];	
	
	//placing particle in center of cell
	for (int i = 0; i < cells[0]; i++) {
		for (int j = 0; j < cells[1]; j++) {
			for (int k = 0; k < cells[2]; k++) {
				particles[arr_ind(i,j,k)][0] = ((double)i + 0.5) * cellSizes[0];
				particles[arr_ind(i,j,k)][1] = ((double)j + 0.5) * cellSizes[1];
				particles[arr_ind(i,j,k)][2] = ((double)k + 0.5) * cellSizes[2];
								
				// velocities
				particles[arr_ind(i,j,k)][3] = 0.0;	
				particles[arr_ind(i,j,k)][4] = 0.0;
				particles[arr_ind(i,j,k)][5] = 0.0;	
		
				// mass
				particles[arr_ind(i,j,k)][6] = 0.0;
			}
		}
	}

	return 0;
}

int toRedshift(double a) {		
	for (int i = 0; i < particle_no; i++) {		
		particles[i][2] += particles[i][5] / (a*H0);	//z-axis line of sight		
		PBC(&particles[i][0], &particles[i][1], &particles[i][2]);
	}	
	return 0;
}

int populateGrid() {
	int grid_x, grid_y, grid_z, index;		
	for (int i = 0; i < particle_no; i++) {
		grid_x = (int) floor((particles[i][0] / volume_limits[0]) * (double) cells[0]);
		grid_y = (int) floor((particles[i][1] / volume_limits[1]) * (double) cells[1]);
		grid_z = (int) floor((particles[i][2] / volume_limits[2]) * (double) cells[2]);		
		index = arr_ind(grid_x, grid_y, grid_z);
		grid[index] += 1.0;		
	}
	WindowFunc_pow = 1.0;
	return 0;
}

//http://background.uchicago.edu/~whu/Courses/Ast321_11/pm.pdf p14
int populateGridCIC() {
	int grid_x, grid_y, grid_z;
	double xc, yc, zc, dx, dy, dz, tx, ty, tz;
	double cellSizes[3];
	cellSizes[0] = volume_limits[0] / (double) cells[0];
	cellSizes[1] = volume_limits[1] / (double) cells[1];
	cellSizes[2] = volume_limits[2] / (double) cells[2];	
	
	for (int i = 0; i < particle_no; i++) {
		grid_x = (int) floor((particles[i][0] / volume_limits[0]) * (double) cells[0]);
		grid_y = (int) floor((particles[i][1] / volume_limits[1]) * (double) cells[1]);
		grid_z = (int) floor((particles[i][2] / volume_limits[2]) * (double) cells[2]);			
		
		xc = cellSizes[0] * ((double) grid_x + 0.5);
		yc = cellSizes[1] * ((double) grid_y + 0.5);
		zc = cellSizes[2] * ((double) grid_z + 0.5);
		
		//must be from 0 to 1
		dx = fabs(particles[i][0] - xc) / cellSizes[0];
		dy = fabs(particles[i][1] - yc) / cellSizes[1];
		dz = fabs(particles[i][2] - zc) / cellSizes[2];	
		
		tx = 1.0 - dx;
		ty = 1.0 - dy;
		tz = 1.0 - dz;
		
		//ensure PBC
		int xplus, yplus, zplus;
		xplus = (grid_x + 1) % cells[0];
		yplus = (grid_y + 1) % cells[1];
		zplus = (grid_z + 1) % cells[2];		

		grid[arr_ind(grid_x, grid_y, grid_z)] += tx*ty*tz;	
		grid[arr_ind(xplus, grid_y, grid_z)] += dx*ty*tz;
		grid[arr_ind(xplus, yplus, grid_z)] += dx*dy*tz;	
		grid[arr_ind(xplus, grid_y, zplus)] += dx*ty*dz;
		grid[arr_ind(xplus, yplus, zplus)] += dx*dy*dz;			
		grid[arr_ind(grid_x, yplus, grid_z)] += tx*dy*tz;		
		grid[arr_ind(grid_x, grid_y, zplus)] += tx*ty*dz;
		grid[arr_ind(grid_x, yplus, zplus)] += tx*dy*dz;		
	}
	WindowFunc_pow = 2.0;
	return 0;
}

int gridIntoOverdensity() {
	//average density (particles per grid box) simply particle_no / size of box
	double N_avg = (double) particle_no / ((double) cells[0]*cells[1]*cells[2]);	
	int index;	
	for (int i = 0; i < cells[0]; i++) {
		for (int j = 0; j < cells[1]; j++) {
			for (int k = 0; k < cells[2]; k++) {	
				//density of particles of a single box = number of particles in that box
				//overdensity = number of particles in box / rho_avg - 1	
				index = arr_ind(i, j, k);				
				overdensity[index][0] = (grid[index] / N_avg) - 1.0;			
			}
		}
	}
	return 0;
}

//gaussian here just means use overdensity fourier coeffs from a gaussian random field (NOT IMPLEMENTED)
int overdensity_pk() {
    double kx, ky, kz, k, logdk, logkmin, logkmax, WindowFunc, mu, kmin, kmax;
    double fundmodes[3], nyquist[3], cellSizes[3];
    bool flag;
    int index;  
    
    cellSizes[0] = volume_limits[0] / (double) cells[0];
	cellSizes[1] = volume_limits[1] / (double) cells[1];
	cellSizes[2] = volume_limits[2] / (double) cells[2];	

    fundmodes[0] = 2.0*pi/volume_limits[0];
    fundmodes[1] = 2.0*pi/volume_limits[1];
    fundmodes[2] = 2.0*pi/volume_limits[2];
    
    nyquist[0] = pi / cellSizes[0];
    nyquist[1] = pi / cellSizes[1];
    nyquist[2] = pi / cellSizes[2];
       
    //should be min(fundmodes), different law for freqs 
    kmin = fundmodes[0] * sqrt(3.0);
    kmax = nyquist[0] * sqrt(3.0);
    
    //get fourier overdensities
	fftw_execute(overdensity_plan);		

    logkmax = log10(kmax);
	logkmin = log10(kmin);
    logdk = init_spectrum_storage_len(spectrum_size, logkmin, logkmax); 

    for (int aa = 0; aa < cells[0]; aa++) {
    	for (int bb = 0; bb < cells[1]; bb++) {
    		for (int cc = 0; cc < cells[2]; cc++) {
    			kx = (double)aa * fundmodes[0];
    			ky = (double)bb * fundmodes[1];
				kz = (double)cc * fundmodes[2];  			
				
				//convention - up to N/2 positive frequencies, N/2 is fc, then -fc+ffund, ...
				if (aa > cells[0]/2) kx -= ((double) cells[0]) * fundmodes[0]; //or -= 2*kx_nyquist
				if (bb > cells[1]/2) ky -= ((double) cells[1]) * fundmodes[1];
				if (cc > cells[2]/2) kz -= ((double) cells[2]) * fundmodes[2];
				
				index = arr_ind(aa, bb, cc);					
				
				overdensity_fourier[index][0] /= ((double)cells[0]*cells[1]*cells[2]);
				overdensity_fourier[index][1] /= ((double)cells[0]*cells[1]*cells[2]);
				
				//Window function (convolution in real -> multip in fourier) to reduce aliased power
				WindowFunc = 1.0;
				if(kx != 0.0) WindowFunc *= sin(pi*kx*0.5/nyquist[0])/(pi*kx*0.5/nyquist[0]);
				if(ky != 0.0) WindowFunc *= sin(pi*ky*0.5/nyquist[1])/(pi*ky*0.5/nyquist[1]);
				if(kz != 0.0) WindowFunc *= sin(pi*kz*0.5/nyquist[2])/(pi*kz*0.5/nyquist[2]);
									
				// Correct for mass assignment of randoms.
				// Cloud in cell pow 2.0. NGP corresponds to WindowFunc rather than WindowFunc^2. 
				overdensity_fourier[index][0] /= pow(WindowFunc, WindowFunc_pow);
				overdensity_fourier[index][1] /= pow(WindowFunc, WindowFunc_pow);
				
				k = sqrt(kx*kx + ky*ky + kz*kz);				
				//unit vector along k dot n, chosen along polar axis kz
				mu = kz / k;
				
				// independence_flag? Stars_and_Stripes.  one half independent. 
				flag = false;				
				if (kz > 0.0) {	flag = true; }
				else if ((kz == 0.0) && (ky > 0.0)) { flag = true; } 
				else if ((kz == 0.0) && (ky == 0.0) && (kx > 0.0)) { flag = true; } 
				
				if (flag && k >= kmin && k < kmax) pk_logbin(index, k, logkmin, logdk, mu);				
											  		
    		}    	
    	}    
	}
								        
	//find the average k for each bin - where the P(k) values are plotted, check for empty bins
	for (int i = 0; i < spectrum_size; i++) {
		if (k_number[i] != 0) {
			bin_k_average[i] /= (double) k_number[i];
		} else {
			toRemove_bins[i] = true;
		}
	}		
	//finally, calculate multipoles (mono, quadro)
	multipole_calc();	
	return 0;
}

int pk_logbin(int index, double k, double logkmin, double logdk, double mu) {
	double fk_squared, coef_real, coef_im, L2;
	int index_log_pk;	

	coef_real = overdensity_fourier[index][0];
	coef_im = overdensity_fourier[index][1];		
	fk_squared = coef_real * coef_real + coef_im * coef_im;	
	
	//plot is VP(k) versus k
	fk_squared *= volume;
	
	//SHOT NOISE, USE CAREFULLY
	//fk_squared -= volume/particle_no;	
		
	L2 = 0.5 * (3.0 * mu * mu - 1.0);							
	index_log_pk = (int) floor((log10(k) - logkmin)/logdk);	
	fk_abs_squares[index_log_pk] += fk_squared;
	
	k_number[index_log_pk] += 1;	
	bin_k_average[index_log_pk] += k;
	
	if (k > bin_k_max[index_log_pk]) {
		bin_k_max[index_log_pk] = k;
	}
	if (k < bin_k_min[index_log_pk]) {
		bin_k_min[index_log_pk] = k;
	}	
	
	binsums_L2[index_log_pk][0] += L2;
	binsums_L2[index_log_pk][1] += L2 * L2;	
	binsums_L2[index_log_pk][2] += L2 * fk_squared;	
	return 0;
}

//monopole + qudrupole
int multipole_calc() {	
	double margin, determinant;	
	margin = 1e-6;
	for (int i = 0; i < spectrum_size; i++) {		
		determinant = ((double) k_number[i]) * binsums_L2[i][1] - binsums_L2[i][0] * binsums_L2[i][0];
	
		//only has a meaningful value for a nonzero determinant
		if (determinant > margin) {			
			power_spectrum_monopole[i] = binsums_L2[i][1]*fk_abs_squares[i] - binsums_L2[i][0]*binsums_L2[i][2];
			power_spectrum_quadrupole[i] = ((double) k_number[i])*binsums_L2[i][2] - binsums_L2[i][0]*fk_abs_squares[i];
			power_spectrum_monopole[i] /= determinant;
			power_spectrum_quadrupole[i] /= determinant;
		} else {
			toRemove_bins[i] = true;
		}
	}		
	return 0;
}


//prints k, monopole (simple calculation), measurements in bin
int just_print_monopole(char filename[]) {
	FILE *pk_out_file = fopen(filename, "w");	
	for (int i = 0; i < spectrum_size; i++) {
		if (k_number[i] > 0) {
			fprintf(pk_out_file, "%lf \t %le \t %d \t %lf \t %lf \n", bin_k_average[i], fk_abs_squares[i]/(double)k_number[i], k_number[i], bin_k_min[i], bin_k_max[i]);
		}
	}
	fclose(pk_out_file);
	return 0;
}


//prints k, monopole, quadrupole, measurements in bin
int pk_to_file_logbin(const char filename[]) {
	FILE *pk_out_file = fopen(filename, "w");
	for (int i = 0; i < spectrum_size; i++) {
		if (toRemove_bins[i] == false) {			
			fprintf(pk_out_file, pk_format, bin_k_average[i], power_spectrum_monopole[i], power_spectrum_quadrupole[i], k_number[i], bin_k_min[i], bin_k_max[i]);
		} 
	}	
	fclose(pk_out_file);
	return 0;
}

int delsq_to_file_logbin(const char filename[]) {
	FILE *delsq_out_file = fopen(filename, "w");	
	double delsq;
	for (int i = 0; i < spectrum_size_final; i++) {				
		if (toRemove_bins[i] == false) {
			delsq = pow(bin_k_average[i], 3.0)*power_spectrum_monopole[i]/(2.0*pi*pi);
			fprintf(delsq_out_file, pk_format, bin_k_average[i], delsq, power_spectrum_quadrupole[i], k_number[i], bin_k_min[i], bin_k_max[i]);
		}	
	}
	fclose(delsq_out_file);
	return 0;
}


