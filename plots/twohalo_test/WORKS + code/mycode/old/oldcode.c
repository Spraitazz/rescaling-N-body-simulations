//	CHANGE TO EITHER COORDINATES OF MIDDLE OF CELL, OR COORDINATES OF CENTRE OF MASS OF CELL
/*
int fourier_coefficient_grid_old(double kx, double ky, double kz, double *coef_real, double *coef_im) {
	double overdensity_value, kdotr, x, y, z, norm, cellSize;
	
	norm = (double) cells_x * cells_y * cells_z;
	
	cellSize = limit_x / cells_x;
	
	*coef_real = 0.;
	*coef_im = 0.;		
	
	for (int i = 0; i < cells_x; i++) {
		for (int j = 0; j < cells_y; j++) {
			for (int k = 0; k < cells_z; k++) {				
				overdensity_value = grid[i + cells_x * j + cells_x * cells_y * k];	
				
				x = i * cellSize;	// mid-point cell. or mean position of particles in that cell. 	
				y = j * cellSize;	
				z = k * cellSize;		
				
				kdotr = kx * x + ky * y + kz * z;				
				
				*coef_real += overdensity_value * cos(kdotr);
				*coef_im   -= overdensity_value * sin(kdotr);
			}
		}
	}
	
	//normalisation	
	*coef_real /= norm;
	*coef_im   /= norm;

	return 0;

}
*/


/*
int pk_linearbin(int index, double k, double kmin, double kmax, double dk) {

	double coef_real, coef_im, fk_squared;
	int index_pk;

	//Can change to check against ksquared values, a bit faster, take sqrt inside only
	if ((k >= kmin-0.000001) && (k <= kmax+0.000001)) { 
										
		coef_real = overdensity_fourier_out[index][0];
		coef_im = overdensity_fourier_out[index][1];				
		fk_squared = coef_real * coef_real + coef_im * coef_im;						
		index_pk = (int) floor((k-kmin) / dk);				
		power_spectrum[index_pk] += fk_squared;
		k_number[index_pk] += 1;			
		
	} 
	
	return 0;
}
*/



int fourier_coefficient_particles(double kx, double ky, double kz, double *coef_real, double *coef_im) {
	double kdotr, x, y, z, norm;
	
	//norm = (double) cells_x * cells_y * cells_z;
	
	*coef_real = 0.;
	*coef_im = 0.;	
	
	for (int i = 0; i < particle_no; i++) {
		x = particles[i][0];
		y = particles[i][1];
		z = particles[i][2];
		kdotr = kx * x + ky * y + kz * z;
		
		*coef_real += cos(kdotr);
		*coef_im   -= sin(kdotr);	
	}
	
	// normalisation	
	//*coef_real /= norm;
	//*coef_im   /= norm;
	
	return 0;
}



//"continuous" case for power spectrum calculation

/*
	volume_per_mode = fundmode_x * fundmode_y * fundmode_z;
    dV_prefactor = 4 * M_PI * dk;
	//dV = dV_prefactor * k_squared;
	//number of modes in given volume of k-space, each mode has fundmode_x * 						fundmode_y * fundmode_z volume
	//no_modes = dV / volume_per_mode;
	
	
	//now divide by number of modes at this k (mean value)
	//Pk = fk_squared / no_modes;	

*/


/*
int coefficients_to_file() {

	FILE *f = fopen("coef.txt", "w");	
	//double k;
	
	for (int i = 0; i < spectrum_size; i++) {		
		//k = kmin + dk / 2.  + (double) i * dk;
		fprintf(f, "%.5lf \t %d\n", power_spectrum[i], k_number[i]);	
	}
	
	fclose(f);
	return 0;
}
*/



/*
int pk_to_file_linearbin() {

	FILE *f = fopen("pk.txt", "w");
	double k; 
		
	for (int i = 0; i < spectrum_size; i++) {
		k = kmin + dk / 2.  + (double) i * dk;	
        fprintf(f, "%.5lf \t %.5lf \t %d\n", k, power_spectrum[i], k_number[i]);
	}

	fclose(f);
	return 0;
}
*/

//add a single power spectrum to the array of power spectra for all sims
int spectra_add(int sim_no) {

	for (int i = 0; i < spectrum_size_final; i++) {
		power_spectra_monopole[sim_no][i] = power_spectrum_monopole[i];
		power_spectra_quadrupole[sim_no][i] = power_spectrum_quadrupole[i];
	}

	return 0;
}

int fill_covariance() {


	
	for (int i = 0; i < spectrum_size_final; i++) {
		for (int j = 0; j < spectrum_size_final; j++) {	
		
			double sum_same[3] = {0.};
			double sum_cross[3] = {0.};	
					
			for (int sim = 0; sim < simulation_number - 1; sim++) {
			    sum_same[0] += power_spectra_monopole[sim][i] * power_spectra_monopole[sim][j];		
			    sum_same[1] += power_spectra_monopole[sim][i] * power_spectra_quadrupole[sim][j];
			    sum_same[2] += power_spectra_quadrupole[sim][i] * power_spectra_quadrupole[sim][j];
			    
				for (int sim2 = sim + 1; sim2 < simulation_number; sim2++) {
					sum_cross[0] += power_spectra_monopole[sim][i] * power_spectra_monopole[sim2][j];
					sum_cross[1] += power_spectra_monopole[sim][i] * power_spectra_quadrupole[sim2][j];
					sum_cross[2] += power_spectra_quadrupole[sim][i] * power_spectra_quadrupole[sim2][j];
				}			
			}
			//one term lost:
			sum_same[0] += power_spectra_monopole[simulation_number - 1][i] * power_spectra_monopole[simulation_number - 1][j];
			sum_same[1] += power_spectra_monopole[simulation_number - 1][i] * power_spectra_quadrupole[simulation_number - 1][j];
			sum_same[2] += power_spectra_quadrupole[simulation_number - 1][i] * power_spectra_quadrupole[simulation_number - 1][j];
			
			covariance_matrix[0][i][j] = (sum_same[0] * (5. / 36.));			
			covariance_matrix[0][i][j] -= (sum_cross[0] / 18.);
			
			covariance_matrix[1][i][j] = (sum_same[1] * (5. / 36.));			
			covariance_matrix[1][i][j] -= (sum_cross[1] / 18.);	
				
			covariance_matrix[2][i][j] = (sum_same[2] * (5. / 36.));			
			covariance_matrix[2][i][j] -= (sum_cross[2] / 18.);					
		}	
	}

	return 0;
}

int covariance_to_file() {

	FILE *cov_P0P0 = fopen("covariance_GR_P0P0", "w");
	FILE *cov_P0P2 = fopen("covariance_GR_P0P2", "w");
	FILE *cov_P2P2 = fopen("covariance_GR_P2P2", "w");
	
	for (int i = 0; i < spectrum_size; i++) {
		for (int j = 0; j < spectrum_size; j++) {
			fprintf(cov_P0P0, "%.5lf \t", covariance_matrix[0][i][j]);
			fprintf(cov_P0P2, "%.5lf \t", covariance_matrix[1][i][j]);
			fprintf(cov_P2P2, "%.5lf \t", covariance_matrix[2][i][j]);
		}	
		fprintf(cov_P0P0, "\n");
		fprintf(cov_P0P2, "\n");
		fprintf(cov_P2P2, "\n");
	}

	fclose(cov_P0P0);
	fclose(cov_P0P2);
	fclose(cov_P2P2);
	return(0);
}


int initOverdensityFT() {
	overdensities_ft = calloc(cells_x * cells_y * cells_z, sizeof(*overdensities_ft));
	for (int i = 0; i < cells_x; i++) {
		for (int j = 0; j < cells_y; j++) {
			for (int k = 0; k < cells_z; k++) {
				overdensities_ft[i + cells_x * j + cells_x * cells_y * k] = (double *) calloc(2, sizeof(**overdensities_ft));
			}	
		}
	}
	return 0;
}

