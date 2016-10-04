#define M_PI 3.14159265358979323846264338327

double randomDouble(double min, double max) {	
	return fmod((double) rand(),(max - min)) + min;
}





int populateParticles() {	
	
	for (int i = 0; i < particle_no; i++) {			
		particles[i][0] = randomDouble(0, limit_x);		
		particles[i][1] = randomDouble(0, limit_y);
		particles[i][2] = randomDouble(0, limit_z);
	}
		
	return 0;

}







int populateGrid() {

	int grid_x, grid_y, grid_z;	

	for (int i = 0; i < particle_no; i++) {
		grid_x = (int) floor((particles[i][0] / limit_x) * cells_x);
		grid_y = (int) floor((particles[i][1] / limit_y) * cells_y);
		grid_z = (int) floor((particles[i][2] / limit_z) * cells_z);		
		
		grid[grid_x + cells_x * grid_y + cells_x * cells_y * grid_z] += 1;		
	}

	return 0;

}






int gridIntoOverdensity() {

	//average density (particles per grid box) simply particle_no / size of box
	//box effectively has unit dimensions
	double rho_avg = (double) particle_no / ((double)(cells_x * cells_y * cells_z));
	for (int i = 0; i < cells_x; i++) {
		for (int j = 0; j < cells_y; j++) {
			for (int k = 0; k < cells_z; k++) {	
				//density of particles of a single box = number of particles in that box
				//overdensity = number of particles in box / rho_avg - 1			 
				grid[i + cells_x * j + cells_x * cells_y * k] = (grid[i + cells_x * j + cells_x * cells_y * 				k] / rho_avg) - 1; 
			}
		}
	}

	return 0;
}







int fourier_coefficient_grid(double kx, double ky, double kz, double *coef_real, double *coef_im) {

	double overdensity_value, kdotr, x, y, z, norm, length_over_cells;
	norm = (double) cells_x * cells_y * cells_z;
	length_over_cells = limit_x / cells_x;
	*coef_real = 0.;
	*coef_im = 0.;		
	//could be different, same for simplicity here
	
	
	for (int i = 0; i < cells_x; i++) {
		for (int j = 0; j < cells_y; j++) {
			for (int k = 0; k < cells_z; k++) {				
				overdensity_value = grid[i + cells_x * j + cells_x * cells_y * k];	
				x = i * length_over_cells;		
				y = j * length_over_cells;	
				z = k * length_over_cells;		
				kdotr = kx * x + ky * y + kz * z;				
				*coef_real += overdensity_value * cos(kdotr);
				*coef_im -= overdensity_value * sin(kdotr);
			}
		}
	}
	
	//normalisation	
	*coef_real /= norm;
	*coef_im /= norm;

	return 0;

}




int fourier_coefficient_particles(double kx, double ky, double kz, double *coef_real, double *coef_im) {
	
	double kdotr, x, y, z;
	*coef_real = 0.;
	*coef_im = 0.;	
	
	for (int i = 0; i < particle_no; i++) {
		x = particles[i][0];
		y = particles[i][1];
		z = particles[i][2];
		kdotr = kx * x + ky * y + kz * z;
		*coef_real += cos(kdotr);
		*coef_im -= sin(kdotr);	
	}
	return 0;
}




//loop from 0 or 1??
int overdensity_pk() {


    double kx, ky, kz, coef_real, coef_im, Pk, k, k_to_index;
    int index, loopmax;
    coef_im = 0.;
    coef_real = 0.;
    loopmax = 3;
    
    //100 comes from size of power_spectrum array, sqrt(3) because when kx = ky = kz = PI, k = PI*sqrt(3)
    k_to_index = 100. / (M_PI * sqrt(3.));
       
	
	//k from -pi to +pi here
	//f(k) = [f(-k)]*
	//just 0 to pi enough
	// k = 0 unnecessary??? only 0,0,0
	//numerical recipes 12.1.8 - saying it's continuous here, must multiply by delta (= 1 here)?
	for (int a = 0; a < loopmax; a++) {
		for (int b = 0; b < loopmax; b++) {
			for (int c = 0; c < loopmax; c++) {	//c++ :D			
				kx = a * M_PI * (double) cells_x / ((double) loopmax * limit_x);
				ky = b * M_PI / ((double) loopmax);
				kz = c * M_PI / ((double) loopmax);
				//printf("in loop, a = %d, b = %d, c = %d", a, b, c);
				
				//fourier_coefficient_grid(kx, ky, kz, &coef_real, &coef_im);
				fourier_coefficient_particles(kx, ky, kz, &coef_real, &coef_im);

				//two contributions f(k) and f(-k)
				//in Pk both result in f(k)f(k)* = f(-k)f(-k)*
				Pk = 2*(coef_real * coef_real + coef_im * coef_im);
				k = sqrt(kx*kx + ky*ky + kz*kz);
				
				index = (int) floor(k * k_to_index);
				power_spectrum[index] += Pk;
			} 
		} 
	} 

	return 0;
}


int pk_to_file() {

	FILE *f = fopen("pk.txt", "w");
	double k, index_to_k;
	char output_pk[10];
	char output_k[6];
	char output_final[19];	
	index_to_k = 100. / (M_PI * sqrt(3.));
	
	for (int i = 0; i < 100; i++) {
	
		//reset
		output_final[0] = '\0';
		output_pk[0] = '\0';
		output_k[0] = '\0';
	
		k = ((double) i) / index_to_k;		
	
		//turn into strings, safe print prevents errors from buffer overflow
		snprintf(output_pk, 10*sizeof(char), "%f", power_spectrum[i]);		
		snprintf(output_k, 10*sizeof(char), "%f", k);
		
		printf("%s", output_pk);
		
		//piece together output
		strcat(output_final, output_k);
		strcat(output_final, " ");
		strcat(output_final, output_pk);
		strcat(output_final, "\n");
		fprintf(f, "%s", output_final);

		
	
	}
	



	fclose(f);
	return 0;

}


