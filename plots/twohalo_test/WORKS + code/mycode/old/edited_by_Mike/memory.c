int initParticles() {
	particles = calloc(particle_no, sizeof(*particles));
	for (int i = 0; i < particle_no; i++) {
		particles[i] = (double *) calloc(3, sizeof(**particles));
	}
	return 0;
}




int initGrid() {
	grid = calloc(cells_x * cells_y * cells_z, sizeof(*grid));	
	return 0;
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


/*
int initPowerSpectrum() {
	power_spectrum = malloc(10
	return 0;
}

*/




/*


int freeGrid(int cells_x, int cells_y, int cells_z, double ***grid_tmp) {

	
	for (int x = 0; x < cells_x; x++) {		
		for (int y = 0; y < cells_y; y++) {
			free(grid_tmp[x][y]);
		}
		free(grid_tmp[x]);
	}
	
	free(grid_tmp);
	
	return 0;

}

*/



