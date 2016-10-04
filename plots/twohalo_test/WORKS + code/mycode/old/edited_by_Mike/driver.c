#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_rng.h>

#include <fftw3.h>

#include "driver.h"
#include "memory.c"
#include "functions.c"


int main(int argc, char **argv) {	
	//seed rand()
        srand((unsigned)time(NULL)); // replace by gsl.
	
	//define volume for particles
	limit_x = 100.;
	limit_y = 100.;
	limit_z = 100.;
	
	//specify particle number, particle array
	particle_no = 10000;
	//particles = malloc(particle_no * 3 * sizeof(double));
	
	//define grid	
	cells_x = 100;
	cells_y = 100;
	cells_z = 100;

	// Random variable generation.
	gsl_rng_env_setup();

	gsl_ran_T = gsl_rng_default;

	gsl_ran_r = gsl_rng_alloc(gsl_ran_T);

	printf("\n\nGSL random number test: %le \n\n", gsl_rng_uniform(gsl_ran_r));

	// FFTw test	
	int j;

	fftw_plan     p, iplan;
	fftw_complex *in, *out;

	in                 = (fftw_complex*) fftw_malloc(10*10*10*sizeof(fftw_complex));
	out                = (fftw_complex*) fftw_malloc(10*10*10*sizeof(fftw_complex));

	for(j=0; j<10*10*10; j++) in[j][0] = gsl_rng_uniform(gsl_ran_r);
	for(j=0; j<10*10*10; j++) in[j][1] = 0.0;
	
	p                  = fftw_plan_dft_3d(10, 10, 10, in, out, FFTW_FORWARD,  FFTW_ESTIMATE);
	
	fftw_execute(p);

	for(j=0; j<20; j++) printf("%le \t %le \n", out[j][0], out[j][1]);
	
	//  double* grid;
	/*
	//create particles
	initParticles();
	populateParticles();
	
	//set grid memory, populate
	initGrid();
	populateGrid();
	
	//turn grid into array of overdensities
	gridIntoOverdensity();
	
	//initialize overdensity_ft array 
	//initOverdensityFT();
	
	//printf("\nstarted power spectrum calc \n");
	overdensity_pk();
	//printf("done, power_spectrum[10]: %f\n", power_spectrum[10]);
	pk_to_file();
	
	//printf("file ready\n\n");
	for (int i = 0; i < 100; i++) {
	  // printf("ps[%d]: %f\n", i, power_spectrum[i]);
	
	}
	*/
	
	return 0;

}
