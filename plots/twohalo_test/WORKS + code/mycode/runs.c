int basic_run(char outFile[]) {
	//set stuff, particle no etc.
	
	//overlay grid over particles
	populateGridCIC();
	
	//grid - > overdensity_in
	gridIntoOverdensity();
	
	//overdensity_in -> FFT -> overdensity_fourier_out -> bin by |k| -> / by # modes per bin	
	overdensity_pk(false);		

	pk_to_file_logbin(outFile);		
	//delsq_to_file_logbin(outFile);
	
	
	return 0;
}

int pre_covariance_run(const char inputData[], const char input_format[], const char outFile[], bool save_files, bool fold, int skipLines) {
	
	//reads given file, places all particle positions in particles array	
	read_full_xyz(inputData, skipLines, input_format);	
	toRedshift();	
	if (fold) {
		Jenkins_foldfactor = 4.0;
		Jenkins_fold_volume();		
	}

	populateGridCIC();	
	gridIntoOverdensity();	
	overdensity_pk(false);		
	if (save_files) pk_to_file_logbin(outFile);		
	if (fold) Jenkins_restore();	
	
	clearMemory();	
	return 0;	
}

int haloes_measure_Pk(char outFile[], double foldfactor, bool CIC) {

	bool folded = false;
	
	if (foldfactor < 1.0) {
		printf("what are you doing? foldfactor: %lf\n", foldfactor);
		exit(0);
	} else if (foldfactor > 1.0) {
		Jenkins_foldfactor = foldfactor;
		folded = true;
		Jenkins_fold_volume();		
	}	
	
	if (CIC) {
		populateGridCIC();
	} else {
		populateGrid();	
	}

	gridIntoOverdensity();	
	overdensity_pk(false);
	pk_to_file_logbin(outFile);	
	
	if (folded)	Jenkins_restore();

	clearGrid();
	overdensity_clearMemory();	
	free_spectrum_storage();
	return 0;
}

//I CALL GAUSSIAN RANDOMS THOSE THAT ARE GIVEN A VALUE FROM A DISTRIBUTION OR OTHERWISE ASSIGNED
int randoms_fullRun(int particle_number, const char outFile[], bool fold, bool test_mode, bool gaussian) {
	particle_no = particle_number;
	
	if (gaussian == false) {	
	initParticles();
	populateParticles();
	}
	//printf("random particle [10][0]: %f\n", particles[10][0]);
	if (fold) {
		Jenkins_foldfactor = 2.0;
		Jenkins_fold_volume();		
	}
	
	populateGrid();
	gridIntoOverdensity();
	overdensity_pk(gaussian);
	
	if (test_mode) {
		//send to memory, pick up at tests.c
	
	} else {
		pk_to_file_logbin(outFile);
	}

	//clean up
	clearMemory();	
	volume_limits[0] = 1500.0;  // h^-1 Mpc 
	volume_limits[1] = 1500.0;
	volume_limits[2] = 1500.0;	
	
	return 0;
}

int noise_randoms(int particle_number, const char outFile[], bool fold, bool test_mode) {
	randoms_fullRun(particle_number, outFile, fold, test_mode, false);
	return 0;
}

int Gaussian_randoms(int particle_number, const char outFile[], bool fold) {
	randoms_fullRun(particle_number, outFile, fold, false, true);
	return 0;
}





