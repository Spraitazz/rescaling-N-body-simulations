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
	//toRedshift();	
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

int haloes_measure_Pk(char outFile[], double foldfactor, int grid_func) {

	bool folded = false;
	
	if (foldfactor < 1.0) {
		printf("what are you doing? foldfactor: %lf\n", foldfactor);
		exit(0);
	} else if (foldfactor > 1.0) {
		Jenkins_foldfactor = foldfactor;
		folded = true;
		Jenkins_fold_volume();		
	}	
	
	if (grid_func == CIC) {
		populateGridCIC();
	} else if (grid_func == NGP) {
		populateGrid();	
	} else {
		printf("grid func neither NGP nor CIC in haloes_measure_Pk() \n");
		exit(0);
	}

	gridIntoOverdensity();	
	overdensity_pk();
	pk_to_file_logbin(outFile);	
	
	if (folded)	Jenkins_restore();

	clearGrid();
	overdensity_clearMemory();	
	free_spectrum_storage();
	return 0;
}

/*
int recursive_fold_automatic(char outFile[], double kmin, double kmax, int foldPoints, int** foldings, double k_folds[]) {
	printf("folding from kmin: %lf to kmax: %lf \n", kmin, kmax);
	char tmp_unfolded[100], tmp_folded[100], tmp_combined[100], tmp[100];
	FoldingInformation info;
	double kmax_cur, foldfactor;
	bool max_reached;
	
	//CHECK FOR KMIN - IF TOO BIG, START STRAIGHT FROM FOLDED

	//if (foldings == NULL) {
	//	int** new_foldings;
	//	calloc2D_int(&new_foldings, 10, 2);
	
	if (foldings == NULL) {
	
		
	
	}
	
	sprintf(tmp_unfolded, "%s/tmp1.dat", temp_directory);
	sprintf(tmp_folded, "%s/tmp2.dat", temp_directory);
	sprintf(tmp_combined, "%s/tmp3.dat", temp_directory);
	haloes_measure_Pk(tmp_unfolded, 1.0, false);
	haloes_measure_Pk(tmp_folded, 2.0, false);
	info = combine_folded(tmp_unfolded, tmp_folded, pk_format, tmp_combined, kmin, kmax, -1.0, true);
	foldings[0][0] = info.unfold_last_ind;
	foldings[0][1] = info.fold_first_ind;
	max_reached = info.max_reached;	
	
	foldfactor = 2.0;
	while (!max_reached) {
		sprintf(tmp, "%s", tmp_unfolded);
		sprintf(tmp_unfolded, "%s", tmp_combined); //now the previously combined file is "unfolded"	
		sprintf(tmp_combined, "%s", tmp); //output
		
		foldfactor *= 2.0;
		
		haloes_measure_Pk(tmp_folded, foldfactor, false);
		info = combine_folded(tmp_unfolded, tmp_folded, pk_format, tmp_combined, kmin, kmax, -1.0, true);
		max_reached = info.max_reached;

	}
	
	rename(tmp_combined, outFile);	
	remove(tmp_unfolded);
	remove(tmp_folded);

	printf("done folding \n");
	return 0;
}
*/





