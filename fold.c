int Jenkins_fold_volume() {	

	calloc2D_double(&particles_save, particle_no, 7); //1st and 2nd dimensions
	volume_limits[0] /= Jenkins_foldfactor;
	volume_limits[1] /= Jenkins_foldfactor;
	volume_limits[2] /= Jenkins_foldfactor;		
	
	for (int i = 0; i < particle_no; i++) {
		memcpy(particles_save[i], particles[i], 7 * sizeof(*particles_save[i]));
	
		particles[i][0] = fmod(particles[i][0], volume_limits[0]);
		particles[i][1] = fmod(particles[i][1], volume_limits[1]);
		particles[i][2] = fmod(particles[i][2], volume_limits[2]);	
	}
	return 0;
}

int Jenkins_restore() {
	for (int i = 0; i < particle_no; i++) {		
		memcpy(particles[i], particles_save[i], 7 * sizeof(*particles[i]));
	}
	
	volume_limits[0] *= Jenkins_foldfactor;
	volume_limits[1] *= Jenkins_foldfactor;
	volume_limits[2] *= Jenkins_foldfactor;	
	
	Jenkins_foldfactor = 1.0;	
	free2D_double(&particles_save, particle_no); //only second dimension needed for inner frees	
	return 0;
}

FoldingInformation combine_folded(const char unfolded_path[], const char folded_path[], char pk_format[], char outFile[], double kmin, double kmax, double k_fold, bool print) {

	FILE *unfolded_file = fopen(unfolded_path, "r");
	FILE *folded_file = fopen(folded_path, "r");
	FILE *output = fopen(outFile, "w");	
	
	double min_diff, k_diff, k_unfolded, k_folded, mono_unfolded;
	double mono_folded, perc_diff_k, perc_diff_mono;
	int len_folded, len_unfolded, min_k_ind_unfolded, min_k_ind_folded;
	int min_diff_ind, aliasing_index_unfolded, unfold_last_ind, fold_conn_ind;	
	int aliasing_index_folded;
	bool max_reached;

	//CHECK FOR NULL FILES
	if (unfolded_file == NULL) {
		printf("error opening %s\n", unfolded_path); 
		exit(0);
	}
	if (folded_file == NULL) {
		printf("error opening %s\n", folded_path);
		exit(0);
	}
	if (output == NULL) {
		printf("error opening %s\n", outFile);
		exit(0);
	}
	if (kmin < 0.0) {	
		printf("kmin can only be more or equal to 0.0\n");
		exit(0);	
	}
	
	//count lines
	len_folded = countLines(folded_file);
	len_unfolded = countLines(unfolded_file);

	//[k][mono][quadro]
	double** folded_values;
	double** unfolded_values;
	calloc2D_double(&folded_values, len_folded, 6);
	calloc2D_double(&unfolded_values, len_unfolded, 6);	
	
	//read into arrays from the respective files
	for (int i = 0; i < len_folded; i++) {
		ph = fscanf(folded_file, pk_format, &folded_values[i][0], &folded_values[i][1], &folded_values[i][2], &folded_values[i][3], &folded_values[i][4], &folded_values[i][5]);		
	}	
	
	for (int i = 0; i < len_unfolded; i++) {
		ph = fscanf(unfolded_file, pk_format, &unfolded_values[i][0], &unfolded_values[i][1], &unfolded_values[i][2], &unfolded_values[i][3], &unfolded_values[i][4], &unfolded_values[i][5]);
	}

	if (kmax > folded_values[len_folded - 1][0]) {
		printf("kmax: %lf, more than maximum in folded file \n", kmax);	
		exit(0);	
	}
	
	//-----------------------------------------------------------------------------------------
	
	//IMPORTANT: indices for connecting unfolded to folded
	
	//find the k index from which to begin in the unfolded file	
	int v = 0;
	while (!(unfolded_values[v][0] >= kmin)) {
		v+=1;
	}
	min_k_ind_unfolded = v;	
	
	//folded file start of aliasing
	v = len_folded - 1;	
	while (folded_values[v-1][1] < folded_values[v][1]) {
		v -= 1;
	}
	aliasing_index_folded = v+1;	
	
	//maximum is input, check for finish of file before aliasing kicks in
	if (kmax > 0.0) {
		v = 0;
		while (folded_values[v][0] < kmax) {
			v += 1;
		}
		
		aliasing_index_folded = v;
		
		//finish before aliasing
		if (v < aliasing_index_folded) {
			aliasing_index_folded = v;
			max_reached = true;
		}
	}
	
	//unfolded file start of aliasing
	v = len_unfolded - 1;	
	while (unfolded_values[v-1][1] < unfolded_values[v][1]) {
		v -= 1;
	}
	aliasing_index_unfolded = v+1;
	

	
	//--------------------------------------------------------
	
	if (k_fold < 0.0) {
		//finding fold point automatically		
	
		//find the k index from which to begin in the folded file
		//SEND JENKINS FOLD FACTOR TO THIS FUNCTION
		v = 0;
		//double foldfactor = folded_values[0][0] / unfolded_values[0][0];
		//printf("foldfactor: %lf \n", foldfactor);
		while(!(folded_values[v][0] >= kmin*2.0)) {
			v+=1;
		}
		min_k_ind_folded = v;
	
		//find the minimum k and monopole differences combining pairs from the folded/unfolded files	
		min_diff = 1e10;

		//looping through unfolded file
		for (int i = min_k_ind_unfolded; i < aliasing_index_unfolded-1; i++) {
			k_unfolded = unfolded_values[i][0];
			mono_unfolded = unfolded_values[i][1];
			
			//2x folded file
			for (int j = min_k_ind_folded; j < aliasing_index_folded-1; j++) {
				k_folded = folded_values[j][0];
				mono_folded = folded_values[j][1];
				
				//make sure k is going up from unfolded to folded
				if (k_folded > k_unfolded) {
					perc_diff_k = (k_folded - k_unfolded) / k_unfolded;
					perc_diff_mono = fabs((mono_folded - mono_unfolded) / mono_unfolded);
					
					//could POTENTIALLY just make sure monopole decreases with increasing k here (now more general)
					if (perc_diff_mono*perc_diff_k < min_diff && perc_diff_k < 0.5) {
											
						min_diff = perc_diff_mono*perc_diff_k;
						//best connection indices
						unfold_last_ind = i;
						fold_conn_ind = j;						
						
					} 
				}
			}				
		}	
	} else {
		
		//folding point set, just find the necessary indices for connecting
		if (k_fold < kmin || k_fold > kmax) {
			printf("cannot fold below kmin OR above kmax\n");
			exit(0);
		} else {
			//just find the last index where k < k_fold in the unfolded file,
			//then jump to where k > k_fold in the folded file
			v = 0;
			while (unfolded_values[v][0] < k_fold) {
				v+=1;
			}
			unfold_last_ind = v;
			
			v = 0;
			while(folded_values[v][0] < k_fold) {
				v+=1;
			}
			fold_conn_ind = v;		
		}	
	}	
	
	int lines_out = unfold_last_ind - min_k_ind_unfolded + aliasing_index_folded - fold_conn_ind;
	
	//INFO TO SEND TO MULTI-FILE FUNCTION
	FoldingInformation toReturn = {
	.folded_values = folded_values,
	.unfolded_values = unfolded_values,
	.kmin = unfolded_values[min_k_ind_unfolded][0],
	.kmax = folded_values[aliasing_index_folded-1][0],
	.unfold_last_ind = unfold_last_ind,
	.fold_first_ind = fold_conn_ind,
	.min_k_ind = min_k_ind_unfolded,
	.lines = lines_out,
	.aliasing_index = aliasing_index_folded,
	.max_reached = max_reached
	};
	
	//-----------------------------------------------------------------------	

	//if printing a single example
	if (print) {
		
		double final_values[lines_out][6];
	
		for (int i = min_k_ind_unfolded; i < unfold_last_ind; i++) {

			for (int k = 0; k < 6; k++) {
				final_values[i - min_k_ind_unfolded][k] = unfolded_values[i][k];
			}
		}		

		for (int i = fold_conn_ind; i < aliasing_index_folded; i++) {					
			for (int k = 0; k < 6; k++) {
				final_values[(i - min_k_ind_unfolded) + unfold_last_ind - fold_conn_ind][k] = folded_values[i][k];	
			}
		}

		for (int i = 0; i < lines_out; i++) {
			if (final_values[i][1] > 1e-6) {
				fprintf(output, pk_format, final_values[i][0], final_values[i][1], final_values[i][2], (int)final_values[i][3], final_values[i][4], final_values[i][5]);
			}
		}	
	}	
	
	fclose(unfolded_file);
	fclose(folded_file);
	fclose(output);
	free2D_double(&folded_values, len_folded);
	free2D_double(&unfolded_values, len_unfolded);
	return toReturn;
}



//THIS PROBABLY DOESNT WORK
/* 
int combine_multiple(char paths_unfolded[], char paths_folded[], char outFiles[], double kmin, double kmax, double k_fold, int simNo, int first_index) {

	char pu[100];
	char pf[100];
	char out[100];
	FoldingInformation allInfo[simNo];
	FoldingInformation thisInfo;
			
	//collect info for all sims - where the transition from unfolded to folded happens
	for (int i = first_index; i < simNo + first_index; i++) {

		sprintf(pu, paths_unfolded, i);
		sprintf(pf, paths_folded, i);
		sprintf(out, outFiles, i);		
		//get indices, dont print
		//printf("%s\n", pu);
		thisInfo = combine_folded(pu, pf, out, kmin, kmax, k_fold, false);
		allInfo[i - first_index] = thisInfo;	
	}
	
	//sum how many sims jump from [ind] index in the unfolded file
	//[bin][0 - how many counts of this starting index, 1 - which index it connects to, 2 - kmin ind, 3 - foldover]
	int index_popularity[100][4];
	int ind;

	for (int i = 0; i < simNo; i++) {
		thisInfo = allInfo[i];
		ind = thisInfo.unfold_last_ind;
		index_popularity[ind][0] += 1;
		index_popularity[ind][1] = thisInfo.fold_first_ind;
		index_popularity[ind][2] = thisInfo.min_k_ind;		
		index_popularity[ind][3] = thisInfo.foldover_index;	
	}
	
	//find which [ind] index has greatest number of unfolded files ending there
	int ind_end_unfold, ind_start_fold, maxCount, maxAt, min_k_ind, foldover_index;
	maxCount = 0;
	for (int i = 0; i < 100; i++) {
		if (index_popularity[i][0] > maxCount) {
			maxCount = index_popularity[i][0];
			maxAt = i;
		}	
	}

	
	//set final variables
	ind_end_unfold = maxAt;
	ind_start_fold = index_popularity[maxAt][1];
	min_k_ind = index_popularity[maxAt][2];
	foldover_index = index_popularity[maxAt][3];
	
	FILE *final_output;
	//print the necessary files	
	for (int sim = 0; sim < simNo; sim++) {	
		thisInfo = allInfo[sim];
		sprintf(out, outFiles, sim + first_index);		
		final_output = fopen(out, "w");	
		double final_values[200][4];
		for (int i = min_k_ind; i < ind_end_unfold; i++) {
			for (int k = 0; k < 4; k++) final_values[i - min_k_ind][k] = thisInfo.unfolded_values[i][k];
		}	
		
		//printf("ind start fold: %d, foldover: %d, ind_end_unfold: %d\n", ind_start_fold, foldover_index, ind_end_unfold);		

		for (int i = ind_start_fold; i < foldover_index; i++) {
			for (int k = 0; k < 4; k++) {			
				final_values[ind_end_unfold - min_k_ind + i - ind_start_fold][k] = thisInfo.folded_values[i][k];
			}	
		}	
		
		for (int i = 0; i < 200; i++) {
			if (final_values[i][1] > 0.000001) {
				fprintf(final_output, "%.5lf \t %.5lf \t %.5lf \t %d\n", final_values[i][0], final_values[i][1], final_values[i][2], (int) final_values[i][3]);
			}
		}	
		fclose(final_output);	
	}
	return 0;
}

*/
