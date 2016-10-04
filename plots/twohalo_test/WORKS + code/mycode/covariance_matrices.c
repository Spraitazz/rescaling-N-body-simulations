/*
Takes the input file path string, which should contain %d for iteration,
the output path string, which should also contain %d, the number of simulations
and the index (usually 0) of the first simulation.
INPUT FILE MUST HAVE THE FORMAT "%lf \t %lf \t %lf \t %d \n" indicating k, monopole,
quadrupole and bin number, respectively.
Generates covariance matrix and outputs it to the outPath, where the k-values
increase down and to the right in the file.
*/
int generate_covariance(char inPath[], char outPath[], int sim_number, int startIndex) {

	//check indices
	if (sim_number < 1 || startIndex < 0) {
		printf("wrong indices\n");
		return 1;
	}	
	
	FILE *inFile, *outFile;
	char in_filename[100];
	char out_filename[100];
	int lineno, binCount;		
	double sum_mono, sum_quadro;
	//Format: multipoles[sim no][bin no][mono/quadro]
	double multipoles[sim_number][100][2];	
	double temp[3] = {0.0, 0.0, 0.0};
	
	//collect multipoles to array
	for (int i = startIndex; i < sim_number + startIndex; i++) {
		sprintf(in_filename, inPath, i);
		inFile = fopen(in_filename, "r");
		//file exists check
		if (inFile == NULL) {
			printf("invalid covariance matrix input file %s\n", in_filename);
			return 1;		
		}
		lineno = countLines(inFile);
		rewind(inFile);
		//all files should have the same number of lines (can also check k-vals are equal)
		if (i > startIndex) {
			if (lineno != binCount) {
				printf("Input files have different lengths, file: %s\n", in_filename);
			}
		}	
		binCount = lineno;	
		
		for (int k = 0; k < lineno; k++) {
			ph = fscanf(inFile, "%*f \t %lf \t %lf \t %*d", &multipoles[i - startIndex][k][0], 
			&multipoles[i - startIndex][k][1]);		
		}		
		fclose(inFile);
	}	
	
	//variables that depended on binCount
	double averages_multipoles[binCount][2];
	//Matr(i,j) Format: [mono-mono, mono-quadro, quadro-quadro][i][j]
	double cov_matrix[3][binCount][binCount];
	double cov_matrix_final[3][binCount][binCount];
	
	//averages of multipoles and quadrupoles for all sims
	for (int bin = 0; bin < binCount; bin++) {
		sum_mono = 0.0;
		sum_quadro = 0.0;
		for (int sim = 0; sim < sim_number; sim++) {
			sum_mono += multipoles[sim][bin][0];
			sum_quadro += multipoles[sim][bin][1];
		}
		averages_multipoles[bin][0] = sum_mono / (double)sim_number;
		averages_multipoles[bin][1] = sum_quadro / (double)sim_number;
	}		

	//calculate covariance entries
	for (int i = 0; i < binCount; i++) {
		for (int j = 0; j < binCount; j++) {				
			temp[0] = 0.0;	
			temp[1] = 0.0;
			temp[2] = 0.0;
			for (int sim = 0; sim < sim_number; sim++) {				
				temp[0] += ((multipoles[sim][i][0] - averages_multipoles[i][0]) * (multipoles[sim][j][0] - averages_multipoles[j][0]));
			    temp[1] += ((multipoles[sim][i][0] - averages_multipoles[i][0]) * (multipoles[sim][j][1] - averages_multipoles[j][1]));
			    temp[2]	+= ((multipoles[sim][i][1] - averages_multipoles[i][1]) * (multipoles[sim][j][1] - averages_multipoles[j][1]));				   	
			}		
			cov_matrix[0][i][j] = temp[0] / (double) sim_number;
			cov_matrix[1][i][j] = temp[1] / (double) sim_number;
			cov_matrix[2][i][j] = temp[2] / (double) sim_number;
		}	
	}		

	//normalize by STD DEVS
	double sigma_j, sigma_k;
	for (int poleOrder = 0; poleOrder < 3; poleOrder++) {
		for (int j = 0; j < binCount; j++) {
			sigma_j = sqrt(fabs(cov_matrix[poleOrder][j][j]));
			for (int k = 0; k < binCount; k++) {				
				sigma_k = sqrt(fabs(cov_matrix[poleOrder][k][k]));				
				if (sigma_j == 0.0 || sigma_k == 0.0) {
					printf("zero variance error in cov matrices, [%d][%d][%d]\n", poleOrder, j, k);
					return 1;
				} else {
					cov_matrix_final[poleOrder][j][k] = cov_matrix[poleOrder][j][k] / (sigma_j * sigma_k);
				}
			}
		}	
	}	
	
	//print to files, 0 - mono/mono, 1 - mono/quadro, 2 - quadro/quadro
	for (int i = 0; i < 3; i++) {
		sprintf(out_filename, outPath, i);
		outFile = fopen(out_filename, "w");		
		for (int j = 0; j < binCount; j++) {
			for (int k = 0; k < binCount; k++) {
				fprintf(outFile, "%f \t", cov_matrix_final[i][j][k]);
			}
			fprintf(outFile, "\n");
		}	
		fclose(outFile);
	}
	return 0;
}
