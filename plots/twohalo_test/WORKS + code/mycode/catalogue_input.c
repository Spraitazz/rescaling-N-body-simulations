int read_full_xyz(const char filename[], int skipLines, const char format[]) {

	FILE *catalogue_input = fopen(filename, "r");
	
	if (catalogue_input == NULL) {
		printf("error opening %s\n", filename);
		return 1;
	} else {
		//count lines
		int num_lines = countLines(catalogue_input);
	
		//ALLOCATE MEMORY TO PARTICLES
		particle_no = num_lines - skipLines;
		initParticles();
		
		//skip first [skipLines] lines, due to the format of file
		for (int i = 0; i < skipLines; i++) {
			ph = fscanf(catalogue_input, "%*[^\n]\n");
		}
		
		//x, y, z in Mpc/h
		//vx, vy, vz in km/s, \\\\\\%* ignore flag
		//mass in solar masses, 	
		for (int i = 0; i < num_lines-skipLines; i++) {
			ph = fscanf(catalogue_input, format, &particles[i][0], &particles[i][1], &particles[i][2],
			&particles[i][3], &particles[i][4], &particles[i][5], &particles[i][6]);	
		}	
	}
	
	fclose(catalogue_input);
	return 0;
}
