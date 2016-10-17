Particle_Catalogue* input_catalogue_file(const char filename[], int skipLines, const char format[]) {

	FILE *catalogue_input = fopen(filename, "r");
	Particle_Catalogue* this_cat;
	
	if (catalogue_input == NULL) {
		printf("\nerror opening %s in input_catalogue_file()\n", filename);
		exit(0);
	} else {
		//count lines
		int num_lines = countLines(catalogue_input);
	
		//ALLOCATE MEMORY TO PARTICLES
		particle_no = num_lines - skipLines;
		Particle* particles = initParticles(particle_no);		
		
		//skip first [skipLines] lines, due to the format of file
		for (int i = 0; i < skipLines; i++) {
			ph = fscanf(catalogue_input, "%*[^\n]\n");
		}
		
		//create catalogue
		this_cat = malloc(sizeof(*this_cat));
		this_cat->particle_no = particle_no;
		
		//x, y, z in (Mpc/h), then vx, vy, vz in (km/s), then mass in (solar masses); \\\\\\%* ignore flag		
		for (int i = 0; i < num_lines-skipLines; i++) {
			Particle *particle = malloc(sizeof(*particle));

			ph = fscanf(catalogue_input, format, &(particle->x), &(particle->y), &(particle->z),
			&(particle->vx), &(particle->vy), &(particle->vz), &(particle->mass));
			if (ph != 7) {
				printf("format in input catalogue mismatch? Read %d values in row %d \n", ph, i);
				exit(0);
			}
			
			particles[i] = *particle;	
		}
		
		this_cat->particles = particles;
		
		printf("\ncatalogue from %s input successfully \n", filename);	
	}
	
	fclose(catalogue_input);
	return this_cat;
}
