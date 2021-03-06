Particle_Catalogue* populateHaloes_oneMass(int halo_no, double halo_mass) {	
	Particle_Catalogue* toReturn = malloc(sizeof(*toReturn));
	toReturn->particle_no = halo_no;
	Particle* particles = initParticles(halo_no);
	for (int i = 0; i < halo_no; i++) {	
		//positions		
		particles[i].x = randomDouble(0.0, volume_limits[0]);	
		particles[i].y = randomDouble(0.0, volume_limits[1]);
		particles[i].z = randomDouble(0.0, volume_limits[2]);
		
		// velocities
		particles[i].vx = 0.0;	
		particles[i].vy = 0.0;
		particles[i].vz = 0.0;	
		
		// mass
		particles[i].mass = halo_mass;
	}	
	toReturn->particles = particles;			
	return toReturn;
}

Particle_Catalogue* combine_catalogues(Particle_Catalogue* cat1, Particle_Catalogue* cat2) {
	Particle_Catalogue* toReturn = malloc(sizeof(*toReturn));
	Particle* particles = initParticles(cat1->particle_no + cat2->particle_no);
	
	for (int i = 0; i < cat1->particle_no; i++) {
		memcpy(&(particles[i]), &(cat1->particles[i]), sizeof(particles[i]));	
	}
	for (int i = cat1->particle_no; i < cat1->particle_no + cat2->particle_no; i++) {
		memcpy(&(particles[i]), &(cat2->particles[i - cat1->particle_no]), sizeof(particles[i]));	
	}
	toReturn->particles = particles;
	return toReturn;
}



Particle_Catalogue* populateHaloes_twoMasses(double smallMass, double bigMass, int halo_no, int no_big) {
	Particle_Catalogue* toReturn = malloc(sizeof(*toReturn));
	toReturn->particle_no = halo_no;
	Particle* particles = initParticles(halo_no);	
	for (int i = 0; i < halo_no; i++) {	
		//positions		
		particles[i].x = randomDouble(0.0, volume_limits[0]);	
		particles[i].y = randomDouble(0.0, volume_limits[1]);
		particles[i].z = randomDouble(0.0, volume_limits[2]);
		
		// velocities
		particles[i].vx = 0.0;	
		particles[i].vy = 0.0;
		particles[i].vz = 0.0;	
		
		// mass
		if (i < no_big) {
			particles[i].mass = bigMass;
		} else {
			particles[i].mass = smallMass;
		}
	}		
	toReturn->particles = particles;			
	return toReturn;	
}



Particle_Catalogue* populateParticles_crystalLattice(int particle_no) {
	if (particle_no != cells[0]*cells[1]*cells[2]) {
		printf("for lattice particle number should equal cell number\n");
		exit(0);
	}	
	double cellSizes[3];
	cellSizes[0] = volume_limits[0] / (double) cells[0];
	cellSizes[1] = volume_limits[1] / (double) cells[1];
	cellSizes[2] = volume_limits[2] / (double) cells[2];
	
	Particle_Catalogue* toReturn = malloc(sizeof(*toReturn));
	toReturn->particle_no = particle_no;
	Particle* particles = initParticles(particle_no);	
	
	//placing particle in center of cell
	for (int i = 0; i < cells[0]; i++) {
		for (int j = 0; j < cells[1]; j++) {
			for (int k = 0; k < cells[2]; k++) {

				particles[arr_ind(i,j,k)].x = ((double)i + 0.5) * cellSizes[0];
				particles[arr_ind(i,j,k)].y = ((double)j + 0.5) * cellSizes[1];
				particles[arr_ind(i,j,k)].z = ((double)k + 0.5) * cellSizes[2];
								
				// velocities
				particles[arr_ind(i,j,k)].vx = 0.0;	
				particles[arr_ind(i,j,k)].vy = 0.0;
				particles[arr_ind(i,j,k)].vz = 0.0;	
		
				// mass
				particles[arr_ind(i,j,k)].mass = 0.0;
			}
		}
	}
	toReturn->particles = particles;
	return toReturn;
}

/*

int populateHaloes_diffMass_lattice(double smallMass, double bigMass, int no_big) {	

	if (particle_no != cells[0]*cells[1]*cells[2]) {
		printf("for lattice particle number should equal cell number\n");
		exit(0);
	}	
	
	double cellSizes[3];
	cellSizes[0] = volume_limits[0] / (double) cells[0];
	cellSizes[1] = volume_limits[1] / (double) cells[1];
	cellSizes[2] = volume_limits[2] / (double) cells[2];
	int p_no = 0;
	
	//placing particle in center of cell
	for (int i = 0; i < cells[0]; i++) {
		for (int j = 0; j < cells[1]; j++) {
			for (int k = 0; k < cells[2]; k++) {
				particles[arr_ind(i,j,k)][0] = ((double)i + 0.5) * cellSizes[0];
				particles[arr_ind(i,j,k)][1] = ((double)j + 0.5) * cellSizes[1];
				particles[arr_ind(i,j,k)][2] = ((double)k + 0.5) * cellSizes[2];
								
				// velocities
				particles[arr_ind(i,j,k)][3] = 0.0;	
				particles[arr_ind(i,j,k)][4] = 0.0;
				particles[arr_ind(i,j,k)][5] = 0.0;			
				
				// mass
				if (p_no < no_big) {
					particles[arr_ind(i,j,k)][6] = bigMass;
				} else {
					particles[arr_ind(i,j,k)][6] = smallMass;
				}
				
				p_no += 1;

			}
		}
	}
	return 0;
}

int populateHaloes_randomMass(double M_min, double M_max) {	
	for (int i = 0; i < halo_no; i++) {	
		//positions		
		particles[i][0] = randomDouble(0.0, volume_limits[0]);	
		particles[i][1] = randomDouble(0.0, volume_limits[1]);
		particles[i][2] = randomDouble(0.0, volume_limits[2]);
		
		// velocities
		particles[i][3] = 0.0;	
		particles[i][4] = 0.0;
		particles[i][5] = 0.0;		
		
		double mass = pow(10.0, 12.0 + gsl_ran_gaussian(gsl_ran_r, 0.5));			
		
		particles[i][6] = mass;
		
	}				
	return 0;
}



int populateHaloes_catalogueMass(char inFile[]) {
	FILE *catalogue_input = fopen(inFile, "r");
	int num_lines = countLines(catalogue_input);
	
	particle_no = num_lines-1;
	initParticles();
	
	ph = fscanf(catalogue_input, "%*[^\n]\n");
	printf("adding %d particles \n", particle_no);
	
	for (int i = 0; i < particle_no; i++) {
		//positions		
		particles[i][0] = randomDouble(0.0, volume_limits[0]);	
		particles[i][1] = randomDouble(0.0, volume_limits[1]);
		particles[i][2] = randomDouble(0.0, volume_limits[2]);
		
		// velocities
		particles[i][3] = 0.0;	
		particles[i][4] = 0.0;
		particles[i][5] = 0.0;	
		
		// mass
		ph = fscanf(catalogue_input, "%*f \t %*f \t %*f \t %*f \t %*f \t %*f \t %le", &particles[i][6]);	
	}	

	mass_min_max();
}
*/

