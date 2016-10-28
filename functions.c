double randomDouble(double min, double max) {	
	return (gsl_rng_uniform(gsl_ran_r) * (max - min)) + min;
}

//Ensures PBC given current volume
int PBC(double *xx, double *yy, double *zz) {
	while (*xx < 0.0) {
		*xx += volume_limits[0];		
	}
	while (*xx >= volume_limits[0]) {
		*xx -= volume_limits[0];
	}
	while (*yy < 0.0) {
		*yy += volume_limits[1];
	}
	while (*yy >= volume_limits[1]) {
		*yy -= volume_limits[1];
	}
	while (*zz < 0.0) {
		*zz += volume_limits[2];
	}
	while (*zz >= volume_limits[2]) {
		*zz -= volume_limits[2];
	}
	return 0;
}

bool equalsDouble(double x, double y) {
	if (x - 1e-10 < y && x + 1e-10 > y) {
		return true;
	} else {
		return false;
	}

}

//Does what it says on the tin, file starts rewinded? :( can find current line in function, revert afterwards :)
int countLines(FILE *f) {
	rewind(f);
	int num_lines = 0;
	int ch;	
	while(!feof(f)) {
		ch = fgetc(f);
		if(ch == '\n') {
			num_lines++;
		}
	}	
	rewind(f);
	return num_lines;
}

//to match fftw3 convention, using row-major format
//[x][y][z] = [z + cells_z * y + cells_z * cells_y * x]
//http://www.fftw.org/doc/Row_002dmajor-Format.html#Row_002dmajor-Format
//ONLY FOR 3D ARRAY OF SIZE cells_x * cells_y * cells_z
int arr_ind(int xx, int yy, int zz) {
	return zz + cells[2] * yy + cells[2] * cells[1] * xx;
}

//same as above, but for displacements array
int arr_ind_displ(int xx, int yy, int zz) {
	return zz + cells_displ[2] * yy + cells_displ[2] * cells_displ[1] * xx;
}

int NaN_check_double(double x, char err[]) {
	if (isnan(x) || x != x) {
		printf("NaN value. %s \n", err);
	}
	return 0;
}

char* prep_path(char tail[]) {
	char* toReturn = malloc(100 * sizeof(*toReturn));
	sprintf(toReturn, "%s%s", home_directory, tail);
	return toReturn;
}

//finds the min and max mass of a catalogue, places in global variables
int mass_min_max(Particle_Catalogue *catalogue) {
	Mmax = Mmin = catalogue->particles[0].mass;
	for (int i = 1; i < catalogue->particle_no; i++) {
		if (catalogue->particles[i].mass < Mmin) {
			Mmin = catalogue->particles[i].mass;
		} else if (catalogue->particles[i].mass > Mmax) {
			Mmax = catalogue->particles[i].mass;
		}
	}
	return 0;
}

double z_to_a(double z) {
	return 1.0/(z + 1.0); 
}

double a_to_z(double a) {
	return 1.0/a - 1.0;
}

double H_a(double a, Parameters *params) {
	return params->H0*sqrt(params->omega_v_0 + params->omega_m_0*pow(a,-3.0) + params->omega_r_0*pow(a,-4.0) + (1.0 - params->omega)*pow(a,-2.0));	
}

double rho_bar_z(Parameters *params) {
	//rintf("rho cosm: %lf, z: %lf, a: %lf, omega m: %lf \n", rho_cosm, params->z, z_to_a(params->z), params->omega_m_0);
	return rho_cosm * pow(z_to_a(params->z), -3.0) * params->omega_m_0;
}

double m_to_R(double mass, Parameters *params) {
	return pow(3.0*mass/(4.0*pi*rho_bar_z(params)), 1.0/3.0);
}

double R_to_m(double R, Parameters *params) {
	//printf("R: %lf, rho bar: %lf, om m : %lf \n", R, rho_bar_z(params), params->omega_m_0);
	return (4.0/3.0)*pi*pow(R,3.0)*rho_bar_z(params);
}

double omega_m_z(double z, Parameters *params) {
	double a_cur = z_to_a(z);
	double H = H_a(a_cur, params);
	return params->omega_m_0*pow(a_cur, -3.0)*pow(params->H0/H, 2.0);
}

//Bullock et al. 2006?? or 2004?
double del_nl_z(double z, Parameters *params) {
	double om_m = omega_m_z(z, params);
	return (18.0*pi*pi + 82.0*(om_m - 1.0) - 39.0*pow(om_m - 1.0, 2.0)) / om_m; 
}

//a single struct to pass around binning stuff
BinInfo* prep_bins(double xmin, double xmax, int bins, int bin_type) {
	BinInfo* toReturn = malloc(sizeof(*toReturn));
	toReturn->bins = bins;
	if (bin_type == LOG_BIN) {
		double logxmin = log(xmin);
		double logxmax = log(xmax);
		double logdx = (logxmax - logxmin) / (double)(bins-1);
		toReturn->xmin = logxmin;
		toReturn->xmax = logxmax;
		toReturn->dx = logdx;
	} else if (bin_type == NORMAL_BIN) {
		double dx = (xmax - xmin) / (double)(bins-1);
		toReturn->xmin = xmin;
		toReturn->xmax = xmax;
		toReturn->dx = dx;
	} else {
		printf("bin types can only be NORMAL_BIN and LOG_BIN in prep_bins() \n");
		exit(0);
	}
	toReturn->bin_type = bin_type;
	return toReturn;
}

int x_to_bin(BinInfo *bin_info, double x) {
	int index;
	if (bin_info->bin_type == LOG_BIN) {
		index = (int) floor((log(x) - bin_info->xmin) / bin_info->dx);
	} else {
		index = (int) floor((x - bin_info->xmin) / bin_info->dx);
	}
	return index;
}

double bin_to_x(BinInfo *bin_info, int index) {
	double value = (double)index * bin_info->dx + bin_info->xmin;
	if (bin_info->bin_type == LOG_BIN) value = exp(value);	
	return value;
}




int print_anything(char outFile[], double xmin, double xmax, int bins, double (*y_x)(double)) {
	double x_cur, y_cur, dx;
	dx = (xmax - xmin) / (double)bins;
	FILE *f = fopen(outFile, "w");
	for (int i = 0; i <= bins; i++) {
		x_cur = xmin + (double)i *dx;
		y_cur = (*y_x)(x_cur);
		fprintf(f, "%le \t %le \n", x_cur, y_cur);	
	}
	fclose(f);
	return 0;
}

void prep_out_path(char **holder, char folder[], char subfolder[], char filename[]) {
	char final_str[200];
	//sprintf(final_str, "%s%s%s

}

double splint_generic(Spline *spline, double x) {
	double y, margin;
	SplineInfo *info = spline->splineInfo;
	margin = (info->xmax - info->xmin)*1e-6;
	//allowing small deviation to counter bugs with the result ending up 0 inside the range
	if (info->xmin - margin < x && x <= info->xmin) {
		splint(info->x_vals, info->y_vals, info->coeffs, info->lines, info->xmin, &y);
	} else if (info->xmax <= x && x < info->xmax + margin) {
		splint(info->x_vals, info->y_vals, info->coeffs, info->lines, info->xmax, &y);
	} else if (x <= info->xmin - margin || x >= info->xmax + margin) {
		y = 0.0;
	} else {
		splint(info->x_vals, info->y_vals, info->coeffs, info->lines, x, &y);
	}
	return y;
}

int print_delsq(Spline *Pk_spline, char filename[], BinInfo* binInfo, bool scaled) {
	
	FILE *f = fopen(filename, "w");	
	double delsq, k_cur;
		
	for (int i = 0; i < binInfo->bins; i++) {
		k_cur = bin_to_x(binInfo, i);		
		if (scaled) {
			delsq = delsq_lin(Pk_spline, k_cur*s);
		} else {
			delsq = delsq_lin(Pk_spline, k_cur);
		}
		fprintf(f, "%lf \t %le \n", k_cur, delsq);
	}
	
	fclose(f);
	return 0;
}

