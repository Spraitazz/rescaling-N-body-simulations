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

double H_a(double a) {
	return H0*sqrt(omega_v_0 + omega_m_0*pow(a,-3.0) + omega_r_0*pow(a,-4.0) + (1.0 - omega)*pow(a,-2.0));	
}

double rho_bar_z(double z) {
	return rho_cosm*pow(z_to_a(z), -3.0)*omega_m_0;
}

double m_to_R(double mass, double z) {
	return pow(3.0*mass/(4.0*M_PI*rho_bar_z(z)), 1.0/3.0);
}

double R_to_m(double R, double z) {
	return (4.0/3.0)*pi*pow(R,3.0)*rho_bar_z(z);
}

double omega_m_z(double z) {
	double a_cur = z_to_a(z);
	double H = H_a(a_cur);
	return omega_m_0*pow(a_cur, -3.0)*pow(H0/H, 2.0);
}

//Bullock et al. 2006?? or 2004?
double del_nl_z(double z) {
	return (18.0*pi*pi + 82.0*(omega_m_z(z) - 1.0) - 39.0*pow(omega_m_z(z) - 1.0, 2.0)) / omega_m_z(z); 
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

//generic function to create struct with spline coefficients, x, y values, xmin, xmax, lines
SplineInfo* input_spline_values(int lines, double* xs, double* ys) {
	double xmin, xmax;	
	bool x_mono, x_inc;
	double* x_vals = malloc(lines * sizeof(*xs));
	double* y_vals = malloc(lines * sizeof(*ys));
	memcpy(x_vals, xs, lines * sizeof(*xs));
	memcpy(y_vals, ys, lines * sizeof(*ys));	
	double* coeffs = malloc((size_t)lines * sizeof(*coeffs));	
	int dupcount = 0;
	x_mono = true;	
	x_inc = false;	
	
	//monotoneity check for x
	for (int i = 0; i < lines-1; i++) {
		if (x_vals[i+1] > x_vals[i]) {
			x_inc = true;
		} else if (x_vals[i+1] < x_vals[i]) {
			//x decreased going forward, check monotonous
			if (x_inc == true) {
				x_mono = false;
				printf("the x-values of the function given are not monotonous\n");
				printf("%le followed by %le \n", x_vals[i], x_vals[i+1]);
				exit(0);
			}
		} else {
			//x stayed the same -> nonsense??
			if (dupcount < 2) {
				printf("duplicate x-values in function given, check index %d\n", i);
				dupcount += 1;	
			}			
		}
	}
		
	if (x_inc) {
		xmin = x_vals[0];
		xmax = x_vals[lines-1];
	} else {
		xmin = x_vals[lines-1];
		xmax = x_vals[0];
		//also reverse, spline/splint wants x in increasing order
		double temp, temp2;		
		for (int i = 0; i < lines/2; i++) {
			temp = x_vals[i];
			temp2 = y_vals[i];
			x_vals[i] = x_vals[lines - i - 1];
			y_vals[i] = y_vals[lines - i - 1];
			x_vals[lines - i - 1] = temp;
			y_vals[lines - i - 1] = temp2;
		}	
	}	
	
	//spline to get the coefficients
	spline(x_vals, y_vals, lines, 1.0e31, 1.0e31, coeffs);	
	//store spline information, together with coefficients, and return it	
	SplineInfo* toReturn = malloc(sizeof(*toReturn));
	toReturn->x_vals = x_vals;
	toReturn->y_vals = y_vals;
	toReturn->coeffs = coeffs;
	toReturn->lines = lines;
	toReturn->xmin = xmin;
	toReturn->xmax = xmax;
	return toReturn;		
}

//inputs generic spline, given a file and a format reverse means y in terms of x
//format can only contain 2 counts of %lf.
SplineInfo* input_spline_file(char inFile[], const char format[], int reverse) {
	
	FILE *f = fopen(inFile, "r");
	if (f == NULL) {
		printf("wrong file in spline input: %s \n", inFile);
		exit(0);
	}
	
	double xmin, xmax;
	int lines = countLines(f);
	SplineInfo* toReturn = malloc(sizeof(*toReturn));	
	double* x_vals = malloc((size_t)lines * sizeof(*x_vals));
	double* y_vals = malloc((size_t)lines * sizeof(*y_vals));

	for (int i = 0; i < lines; i++) {
		ph = fscanf(f, format, &x_vals[i], &y_vals[i]);
	}
	
	if (reverse == REVERSED) {
		toReturn = input_spline_values(lines, y_vals, x_vals);
	} else if (reverse == NORMAL) {
		toReturn = input_spline_values(lines, x_vals, y_vals);
	} else {
		printf("expected either NORMAL or REVERSED for input_spline_file() \n");
		exit(0);
	}
	
	free(x_vals);
	free(y_vals);		
	return toReturn;		
}

//sets the value at the pointer "reversed" to the reverse x = f(y) function of the spline "current"
SplineInfo* reverse_spline(SplineInfo *current) {
	SplineInfo *reversed = input_spline_values(current->lines, current->y_vals, current->x_vals);
	return reversed;
}

//prepares a spline from a function of the type double y(double x){}
SplineInfo* prep_spline_generic(BinInfo *binInfo, int reverse, double (*func)(double)) {
	double x_cur, y_cur;
	SplineInfo* toReturn;
	double* xs = malloc((size_t)binInfo->bins * sizeof(*xs));
	double* ys = malloc((size_t)binInfo->bins * sizeof(*ys));	
	
	for (int i = 0; i < binInfo->bins; i++) {
		x_cur = bin_to_x(binInfo, i);
		y_cur = (*func)(x_cur);
		xs[i] = x_cur;
		ys[i] = y_cur;
	}	
	
	if (reverse == REVERSED) {
		toReturn = input_spline_values(binInfo->bins, ys, xs);
	} else if (reverse == NORMAL) {
		toReturn = input_spline_values(binInfo->bins, xs, ys);
	} else {
		printf("prep_spline_generic expected REVERSED or NORMAL \n");
		exit(0);
	}	
	return toReturn;
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

double splint_generic(SplineInfo *info, double x) {
	double y, margin;
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

int print_delsq(SplineInfo_Extended *pk_spline, char filename[], double kmin, double kmax, int bins, bool scaled) {
	FILE *delsq_out_file = fopen(filename, "w");	
	double delsq, k_cur, logkmin, logdk, logk;
	logkmin = log10(kmin);
	logdk = (log10(kmax) - logkmin) / (double)(bins - 1);
	for (int i = 0; i < bins; i++) {
		logk = logkmin + (double)i * logdk;
		k_cur = pow(10.0, logk);		
		if (scaled) {
			delsq = delsq_lin(pk_spline, k_cur*s);
		} else {
			delsq = delsq_lin(pk_spline, k_cur);
		}
		fprintf(delsq_out_file, "%lf \t %le \n", k_cur, delsq);
	}
	fclose(delsq_out_file);
	return 0;
}

