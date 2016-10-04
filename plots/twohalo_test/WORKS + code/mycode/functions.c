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

//finds the min and max mass of a catalogue, places in global variables
int mass_min_max() {
	Mmin = particles[0][6];
	Mmax = particles[0][6];
	for (int i = 1; i < particle_no; i++) {
		if (particles[i][6] < Mmin) {
			Mmin = particles[i][6];
		} else if (particles[i][6] > Mmax) {
			Mmax = particles[i][6];
		}
	}
	return 0;
}

double rho_bar_z(double z) {
	return rho_cosm*z_to_a(z)*omega_m;
}

double m_to_R(double mass, double z) {
	return pow(3.0*mass/(4.0*M_PI*rho_bar_z(z)), 1.0/3.0);
}

double R_to_m(double R, double z) {
	return (4.0/3.0)*pi*pow(R,3.0)*rho_bar_z(z);
}

//a single struct to pass around binning stuff
BinInfo prep_bins(double xmin, double xmax, int bins, bool logbin) {
	BinInfo toReturn;
	toReturn.bins = bins;
	if (logbin) {
		double logxmin = log(xmin);
		double logxmax = log(xmax);
		double logdx = (logxmax - logxmin) / (double)(bins-1);
		toReturn.xmin = logxmin;
		toReturn.xmax = logxmax;
		toReturn.dx = logdx;
		toReturn.logbin = true;	
	} else {
		double dx = (xmax - xmin) / (double)(bins-1);
		toReturn.xmin = xmin;
		toReturn.xmax = xmax;
		toReturn.dx = dx;
		toReturn.logbin = false;
	}
	return toReturn;
}

int x_to_bin(BinInfo *bin_info, double x) {
	int index;
	if (bin_info->logbin) {
		index = (int) floor((log(x) - bin_info->xmin) / bin_info->dx);
	} else {
		index = (int) floor((x - bin_info->xmin) / bin_info->dx);
	}
	return index;
}

double bin_to_x(BinInfo *bin_info, int index) {
	double value = (double)index * bin_info->dx + bin_info->xmin;
	if (bin_info->logbin) value = exp(value);	
	return value;
}

//generic function to create struct with spline coefficients, x, y values, xmin, xmax, lines
SplineInfo input_spline_values(int lines, double* x_vals, double* y_vals) {
	double xmin, xmax;	
	bool x_mono, x_inc;
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
	SplineInfo toReturn = {x_vals, y_vals, coeffs, lines, xmin, xmax};	
	return toReturn;		
}

//inputs generic spline, given a file and a format reverse means y in terms of x
//format can only contain 2 counts of %lf.
SplineInfo input_spline_file(char inFile[], const char format[], bool reverse) {	
	FILE *f = fopen(inFile, "r");
	if (f == NULL) {
		printf("wrong file in spline input\n");
		exit(0);
	}
	double xmin, xmax;
	int lines = countLines(f);
	SplineInfo toReturn;	
	double* x_vals = malloc((size_t)lines * sizeof(*x_vals));
	double* y_vals = malloc((size_t)lines * sizeof(*y_vals));
	double* coeffs = calloc((size_t)lines, sizeof(*coeffs));

	for (int i = 0; i < lines; i++) {
		ph = fscanf(f, format, &x_vals[i], &y_vals[i]);
	}
	
	if (reverse) {
		toReturn = input_spline_values(lines, y_vals, x_vals);
	} else {
		toReturn = input_spline_values(lines, x_vals, y_vals);
	}	
	return toReturn;		
}

//sets the value at the pointer "reversed" to the reverse x = f(y) function of the spline "current"
int reverse_spline(SplineInfo *current, SplineInfo *reversed) {
	*reversed = input_spline_values(current->lines, current->y_vals, current->x_vals);
	return 0;
}

//prepares a spline from a function of the type double y(double x){}
SplineInfo prep_spline_generic(double xmin, double xmax, int bins, bool logbin, bool reverse, double (*func)(double)) {
	double x_cur, y_cur, dx, logx, logdx, logxmin, logxmax;
	SplineInfo toReturn;
	double* xs = malloc((size_t)bins * sizeof(*xs));
	double* ys = malloc((size_t)bins * sizeof(*ys));	
	
	if (logbin) {
		logxmin = log(xmin);
		logxmax = log(xmax);
		logdx = (logxmax - logxmin) / (double)(bins - 1);
		for (int i = 0; i < bins; i++) {
			logx = logxmin + (double)i * logdx;
			x_cur = exp(logx);
			y_cur = (*func)(x_cur);
			xs[i] = x_cur;
			ys[i] = y_cur;
		}	
	} else {	
		dx = (xmax - xmin) / (double)(bins - 1);	
		for (int i = 0; i < bins; i++) {		
			x_cur = xmin + (double)i * dx;			
			y_cur = (*func)(x_cur);
			xs[i] = x_cur;
			ys[i] = y_cur;				
		}
	}
	if (reverse) {
		toReturn = input_spline_values(bins, ys, xs);
	} else {
		toReturn = input_spline_values(bins, xs, ys);
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
	double y;
	if (x < info->xmin || x > info->xmax) {
		y = 0.0;
	} else {
		splint(info->x_vals, info->y_vals, info->coeffs, info->lines, x, &y);
	}
	return y;
}

double splint_generic_VPk(SplineInfo *info, double x) {
	double y;
	if ((x < info->xmin) || (x > info->xmax)) {
		y = 0.0;
	} else {
		splint(info->x_vals, info->y_vals, info->coeffs, info->lines, x, &y);
	}
	return y/volume;
}
/*
double splint_generic_extended(SplineInfo info, double x) {
	double y;
	int bin1, bin2;
	for (int i = 0; i < info.lines-1; i++) {
		if (info.x_vals[i] < x && x < info.x_vals[i+1]) {
		//	bin1 = 
		} else if (info.x_vals[i] > x && x > info.x_vals[i+1]) {
		
		}
	}
	if (x < info.xmin) {
		splint(info.x_vals, info.y_vals, info.coeffs, info.lines, info.xmin, &y);
	} else if (x > info.xmax) {
		splint(info.x_vals, info.y_vals, info.coeffs, info.lines, info.xmax, &y);
	} else {
		splint_generic(info, x);
	}
	return y;
}
*/



/*
typedef struct {
	double* x_vals;
	double* fx_vals;
	double y;
} func_vals_packed;

typedef struct {
	int param_no;
	int* bins;
	double** xyz_vals;
	double** fxyz_vals;
	double** coeffs;
	double (*splint)(double* param_values);
} SplineInfo_Ndim;

int test_splint_Ndim() {
	SplineInfo_Ndim sp;
	double arr[2] = {0.0, 1.0};
	sp.splint(arr);

	return 0;
}


double* test_function(double xmin, double xmax, double xbins, double y) {
	//f(x, y) = x^2 + y^2
	double dx = (xmax - xmin) / (double)xbins;
	double x_cur;
	double* results = calloc((size_t)xbins, sizeof(*results));
	for (int i = 0; i < xbins; i++) {
		x_cur = xmin + (double)i * dx;
		results[i] = x_cur*x_cur + y*y;
	}
	return results;
}
*/


//typedef double (*spline_function)(double);
//spline_function* spline_funcs_Dplus = calloc(2, sizeof(spline_function));
//spline_function[0] 
/*
typedef struct {
	int param_no;
	int* bins;
	double** xy_vals;
	double** fxy_vals;
	double** coeffs;
	double (*splint)(double* param_values);
} SplineInfo_2D;

SplineInfo_2D 2Dspline(int* bins, double** xy_vals, double** fxy_vals) {
	double** coeffs;
	//2 parameters - x and y
	//spline_coeffs(int lines, double* x_vals, double* y_vals, bool reverse)
	for (int x = 0; x < bins[0]; x++) {
		for (int y = 0; y < bins[1]; y++) {
			coeffs[x] = spline_coeffs(bins[1], xy_vals[x], 
		}
	
	}
	
	

}
*/

/*
SplineInfo_Ndim spline_Ndim(int param_no, int cur_param, int* bins, double** xyz_vals, double** fzyx_vals) {
	//n params -> n functions 
	SplineInfo_Ndim toReturn;
	double** coeffs;
	double* inputs;
	//spline_coeffs(int lines, double* x_vals, double* y_vals, bool reverse)
	for (int par = 0; par < param_no; par++) {
		//param i is binned into bins[i] bins
		inputs[par] = xyz_vals[par][
		for (int bin = 0; bin < bins[par]; bin++) {
			//if f = f(x, y), xyz_vals[0][bin] gives f(x = xbin, y = y0)
			coeffs[par] = spline_coeffs(bins[par], xyz_vals[par]);
		
		}
		
	
	}


	return toReturn;
}*/

//takes away noise power from a given spectrum, makes new file
//NOT USED
/*
int minus_noise(char full_spectrum[], char noise_spectrum[], char outFile[]) {
	FILE *spectrum = fopen(full_spectrum, "r");
	FILE *noise = fopen(noise_spectrum, "r");
	FILE *out = fopen(outFile, "w");
	int lines = countLines(spectrum);
	rewind(spectrum);
	double k, k_noise, mono, mono_noise, quadro, quadro_noise;
	int bins, bins_noise;
	for (int i = 0; i < lines; i++) {
		fscanf(spectrum, "%lf \t %lf \t %lf \t %d\n", &k, &mono, &quadro, &bins);
		fscanf(noise, "%lf \t %lf \t %lf \t %d\n", &k_noise, &mono_noise, &quadro_noise, &bins_noise);
		fprintf(out, "%.5lf \t %.15lf \t %.15lf \t %d\n", k, mono-mono_noise, quadro, bins);
	}
	fclose(spectrum);
	fclose(noise);
	fclose(out);	
	return 0;
} 
*/
/*
SplineInfo input_spline_variance(char inFile[], const char format[], bool reverse) {
	double xmin, xmax;
	int lines;
	SplineInfo toReturn;	
	FILE *f = fopen(inFile, "r");	
	if (f == NULL) {
		printf("wrong file in spline input\n");
		exit(0);
	}
	lines = countLines(f);	
	double* x_vals = calloc((size_t)lines, sizeof(*x_vals));
	double* y_vals = calloc((size_t)lines, sizeof(*y_vals));
	double* coeffs = calloc((size_t)lines, sizeof(*coeffs));	

	for (int i = 0; i < lines; i++) {
		ph = fscanf(f, format, &x_vals[i], &y_vals[i]);
	}
	
	for (int i = 0; i < lines; i++) {
		y_vals[i] *= Dplus_test; //TEST MODE - becomes 0.25*sigma_squared at each r, should return z ~ 1.2
	}
	
	if (reverse) {
		toReturn = input_spline_values(lines, y_vals, x_vals);
	} else {
		toReturn = input_spline_values(lines, x_vals, y_vals);
	}
	
	return toReturn;		
}
*/

/*
double* spline_coeffs(int lines, double* x_vals, double* y_vals, bool reverse) {
	double xmin, xmax;
	double* coeffs = calloc((size_t)lines, sizeof(*coeffs));
	bool x_mono, y_mono, x_inc, y_inc;
	int nonmono_ind;
	x_mono = true;	
	y_mono = true;
	x_inc = false;
	y_inc = false;	
	
	//monotoneity check for x
	for (int i = 0; i < lines-1; i++) {
		if (x_vals[i+1] > x_vals[i]) {
			x_inc = true;
		} else if (x_vals[i+1] < x_vals[i]) {
			//x decreased going forward, check monotonous
			if (x_inc == true) {
				x_mono = false;
			}
		} else {
			//x stayed the same -> nonsense??
			printf("duplicate x-values in function given, check index %d", i);
			//exit(0);		
		}
	}
	//x MUST BE monotonous	
	if (x_mono == false && reverse == false) {
		printf("the x-values of the function given are not monotonous\n");
		exit(0);
	}
	
	//same check for y
	for (int i = 0; i < lines-1; i++) {
		if (y_vals[i+1] > y_vals[i]) {
			y_inc = true;
		} else if (y_vals[i+1] < y_vals[i]) {
			if (y_inc == true) {
				y_mono = false;
				nonmono_ind = i;
			}
		} else {
			//printf("duplicate y-values in function given\n");
		}
	}
	
	//y MUST BE Monotonous (reversed case)
	if (y_mono == false && reverse == true) {
		printf("the y-values of the function given are not monotonous, index %d\n", nonmono_ind);
		exit(0);
	}
	
	double temp, temp2;		
	if (reverse) {
		//x = f(y)
		//replace y with x values
		for (int i = 0; i < lines; i++) {
			temp = x_vals[i];
			x_vals[i] = y_vals[i];
			y_vals[i] = temp;
		}
		//make sure x-values (independent variable) are increasing, dont forget to reverse y-vals also	
		if (y_inc == false) {			
			for (int i = 0; i < lines/2; i++) {
				temp = x_vals[i];
				temp2 = y_vals[i];
				x_vals[i] = x_vals[lines - i - 1];
				y_vals[i] = y_vals[lines - i - 1];
				x_vals[lines - i - 1] = temp;
				y_vals[lines - i - 1] = temp2;
			}
		}
	} else {
		//non-reversed case ( y = f(x) )
		if (x_inc == false) {
			printf("reversing \n");
			for (int i = 0; i < lines/2; i++) {
				temp = x_vals[i];
				temp2 = y_vals[i];
				x_vals[i] = x_vals[lines - i - 1];
				y_vals[i] = y_vals[lines - i - 1];
				x_vals[lines - i - 1] = temp;
				y_vals[lines - i - 1] = temp2;
			}
		}	
	}	
	
	//here "x" is the independent variable
	xmin = x_vals[0];
	xmax = x_vals[0];
	for (int i = 0; i < lines; i++) {
		if (x_vals[i] < xmin) {
			xmin = x_vals[i];
		} else if (x_vals[i] > xmax) {
			xmax = x_vals[i];
		}
	}
	
	//spline to get the coefficients
	spline(x_vals, y_vals, lines, 1.0e31, 1.0e31, coeffs);	

	return coeffs;
}
*/
