//generic function to create struct with spline coefficients, x, y values, xmin, xmax, lines
Spline* input_spline_values(int lines, double* xs, double* ys, int model) {
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
	return extend_spline(toReturn, model);	
}

Spline* extend_spline_model(SplineInfo *spline) {
	
	Spline *model_spline = malloc(sizeof(*model_spline));
	double pk_loA, pk_hiA, pk_lon, pk_hin, xmin, xmax, xend_min, xstart_max;
		
	xmin = spline->xmin;
	xmax = spline->xmax;
	xend_min = exp(log(xmin) + 0.2 * (log(xmax) - log(xmin))); //10% of range, can vary this
	xstart_max = xmax - 0.1 * (xmax - xmin);	
	
	powerlaw_regression(spline->lines, xmin, xend_min, 1.0, spline->x_vals, spline->y_vals, &pk_loA, &pk_lon);    
    powerlaw_regression(spline->lines, xstart_max, xmax, 1.0, spline->x_vals, spline->y_vals, &pk_hiA, &pk_hin);
    
    model_spline->splineInfo = spline;
	model_spline->pk_hiA = pk_hiA;
	model_spline->pk_loA = pk_loA;
	model_spline->pk_hin = pk_hin;
	model_spline->pk_lon = pk_lon;
	model_spline->model = MODEL;
	return model_spline;
}

Spline* extend_spline(SplineInfo *splineInfo, int model) {
	if (model == MEASURED) {
		Spline* extended_spline = malloc(sizeof(*extended_spline));
		extended_spline->splineInfo = splineInfo;
		extended_spline->pk_hiA = 0.0;
		extended_spline->pk_loA = 0.0;
		extended_spline->pk_hin = 0.0;
		extended_spline->pk_lon = 0.0;
		extended_spline->model = MEASURED;
		return extended_spline;
	} else if (model == MODEL) {
		return extend_spline_model(splineInfo);
	} else {
		printf("can only be REAL or MODEL in extend_spline \n");
		exit(0);
	}
}

//uses spline from what's available + power laws at low and high k. (all here return VP(k))
double splint_Pk_model(Spline *pk_model, double k) {
	double Pk;
	if (pk_model->model != MODEL) {
		printf("asking to splint a model, but giving a non-model extended SplineInfo\n");
		exit(0);
	} else {
		if (k < pk_model->splineInfo->xmin) {
			Pk = pk_model->pk_loA*pow(k, 3.0 + pk_model->pk_lon);		
		} else if (k > pk_model->splineInfo->xmax) {
			Pk = pk_model->pk_hiA*pow(k, 3.0 + pk_model->pk_hin);		
		} else {
			Pk = splint_generic(pk_model, k);
		}
	}
	return Pk/volume;
}

double splint_Pk(Spline *Pk_spline, double k) {
	if (Pk_spline->model == MODEL) {
		return splint_Pk_model(Pk_spline, k);
	} else {
		return splint_generic(Pk_spline, k)/volume;
	}
}

//inputs generic spline, given a file and a format reverse means y in terms of x
//format can only contain 2 counts of %lf.
Spline* input_spline_file(char inFile[], const char format[], int reverse, int model) {
	
	FILE *f = fopen(inFile, "r");
	if (f == NULL) {
		printf("wrong file in spline input: %s \n", inFile);
		exit(0);
	}
	
	double xmin, xmax;
	int lines = countLines(f);
	Spline* toReturn = malloc(sizeof(*toReturn));	
	double* x_vals = malloc((size_t)lines * sizeof(*x_vals));
	double* y_vals = malloc((size_t)lines * sizeof(*y_vals));

	for (int i = 0; i < lines; i++) {
		ph = fscanf(f, format, &x_vals[i], &y_vals[i]);
	}
	
	if (reverse == REVERSED) {
		toReturn = input_spline_values(lines, y_vals, x_vals, model);
	} else if (reverse == NORMAL) {
		toReturn = input_spline_values(lines, x_vals, y_vals, model);
	} else {
		printf("expected either NORMAL or REVERSED for input_spline_file() \n");
		exit(0);
	}
	
	free(x_vals);
	free(y_vals);		
	return toReturn;		
}

//sets the value at the pointer "reversed" to the reverse x = f(y) function of the spline "current"
Spline* reverse_spline(Spline *current) {
	Spline *reversed = input_spline_values(current->splineInfo->lines, current->splineInfo->y_vals, current->splineInfo->x_vals, current->model);
	return reversed;
}


//prepares a spline from a function of the type double y(double x){}
Spline* prep_spline_generic(BinInfo *binInfo, int reverse, double (*func)(double)) {
	double x_cur, y_cur;
	Spline* toReturn;
	double* xs = malloc((size_t)binInfo->bins * sizeof(*xs));
	double* ys = malloc((size_t)binInfo->bins * sizeof(*ys));	
	
	for (int i = 0; i < binInfo->bins; i++) {
		x_cur = bin_to_x(binInfo, i);
		y_cur = (*func)(x_cur);
		xs[i] = x_cur;
		ys[i] = y_cur;
	}	
	
	if (reverse == REVERSED) {
		toReturn = input_spline_values(binInfo->bins, ys, xs, MEASURED);
	} else if (reverse == NORMAL) {
		toReturn = input_spline_values(binInfo->bins, xs, ys, MEASURED);
	} else {
		printf("prep_spline_generic expected REVERSED or NORMAL \n");
		exit(0);
	}	
	return toReturn;
}


int print_spline(Spline *spline, BinInfo *binInfo) {
	for (int i = 0; i < binInfo->bins; i++) {
		printf("x: %lf, y: %lf \n", bin_to_x(binInfo, i), splint_generic(spline, bin_to_x(binInfo, i))); 
	}
	return 0;
}

int print_spline_file(Spline *spline, BinInfo *binInfo, char outFile[]) {
	File *f = fopen(outFile, "w");
	for (int i = 0; i < binInfo->bins; i++) {
		fprintf(f "%lf \t %le \n", bin_to_x(binInfo, i), splint_generic(spline, bin_to_x(binInfo, i))); 
	}
	fclose(f);
	return 0;

}
