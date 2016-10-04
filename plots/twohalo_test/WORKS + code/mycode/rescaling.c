int prep_Pk_model_spline(SplineInfo *spline, SplineInfo_Extended *model_spline) {
	double pk_loA, pk_hiA, pk_lon, pk_hin, xmin, xmax, xend_min, xstart_max;	
	xmin = spline->xmin;
	xmax = spline->xmax;
	xend_min = exp(log(xmin) + 0.2 * (log(xmax) - log(xmin))); //10% of range, can vary this
	xstart_max = xmax - 0.1 * (xmax - xmin);	
	powerlaw_regression(spline->lines, xmin, xend_min, 1.0, spline->x_vals, spline->y_vals, &pk_loA, &pk_lon);    
    powerlaw_regression(spline->lines, xstart_max, xmax, 1.0, spline->x_vals, spline->y_vals, &pk_hiA, &pk_hin);
    model_spline->spline = spline;
	model_spline->pk_hiA = pk_hiA;
	model_spline->pk_loA = pk_loA;
	model_spline->pk_hin = pk_hin;
	model_spline->pk_lon = pk_lon;
	model_spline->model = true;
	return 0;
}

int extend_spline_nonmodel(SplineInfo *spline, SplineInfo_Extended *extended_spline) {
	extended_spline->spline = spline;
	extended_spline->pk_hiA = 0.0;
	extended_spline->pk_loA = 0.0;
	extended_spline->pk_hin = 0.0;
	extended_spline->pk_lon = 0.0;
	extended_spline->model = false;
	return 0;
}

//uses spline from what's available + power laws at low and high k. (all here return VP(k))
double splint_Pk_model(SplineInfo_Extended *pk_model, double k) {
	if (pk_model->model == false) {
		printf("asking to splint a model, but giving a non-model extended SplineInfo\n");
		exit(0);
	}
	double Pk;
	if (k < pk_model->spline->xmin) {
		Pk = pk_model->pk_loA*pow(k, 3.0 + pk_model->pk_lon);		
	} else if (k > pk_model->spline->xmax) {
		Pk = pk_model->pk_hiA*pow(k, 3.0 + pk_model->pk_hin);		
	} else {
    	Pk = splint_generic(pk_model->spline, k);
	}
	//Pk actually was VP(k), need P(k)
	return Pk/volume;
}

double splint_Pk(SplineInfo_Extended *pk_spline, double k) {
	if (pk_spline->model == true) {
		return splint_Pk_model(pk_spline, k);
	} else {
		return splint_generic_VPk(pk_spline->spline, k);
	}
}

int displacements_ZA(int xyz) {

    double kx, ky, kz, k, k_sq, coef_real, coef_im, del_ratio;
    double Pk_model, Pk_real, phase, amplitude, prefactor;
    double Pk_current, Pk_target;
    double fundmodes[3];
    int index, index_neg, Index, negkIndex;  

    fundmodes[0] = 2.0*M_PI/volume_limits[0];
    fundmodes[1] = 2.0*M_PI/volume_limits[1];
    fundmodes[2] = 2.0*M_PI/volume_limits[2];	  
    
    for (int aa = 0; aa < cells_displ[0]; aa++) {
    	for (int bb = 0; bb < cells_displ[1]; bb++) {
    		for (int cc = 0; cc < cells_displ[2]; cc++) {
    			kx = aa * fundmodes[0];
    			ky = bb * fundmodes[1];
				kz = cc * fundmodes[2];  			
				
				//convention - up to N/2 positive frequencies, N/2 is fc, then -fc+ffund, ...
				if (aa > cells_displ[0]/2) kx -= ((double) cells_displ[0]) * fundmodes[0]; //or -= 2*kx_nyquist
				if (bb > cells_displ[1]/2) ky -= ((double) cells_displ[1]) * fundmodes[1];
				if (cc > cells_displ[2]/2) kz -= ((double) cells_displ[2]) * fundmodes[2];
				
				index = arr_ind_displ(aa, bb, cc);	
				k_sq = kx*kx + ky*ky + kz*kz;
				k = sqrt(k_sq);
				
				//USING REAL SIMULATION DATA/ALREADY DONE FFT OF OVERDENSITY
				/*
				overdensity_fourier[index][0] *= exp(-k_sq*R_nl*R_nl/2.0);
				overdensity_fourier[index][1] *= exp(-k_sq*R_nl*R_nl/2.0);				
				coef_real = overdensity_fourier[index][0];
				coef_im = overdensity_fourier[index][1];				
				Pk_current = splint_generic(Pk_spline_current, k);
				Pk_target = splint_generic(Pk_spline_target, k);
				if (Pk_current < 0.0 || Pk_target < 0.0) {
					printf("negative Pk values\n");
					return 1;
				}
				*/
				
				//USING MODEL	
				Pk_model = splint_Pk(&Pk_spline_current, k);						
				phase = randomDouble(0.0, 2.0*M_PI);
				amplitude = sqrt(Pk_model);
				coef_real = amplitude * cos(phase);// * exp(-k_sq*R_nl*R_nl/2.0);
				coef_im = amplitude * sin(phase);// * exp(-k_sq*R_nl*R_nl/2.0);		
								
				switch (xyz) {					
					case 0:
						displacements_fourier[index][0] = kx * coef_im / k_sq;							
						displacements_fourier[index][1] = -kx * coef_real / k_sq;
						break;
					case 1:
						displacements_fourier[index][0] = ky * coef_im / k_sq;
						displacements_fourier[index][1] = -ky* coef_real / k_sq;
						break;
					case 2:
						displacements_fourier[index][0] = kz * coef_im / k_sq;
						displacements_fourier[index][1] = -kz* coef_real / k_sq;
						break;
				}			
								
					
					/*
					if (Pk_current != 0.0) {					
						//df(k) and dv(k) here					
						del_ratio = Pk_target/Pk_current;
						prefactor = sqrt(del_ratio) - 1.0;					
						displacements_fourier[index][0] *= prefactor;	
						displacements_fourier[index][1] *= prefactor; <-------???? GOOD OR dont need im?NOT??		
					} else {
						printf("Pk current = 0.0\n");
						displacements_fourier[index][0] = 0.0;
						displacements_fourier[index][1] = 0.0;
					}	
					*/				
			}    	
    	}    
	}	
	displacements_fourier[0][0] = displacements_fourier[0][1] = 0.0;	

	// Hermitian condition. One hemi-sphere is independent, e.g. k_z >= 0.
    for (int u = cells_displ[2]-1; u >= cells_displ[2]/2; u--) {
        for (int j = 0; j < cells_displ[1]; j++) {
            for (int i = 0; i < cells_displ[0]; i++) {
		        
		        negkIndex = u*cells_displ[1]*cells_displ[0] + j*cells_displ[0] + i;
		        Index = 0;
		           
		        // zero maps to zero on reflection through the origin.
		        if(i!=0) Index += (cells_displ[0] - i);
		        if(j!=0) Index += (cells_displ[1] - j)*cells_displ[0];
		        Index += (cells_displ[2] - u)*cells_displ[1]*cells_displ[0];			
			
				displacements_fourier[negkIndex][0] = displacements_fourier[Index][0];
				displacements_fourier[negkIndex][1] = -1.0*displacements_fourier[Index][1];
				if(negkIndex == Index) displacements_fourier[Index][1] = 0.0; // purely real  
		    		        
			}
		}
	}
   
    // in the k_z=0 plane one semi-circle is independent, k_y>0.        
    for(int j = cells_displ[1]-1; j >= cells_displ[1]/2; j--){
        for(int i = 0; i < cells_displ[0]; i++){
            negkIndex = j*cells_displ[0] + i;                      
            Index = 0;
               
            // zero maps to zero on reflection through the origin.
            if(i!=0)  Index  += (cells_displ[0] - i);
            Index  += (cells_displ[1] - j)*cells_displ[0];

            displacements_fourier[negkIndex][0]  = displacements_fourier[Index][0];
            displacements_fourier[negkIndex][1]  = -1.0*displacements_fourier[Index][1];
            if(negkIndex == Index) displacements_fourier[Index][1] = 0.0;           
        }
    }
   
    // on the line k_z=k_y=0, one half is independent, k_x>=0.
    for(int i = cells_displ[0]-1; i >= cells_displ[0]/2; i--){
        negkIndex = i;                      
        Index = 0;
               
        // zero maps to zero on reflection through the origin.
        Index += (cells_displ[0] - i);

        displacements_fourier[negkIndex][0] = displacements_fourier[Index][0];
        displacements_fourier[negkIndex][1] = -1.0*displacements_fourier[Index][1];
        if(negkIndex == Index) displacements_fourier[Index][1] = 0.0;       
    }    	
	return 0;
}

int restore_disp_ZA() {				
	int grid_x, grid_y, grid_z, index;		
	double bias, dispmax, var_f;
	double** individual_displacements;
	//for ith particle, stores displacement x, displacement y, displacement z
	calloc2D_double(&individual_displacements, particle_no, 3);
	dispmax = 0.0;	
	
	/*	
	double* fsq_x = malloc((size_t)cells_displ[0]*cells_displ[1]*cells_displ[2] * sizeof(*fsq_x));
	double* fsq_y = malloc((size_t)cells_displ[0]*cells_displ[1]*cells_displ[2] * sizeof(*fsq_y));
	double* fsq_z = malloc((size_t)cells_displ[0]*cells_displ[1]*cells_displ[2] * sizeof(*fsq_z));
	*/
	
	//loop for x, y, z, performing inverse FFT for each, store displacements in individual_displacements
	for (int dimension = 0; dimension < 3; dimension++) {		
		//get the fourier displacements, inverse FFT
		displacements_ZA(dimension);
		fftw_execute(ZA_displacements_plan);
		
		/*
		//GIVES THE SAME ANSWER TO 2 DECIMAL POINTS AS AVERAGING OVER PARTICLES
		for (int i = 0; i < cells_displ[0]*cells_displ[1]*cells_displ[2]; i++) {
			if (dimension == 0) {
				fsq_x[i] = pow(Dplus*displacements[i][0], 2.0);
			} else if (dimension == 1) {
				fsq_y[i] = pow(Dplus*displacements[i][0], 2.0);
			} else {
				fsq_z[i] = pow(Dplus*displacements[i][0], 2.0);
			}		
		}
		*/
				
		//restore particles to "initial" positions, according to ZA		
		for (int i = 0; i < particle_no; i++) {		
			grid_x = (int) floor((particles[i][0] / volume_limits[0]) * (double) cells_displ[0]);
			grid_y = (int) floor((particles[i][1] / volume_limits[1]) * (double) cells_displ[1]);
			grid_z = (int) floor((particles[i][2] / volume_limits[2]) * (double) cells_displ[2]);
			index = arr_ind_displ(grid_x, grid_y, grid_z);	
			individual_displacements[i][dimension] = displacements[index][0];
			
			if (displacements[index][0] > dispmax) {
				dispmax = displacements[index][0];
			} 					
			//to redshift
			if (dimension == 2) {
				//particles[i][2] += fg*Dplus*displacements[index][0]; 			
			}				
			//velocities also scale in ZA according to http://arxiv.org/pdf/1408.1047v2.pdf
			//particles[i][dimension+3] += Dplus*a_primed*H_primed*fg_primed*displacements[index][0];
		}			
		ZA_clearMemory();
	}
	
	/*
	double exp_f_sq_cells = 0.0;
	double exp_sq_f_cells = 0.0;
	double f_sq_cells;
	double var_f_cells;
	for (int i = 0; i < cells_displ[0]*cells_displ[1]*cells_displ[2]; i++) {
		f_sq_cells = fsq_x[i] + fsq_y[i] + fsq_z[i];
		exp_f_sq_cells += f_sq_cells;
		exp_sq_f_cells += sqrt(f_sq_cells);	
	}	
	exp_f_sq_cells /= cells_displ[0]*cells_displ[1]*cells_displ[2];
	exp_sq_f_cells /= cells_displ[0]*cells_displ[1]*cells_displ[2];
	var_f_cells = exp_f_sq_cells - pow(exp_sq_f_cells, 2.0);
	printf("CELLS sqrt var f: %lf \n", sqrt(var_f_cells));
	*/
	
	
	//now getting variance of |f|
	double exp_f_sq = 0.0; //E(|f|^2)
	double exp_sq_f = 0.0; //E(|f|)^2 
	double f_sq;
	for (int i = 0; i < particle_no; i++) {			
		f_sq = pow(individual_displacements[i][0], 2.0) + pow(individual_displacements[i][1], 2.0) + pow(individual_displacements[i][2], 2.0);
		exp_f_sq += f_sq;
		exp_sq_f += sqrt(f_sq);	
	}
	
	exp_f_sq /= particle_no;
	exp_sq_f /= particle_no;
	exp_sq_f = pow(exp_sq_f, 2.0);
	var_f = exp_f_sq - exp_sq_f;
	printf("PARTICLES sqrt var f: %lf \n", sqrt(var_f));
	
	//final loop to apply displacements
	for (int i = 0; i < particle_no; i++) {
		//particles[i][0] += (sigma_exp_smoothed/sqrt(var_f))*individual_displacements[i][0];
		//particles[i][1] += (sigma_exp_smoothed/sqrt(var_f))*individual_displacements[i][1];
		//particles[i][2] += (sigma_exp_smoothed/sqrt(var_f))*individual_displacements[i][2];
		
		particles[i][0] += individual_displacements[i][0];
		particles[i][1] += individual_displacements[i][1];
		particles[i][2] += individual_displacements[i][2];
		
		//bias = b_v(m_to_v(particles[i][6]))/b_eff;
		//particles[i][dimension] += bias*Dplus*displacements[index][0];
		PBC(&particles[i][0], &particles[i][1], &particles[i][2]);		
	}

	
	if (dispmax > volume_limits[0] || dispmax > volume_limits[1] || dispmax > volume_limits[2]) {
		printf("ZA displacements are too big!!\n");
		return 1;
	} else {
		printf("max displacement here: %lf \n", dispmax);
	}
	free2D_double(&individual_displacements, 3);	
	return 0;	
}

int Pk_out(char outFile[]) {
	char out_temp[100];
	char out_temp_folded[100];
	sprintf(out_temp, "%s%s", outFile, "temp");
	sprintf(out_temp_folded, "%s%s", outFile, "temp_folded");
	//halos_run(out_temp, 1.0);
	//halos_run(out_temp_folded, 4.0);
	//combine_folded(out_temp, out_temp_folded, outFile, 0.019, 1.2, 0.45, true);	
	remove(out_temp);
	remove(out_temp_folded);
	return 0;
}

int get_Pk(char sim_path[], SplineInfo *spline) {
	char Pk_outFile[] = "/home/jonas/Testing_GR/code/Jonas/output/rescaling/target_power_temp.dat";	
	read_full_xyz(sim_path, 1, mocks_format);		
	Pk_out(Pk_outFile);
	*spline = input_spline_file(Pk_outFile, spline_format_Pk, false);
	remove(Pk_outFile);			
	return 0;
}

int halo_mass_func() {
	double dm, this_mass;
	int bins, index;
	bins = 50;
	char temp_out[100];
	//sprintf(temp_out, "%s%s", tempdir, "hmf_temp.dat");
	const char format[] = "%le \t %le \n";
	double* n_m = calloc((size_t)bins, sizeof(*n_m));
	double* hmf_m = calloc((size_t)(bins-1), sizeof(*hmf_m));
	
	dm = (Mmax - Mmin)/(double)bins;
	for (int i = 0; i < particle_no; i++) {
		this_mass = particles[i][6];
		index = (int) floor((this_mass - Mmin) / dm);
		if (index == bins) index = bins - 1;		
		n_m[index] += 1.0;
	}
	
	//n = N/V
	for (int i = 0; i < bins; i++) {
		n_m[i] /= volume;
	}
	
	FILE *out = fopen(temp_out, "w");
	for (int i = 0; i < bins-1; i++) {
		hmf_m[i] = ((n_m[i+1] - n_m[i])/dm);
		//take midpoint of bin
		this_mass = Mmin + ((double)i + 0.5)*dm;
		fprintf(out, format, this_mass, hmf_m[i]); 
	}
	
	fclose(out);
	//HMF_spline = input_spline_file(temp_out, format, false);	
	//remove(temp_out);
	return 0;
}

//GR -> F5
//scale box -> ZA -> HOD?
int rescale() {	
	//bool with_bias = false;
	//bool redshift_space = false;	
	//char GR_mono[] = "/home/jonas/Testing_GR/code/Jonas/output/rescaling/GR_monopole.dat";
	//char F5_mono[] = "/home/jonas/Testing_GR/code/Jonas/output/rescaling/F5_monopole.dat";
	
	char Pk_current[] = "/home/jonas/Testing_GR/code/Jonas/linear_matter_pk_sig8_0.593_z_0.75.dat";
	char Pk_target[] = "/home/jonas/Testing_GR/code/Jonas/linear_matter_pk_0.76.dat";
	
	char out_final[] = "/home/jonas/Testing_GR/code/Jonas/output/rescaling/GR_rescaled_test.dat";
	
	char out_model[] = "/home/jonas/Testing_GR/code/Jonas/output/rescaling/model_current.dat";
	char out_model_rescaled[] = "/home/jonas/Testing_GR/code/Jonas/output/rescaling/model_rescaled.dat";
	
	
	//store Pk spline for starting sim
	char sim_current[100];	
	sprintf(sim_current, "%s%s", mock_directory, "AllHalo_GR_1500Mpc_a0.7_1.txt");
	SplineInfo Pk_spline_temp_current;
	//get_Pk(sim_current, &Pk_spline_temp);
	
	
	//s_test = 0.8;
	//Dplus_test = 0.25;
	testmode = false;
	
	
	//Pk spline for target, model
	//SplineInfo Pk_spline_temp_target = input_spline_file(Pk_target, Pk_model_format, false);
	//Pk_spline_current = input_spline_file(Pk_current, Pk_model_format, false);
	//Pk_spline_target = input_spline_file(Pk_target, Pk_model_format, false);			
	
	
	//current simulation data input
	//read_full_xyz(sim_current, 1, mocks_format);
	//min and max masses, MUST COME BEFORE CUMSUMS
	//mass_min_max();
	//printf("mmin: %le,  mmax: %le\n", Mmin,  Mmax);
	//bias_plot();	

	//calc_b_eff();
	
	//get halo mass function
	//halo_mass_func();
	
	
	

	z_current = 0.75;
	z_target = 0.76;
	
	//find optimal s, z	
	//R1_primed = m_to_R(Mmin, z_current);
	//R2_primed = m_to_R(Mmax, z_current);
	//printf("R1': %le, R2': %le \n", R1_primed, R2_primed);
	//R1_primed = 0.1;
	//R2_primed = 10.0;
	//given z in current sim, vary the target z, using initial guess z_primed = z

	vary_z_current = true;
	//SET ZMIN, ZMAX, KMIN, KMAX SOMEWHERE ONCE AND FOR ALL
	//global variables
	k_bins = 300;
	z_bins = 200;
	Dplus_bins = 1000;
	int R_bins = 500;
	k_min = 0.001;
	k_max = 5.0;
	z_min = 0.0;
	z_max = 10.0;
	current_gravity = GR;
	target_gravity = GR;
	
	rescaling_k_bin_info = prep_bins(k_min, k_max, k_bins, true);	//logbin
	rescaling_z_bin_info = prep_bins(z_min, z_max, z_bins, false);
	rescaling_R_bin_info = prep_bins(0.6*R1_primed, 1.5*R2_primed, R_bins, true); //logbin
	
	//for current simulation:
	SplineInfo Pk_current_final;
	SplineInfo* Pk_splines_zBins = malloc((size_t)z_bins * sizeof(*Pk_splines_zBins));
	SplineInfo_Extended* Pk_splines = malloc((size_t)z_bins * sizeof(*Pk_splines));
	if (vary_z_current) {
		//prep_Pk_splines(current_gravity, &Pk_spline_temp_current, true, 0, NULL, &Pk_splines_zBins);
		for (int i = 0; i < z_bins; i++) {
			extend_spline_nonmodel(&Pk_splines_zBins[i], &Pk_splines[i]);
		}
	} else {
		//prep_Pk_splines(current_gravity, &Pk_spline_temp_current, false, z_current, &Pk_current_final, NULL);
		extend_spline_nonmodel(&Pk_current_final, &Pk_spline_current);
	}
	
	
	//for target simulation:
	

	//EdS
	//omega_m = 1.0;
	//omega_v = 0.0;
	//z_to_t_workspace = gsl_integration_workspace_alloc(1000);
	//prep_redshift_spline(0.0, 10.0, 1000);
	//SplineInfo tmp = Dplus_spline(GR, 0.0);
	
	
	//printf("Dplus: %le at redshift 1.6\n", splint_generic(redshift_spline_reversed, 0.005));
	//for (int i = 0; i < Dplus_k_bins; i++) {
		//printf("k bin %d, z: %lf, Dplus: %lf\n", i, 1.1366, splint_generic(Dplus_splines[i], 1.1366));
	//}
	//for (int i = 0; i < OLV_z_bins; i++) {
		//printf("z bin %d, k: 0.01, Pk: %le, Pk model: %le \nk: 1.0, Pk: %le, Pk model: %le \n", i, splint_generic(Pk_splines_zBins[i], 0.01), splint_Pk_model(Pk_splines_zBins_model[i], 0.01), splint_generic(Pk_splines_zBins[i], 1.0), splint_Pk_model(Pk_splines_zBins_model[i], 1.0));
	//}
	//printf("const z. k: 0.01, Pk: %le, Pk model: %le \nk: 1.0, Pk: %le, Pk model: %le\n", splint_generic_VPk(Pk_const_z_spline, 0.0005), splint_Pk_model(Pk_const_z_spline_model, 0.0005), splint_generic_VPk(Pk_const_z_spline, 1.0), splint_Pk_model(Pk_const_z_spline_model, 1.0));
	
	//printf("prepping variance splines\n");

	//prep_variance_splines(R_bins);
	//printf("starting multimin\n\n");
	//dsq_multimin();
	//printf("generating sigma plot data\n");
	//generate_sigma_plots();

	//MODIFY PK_CURRENT, PK_TARGET TO ACCOUNT FOR S, Z RESCALE, FOR USE IN ZA
	
	
	
	
	

	/*	
	volume_limits[0] *= s;
	volume_limits[1] *= s;
	volume_limits[2] *= s;
	volume = volume_limits[0]*volume_limits[1]*volume_limits[2];
	
	for (int i = 0; i < particle_no; i++) {
		PBC(&particles[i][0], &particles[i][1], &particles[i][2]);
	}
	*/
	
	
	/*
	populateGridCIC();
	gridIntoOverdensity();
	overdensity_pk(false);
	pk_to_file_logbin(out_final);
*/
	
	/*
	for (int i = 0; i < particle_no; i++) {
		PBC(&particles[i][0], &particles[i][1], &particles[i][2]);
	}	
	*/
	
	//scale masses, velocities with s
	//does this come before or after HOD??
	//NEED SOME PARAMETERS HERE!!
	//scale_masses();
	//scale_velocities();
	
	
	//Pk_out(out_final);
	
	/*
	
	//move around haloes
	//displacements and velocities! Don't necessarily need to free this yet		
	fftw_allocateMemory_ZA();
	restore_disp_ZA();
	fftw_free(displacements);
	fftw_free(displacements_fourier);
	*/
	
	
	return 0;		
}

/*
int HOD_catalogue() {
	double cos_theta_sat, theta_sat, phi_sat, r_sat, unif_rnd, mean_r_sat;
	double erf_central, mean_N_satellite, mean_N_central, halo_mass;
	double x_sat, y_sat, z_sat, dm_nfw;
	//velocity assignment for satellites variables
	double result_I, error, f_cdm, sigma_sq_sat, sigma_sq_sat_tcdm, vx, vy, vz, sigma;
	int has_central, satellite_number, total_sats, total_centrals, nfw_bins, nfw_index;
	
	//prepare NFW profile cumulative sums
	nfw_bins = 100;
	nfw_cumsums_reversed = malloc((size_t)nfw_bins * sizeof(*nfw_cumsums_reversed));
	nfw_cumsums_catalogue(100, "/home/jonas/Testing_GR/code/Jonas/output/nfw_cumsums/");
	dm_nfw = (Mmax - Mmin)/(double)nfw_bins;
	
	//find out how many centrals and allocate them
	total_centrals = 0;
	for (int i = 0; i < particle_no; i++) {
		halo_mass = particles[i][6];
		erf_central = gsl_sf_erf((log10(halo_mass) - log_Mmin)/sigma_logm);
		mean_N_central = 0.5 * (1.0 + erf_central);		
		has_central = (int)gsl_ran_bernoulli(gsl_ran_r, (unsigned int)mean_N_central);
		if (has_central == 1) {
			particles[total_centrals][0] = particles[i][0];
			particles[total_centrals][1] = particles[i][1];
			particles[total_centrals][2] = particles[i][2];
			particles[total_centrals][3] = particles[i][3];
			particles[total_centrals][4] = particles[i][4];
			particles[total_centrals][5] = particles[i][5];
			particles[total_centrals][6] = particles[i][6];			
			total_centrals += 1;
		}
	}
	
	//either this, or manually zero the rest of memory
	central_no = total_centrals;
	particles = realloc(particles, (size_t)central_no * sizeof(*particles));
	
	//find out how many satellites to allocate for each HALO OR CENTRAL?????
	total_sats = 0;
	int* no_satellites = calloc((size_t)central_no, sizeof(*no_satellites));	
	for (int i = 0; i < central_no; i++) { //<----- OR LOOP THROUGH HALOES
		halo_mass = particles[i][6];
		mean_N_satellite = pow((halo_mass - M0)/M1, alpha);
		satellite_number = (int)gsl_ran_poisson(gsl_ran_r, (unsigned int)mean_N_satellite);			
		no_satellites[i] = satellite_number;
		total_sats += satellite_number;		
	}
	
	//realloc to add memory for satellites
	satellite_no = total_sats;	
	particles = realloc(particles, (size_t)(central_no + satellite_no) * sizeof(*particles));
	//calloc zeros out the values of satellites
	for (int i = central_no; i < central_no + satellite_no; i++) {
		particles[i] = (double*) calloc(7, sizeof(**particles));
	}		
	
	//prepare integral for velocity assignment of satellites		
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);	
	gsl_function F;
  	F.function = &I_integral_HOD;
  	F.params = 0;   
  	//end velocity assignment stuff			

	//add satellites, this counter is to place satellites in the now expanded particles array
	int satellite_index = central_no;
	for (int i = 0; i < central_no; i++) {
		halo_mass = particles[i][6];
		nfw_index = (int) floor((halo_mass - Mmin)/dm_nfw);
		//each central i will have no_satellites[i] as found before
		for (int j = 0; j < no_satellites[i]; j++) {
			//cos theta is uniform random in [-1,1]
			cos_theta_sat = randomDouble(-1.0, 1.0);
			theta_sat = acos(cos_theta_sat);
			
			//phi is uniform random in [0,2pi]
			phi_sat = randomDouble(0.0, 2.0*M_PI);
			
			//selects a cumulative sum in range [cumsum_min, cumsum_max]
			//splint (inverse) then relates that cumsum to some r value between 0 and r_virial
			unif_rnd = randomDouble(nfw_cumsums_reversed[nfw_index].xmin, nfw_cumsums_reversed[nfw_index].xmax);

			r_sat = splint_generic(nfw_cumsums_reversed[nfw_index], unif_rnd);
			
			//track the mean radial displacement
			mean_r_sat += r_sat;			
		
			//actual position is with respect to the i-th central
			x_sat = r_sat * sin(theta_sat) * cos(phi_sat) + particles[i][0];
			y_sat = r_sat * sin(theta_sat) * sin(phi_sat) + particles[i][1];
			z_sat = r_sat * cos(theta_sat) + particles[i][2];			

			PBC(&x_sat, &y_sat, &z_sat);
		
			//place the satellite in particle array		
			particles[satellite_index][0] = x_sat;
			particles[satellite_index][1] = y_sat;
			particles[satellite_index][2] = z_sat;	
			
			//finally, assign velocities			
			gsl_integration_qags(&F, r_sat/rs_global, 1000.0, 0, 1e-5, 1000, w, &result_I, &error); 
	
			f_cdm = log(1.0 + cdm_global) - cdm_global/(1.0 + cdm_global);	
			//sigma squared from the paper with the rest of HOD recipe
			sigma_sq_sat = (G*halo_mass*cdm_global*cdm_global*r_sat/(R_virial_global*R_virial_global*f_cdm))*pow(1.0 + cdm_global*r_sat/R_virial_global, 2.0)*result_I;
			
			//sigma squared from tcdm paper, DIFFERENT
			sigma_sq_sat_tcdm = (G*halo_mass/(3.0*R_virial_global))*cdm_global*(1.0 - 1.0/(1.0+pow(cdm_global,2.0)) - 2.0*log(1.0 + cdm_global)/(1.0 + cdm_global))/(2.0*pow(log(1.0 + cdm_global) - cdm_global/(1.0 + cdm_global),2.0));
			
			double ratio = sigma_sq_sat_tcdm/sigma_sq_sat;
			if (ratio > 1.1 || ratio < 0.9) {
				//printf("ratio: %lf\n", ratio);
			}
	
			//this is the relative velocity to the halo, so add together with the host halo velocity
			sigma = sqrt(sigma_sq_sat);
			vx = gsl_ran_gaussian(gsl_ran_r, sigma_sq_sat);
			vy = gsl_ran_gaussian(gsl_ran_r, sigma_sq_sat);
			vz = gsl_ran_gaussian(gsl_ran_r, sigma_sq_sat);	
			
			particles[satellite_index][3] = particles[i][3] + vx;
			particles[satellite_index][4] = particles[i][4] + vy;
			particles[satellite_index][5] = particles[i][5] + vz;		
			
			satellite_index += 1;
		}		
	}
	
	return 0;
}
*/
