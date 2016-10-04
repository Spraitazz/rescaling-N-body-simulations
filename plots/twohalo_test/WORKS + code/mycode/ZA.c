int ZA_displacements(int dim, bool CIC) {
	double kx, ky, kz, k, k_sq, coef_real, coef_im;
	double phase, amplitude, Pk, WindowFunc;
    int index, Index, negkIndex;
 	double fundmodes[3], nyquist[3], cellSizes[3];  
 	
 	cellSizes[0] = volume_limits[0] / (double) cells_displ[0];
	cellSizes[1] = volume_limits[1] / (double) cells_displ[1];
	cellSizes[2] = volume_limits[2] / (double) cells_displ[2];

    fundmodes[0] = 2.0*pi/volume_limits[0];
    fundmodes[1] = 2.0*pi/volume_limits[1];
    fundmodes[2] = 2.0*pi/volume_limits[2];
    
    nyquist[0] = pi / cellSizes[0];
    nyquist[1] = pi / cellSizes[1];
    nyquist[2] = pi / cellSizes[2];
 
	for (int aa = 0; aa < cells_displ[0]; aa++) {
		for (int bb = 0; bb < cells_displ[1]; bb++) {
			for (int cc = 0; cc < cells_displ[2]; cc++) {
				kx = (double)aa * fundmodes[0];
				ky = (double)bb * fundmodes[1];
				kz = (double)cc * fundmodes[2];  			
			
				//convention - up to N/2 positive frequencies, N/2 is fc, then -fc+ffund, ...
				if (aa > cells_displ[0]/2) kx -= ((double) cells_displ[0]) * fundmodes[0]; //or -= 2*kx_nyquist
				if (bb > cells_displ[1]/2) ky -= ((double) cells_displ[1]) * fundmodes[1];
				if (cc > cells_displ[2]/2) kz -= ((double) cells_displ[2]) * fundmodes[2];
			
				index = arr_ind_displ(aa, bb, cc);					
				k_sq = kx*kx + ky*ky + kz*kz;				
				k = sqrt(k_sq);	
				WindowFunc = 1.0;				
			
				if (k_sq > 0.0) {					
								
					if (CIC) {					
						if(kx != 0.0) WindowFunc *= sin(pi*kx*0.5/nyquist[0])/(pi*kx*0.5/nyquist[0]);
						if(ky != 0.0) WindowFunc *= sin(pi*ky*0.5/nyquist[1])/(pi*ky*0.5/nyquist[1]);
						if(kz != 0.0) WindowFunc *= sin(pi*kz*0.5/nyquist[2])/(pi*kz*0.5/nyquist[2]);						
					}			
								
					Pk = splint_generic_VPk(&Pk_current, k);
					//Pk = splint_Pk(&Pk_current_extended, k);
				
					phase = phases[index];
								
					amplitude = sqrt(Pk);
					//amplitude = sqrt(Pk) * exp(-k_sq*R_nl*R_nl/2.0);
				
					amplitude /= pow(WindowFunc, 2.0);
				
					coef_real = amplitude * cos(phase);
					coef_im = amplitude * sin(phase);	
												
					switch (dim) {					
						case 0:
							displacements_fourier[index][0] = kx * coef_im / k_sq;							
							displacements_fourier[index][1] = -kx * coef_real / k_sq;
							break;
						case 1:						
							displacements_fourier[index][0] = ky * coef_im / k_sq;
							displacements_fourier[index][1] = -ky * coef_real / k_sq;										
							break;
						case 2:						
							displacements_fourier[index][0] = kz * coef_im / k_sq;
							displacements_fourier[index][1] = -kz * coef_real / k_sq;							
							break;
					}	
				}						
			}
		}		
	}
	
	// be safe: intialise k=0 H_k to zero. 
	displacements_fourier[0][0] = 0.0;
	displacements_fourier[0][1] = 0.0;
	
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
            if(i!=0) Index += (cells_displ[0] - i);
            Index += (cells_displ[1] - j)*cells_displ[0];
            displacements_fourier[negkIndex][0] = displacements_fourier[Index][0];
            displacements_fourier[negkIndex][1] = -1.0*displacements_fourier[Index][1];
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

int apply_ZA(bool CIC) {		

	int index, grid_x, grid_y, grid_z;	
	double** individual_displacements;
	calloc2D_double(&individual_displacements, particle_no, 3);	

	//As doing each dimension separately, phases are allocated in advance. IMPORTANT!!
	phases = calloc((size_t)cells_displ[0]*cells_displ[1]*cells_displ[2], sizeof(*phases));
	for (int i = 0; i < cells_displ[0]*cells_displ[1]*cells_displ[2]; i++) {
		phases[i] = randomDouble(0.0, 2.0*pi);
	}
	
	//CIC stuff
	double xc, yc, zc, dx, dy, dz, tx, ty, tz, thisDisp;
	double cellSizes[3];
	cellSizes[0] = volume_limits[0] / (double) cells_displ[0];
	cellSizes[1] = volume_limits[1] / (double) cells_displ[1];
	cellSizes[2] = volume_limits[2] / (double) cells_displ[2];
	//end CIC stuff

	//store displacements for each particle and each dimension in a [][], apply later, check max for sanity
	double dispmax = 0.0;
	double dispmax_im = 0.0;
	for (int dimension = 0; dimension < 3; dimension++) {		
		// get the fourier displacements, inverse FFT
		ZA_displacements(dimension, CIC);		
		fftw_execute(ZA_displacements_plan);	
		
		//restore particles to "initial" positions, according to ZA		
		for (int i = 0; i < particle_no; i++) {		
			grid_x = (int) floor((particles[i][0] / volume_limits[0]) * (double) cells_displ[0]);
			grid_y = (int) floor((particles[i][1] / volume_limits[1]) * (double) cells_displ[1]);
			grid_z = (int) floor((particles[i][2] / volume_limits[2]) * (double) cells_displ[2]);
			index = arr_ind_displ(grid_x, grid_y, grid_z);			
			
			if (CIC) {			
				xc = cellSizes[0] * ((double) grid_x + 0.5);
				yc = cellSizes[1] * ((double) grid_y + 0.5);
				zc = cellSizes[2] * ((double) grid_z + 0.5);
		
				//must be from 0 to 1
				dx = fabs(particles[i][0] - xc) / cellSizes[0];
				dy = fabs(particles[i][1] - yc) / cellSizes[1];
				dz = fabs(particles[i][2] - zc) / cellSizes[2];	
		
				tx = 1.0 - dx;
				ty = 1.0 - dy;
				tz = 1.0 - dz;
		
				//ensure PBC
				int xplus, yplus, zplus;
				xplus = (grid_x + 1) % cells_displ[0];
				yplus = (grid_y + 1) % cells_displ[1];
				zplus = (grid_z + 1) % cells_displ[2];
			
				thisDisp = displacements[arr_ind(grid_x, grid_y, grid_z)][0]*tx*ty*tz + displacements[arr_ind(xplus, grid_y, grid_z)][0]*dx*ty*tz + displacements[arr_ind(xplus, yplus, grid_z)][0]*dx*dy*tz + displacements[arr_ind(xplus, grid_y, zplus)][0]*dx*ty*dz + displacements[arr_ind(xplus, yplus, zplus)][0]*dx*dy*dz + displacements[arr_ind(grid_x, yplus, grid_z)][0]*tx*dy*tz + displacements[arr_ind(grid_x, grid_y, zplus)][0]*tx*ty*dz;	
			
				individual_displacements[i][dimension] = thisDisp;	
			
			} else {
				individual_displacements[i][dimension] = displacements[index][0];
			}			
			    
			// calculate biggest and smallest displacement.     										
			if (displacements[index][0] > dispmax) {
				dispmax = displacements[index][0];
			} 	
			if (displacements[index][1] > dispmax_im) {
				dispmax_im = displacements[index][1];
			}
		}						
		ZA_clearMemory();
	}
	
	//now getting variance of |f|
	double exp_f_sq = 0.0; //E(|f|^2)
	double exp_sq_f = 0.0; //E(|f|)^2 
	double f_sq, var_f;
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
	
	//applying displacements!!
	for (int i = 0; i < particle_no; i++) {		
		particles[i][0] += individual_displacements[i][0];
		particles[i][1] += individual_displacements[i][1];
		particles[i][2] += individual_displacements[i][2];		
		if (redshift_space) {		
			particles[i][2] += fg * individual_displacements[i][2];
		}
		
		PBC(&particles[i][0], &particles[i][1], &particles[i][2]);		
	}
	
	printf("max disp: %lf, max disp im: %lf \n", dispmax, dispmax_im);
		
	free2D_double(&individual_displacements, 3);
	free(phases);
	return 0;
}

int ZA_test() {
	// model spline for Pk
	char Pk_in_model[100];
	sprintf(Pk_in_model, "%s/linear_matter_pk_sig8_0.593_z_0.75.dat", home_directory);
	Pk_current = input_spline_file(Pk_in_model, Pk_model_format, false);		
	
	//PARAMETERS
	set_params();
	redshift_space = false;
	z_current_global = 3.5;	
	
	//binning
	k_bins = 300;	
	Dplus_bins = 1000;
	k_min = 0.99*Pk_current.xmin;
	k_max = 1.01*Pk_current.xmax;	
	z_min = 0.0;
	z_max = 5.0;
	rescaling_k_bin_info = prep_bins(k_min, k_max, k_bins, true); //logbin	
	//rescaling_R_bin_info = prep_bins(0.001, 10.0, 300, true); //logbin
	
	//printf("pk b4: %lf \n", splint_generic(&Pk_current, 0.5));
	prep_Pk_constz(GR, z_current_global, &Pk_current, &rescaling_k_bin_info);
	//printf("pk after: %lf \n", splint_generic(&Pk_current, 0.5));
	
	/*
	//WHY USING REGRESSION EXTENSION WORKS, BUT WITHOUT IT DOESNT????
	prep_Pk_model_spline(&Pk_current, &Pk_current_extended);
	//extend_spline_nonmodel(&Pk_current, &Pk_spline_current);	

	char OLV_out[100];
	sprintf(OLV_out, "%s/output/HOD_rands/OLV.dat", home_directory);
	prep_variances(&rescaling_R_bin_info, OLV_out, &Pk_current_extended);
	variance_spline_current = input_spline_file(OLV_out, "%le \t %le \n", false);
	
	reverse_spline(&variance_spline_current, &variance_spline_current_reversed);	
	R_nl = splint_generic(&variance_spline_current_reversed, 1.0);	
	printf("Rnl is: %le\n", R_nl);	
	
	sigma_exp_smoothed = sqrt(expected_variance_smoothed(R_nl, &Pk_current_extended));
	printf("expected standard deviation of smoothed displacement: %lf \n", sigma_exp_smoothed);	
	
	*/		

	char out_model[100];
	sprintf(out_model, "%s/output/HOD_rands/model_pk.dat", home_directory);
	haloModel_out(out_model, 0.02, 1.15, 200);				
	
	char outFile[100];
	char tmpFile[100];
	char tmpFile2[100];
	sprintf(tmpFile, "%s/output/HOD_rands/unfolded.dat", home_directory);
	sprintf(tmpFile2, "%s/output/HOD_rands/folded.dat", home_directory);	

	//prep_oneHalo();	

	
	for (int i = 0; i < 20; i++) {
		sprintf(outFile, "%s/output/HOD_rands/test/measured_pk_%d.dat", home_directory, i);
		printf("run %d \n", i);
		
		//populateParticles_lattice();
		populateParticles();
		
		ZA_allocateMemory();
		apply_ZA(false);
		ZA_freeMemory();
		
		//halo_no = particle_no; //now particle number will change, stay different for measuring pk, then revert
		//oneHalo();
		
		overdensity_allocateMemory();
		haloes_measure_Pk(tmpFile, 1.0, false);	
		haloes_measure_Pk(tmpFile2, 4.0, false);
		combine_folded(tmpFile, tmpFile2, pk_format_folding, outFile, 0.0, 1.1, 0.3, true);
	
		
		overdensity_freeMemory();
		
		clearParticles();
		remove(tmpFile);
		remove(tmpFile2);
		//particle_no = halo_no;
		//particles = realloc(particles, particle_no * sizeof(*particles));			
		
	}
	
	//gsl_integration_workspace_free(onehalo_workspace);
	

	return 0;
}

int prep_oneHalo() {

	double halo_mass = 5e13;
	set_halo_params(halo_mass);
	char cumsums_out[100];
	sprintf(cumsums_out, "%s/output/HOD_rands/cumsum.dat", home_directory);
	prep_nfw_cumsums(cumsums_out);
	nfw_spline_reverse = input_spline_file(cumsums_out, "%le \t %le \n", true);
	onehalo_workspace = gsl_integration_workspace_alloc(1000);
	return 0;
}

int oneHalo() {

	double cos_theta_sat, theta_sat, phi_sat, r_sat, unif_rnd, mean_r_sat;
	double x_sat, y_sat, z_sat;
	//velocity assignment for satellites variables
	double result_I, error, rho_nfw, f_cdm, sigma_sq_sat, sigma_sq_sat_tcdm, vx, vy, vz, sigma;
	int satellite_number, has_central, particle_no_init, total_sats;
	
	double halo_mass = 5e13;		
	
	//this should be in principle different
	central_no = particle_no;
	
	int mean_N_sats = pow((halo_mass - M0)/M1, alpha);
	mean_r_sat = 0.0;
	total_sats = 0;
	
	int* no_satellites = calloc(particle_no, sizeof(*no_satellites));
		
	for (int i = 0; i < central_no; i++) {
		satellite_number = gsl_ran_poisson(gsl_ran_r, mean_N_sats);
		no_satellites[i] = satellite_number;
		total_sats += satellite_number;
	}
	printf("adding %d satellites\n", total_sats);
	
	satellite_no = total_sats;
	particle_no += total_sats;
	
	//realloc to add memory for satellites
	particles = realloc(particles, particle_no * sizeof(*particles));
	//calloc zeros out the values of satellites
	for (int i = central_no; i < particle_no; i++) {
		particles[i] = (double*) calloc(7, sizeof(**particles));
	}		
	
	//prepare integral for velocity assignment of satellites
	gsl_function F;
  	F.function = &I_integral_HOD;
  	F.params = 0;   
  	//end velocity assignment stuff	
				
	//add satellites, this counter is to place satellites in the now expanded particles array
	int satellite_index = particle_no_init;
	for (int i = 0; i < central_no; i++) {
		//each central i will have no_satellites[i] as found before
		for (int j = 0; j < no_satellites[i]; j++) {

			cos_theta_sat = randomDouble(-1.0, 1.0);
			theta_sat = acos(cos_theta_sat);
			phi_sat = randomDouble(0.0, 2.0*pi);
			
			//selects a cumulative sum in range [cumsum_min, cumsum_max]
			//splint (inverse) then relates that cumsum to some r value between 0 and r_virial
			unif_rnd = randomDouble(nfw_spline_reverse.xmin, nfw_spline_reverse.xmax);

			r_sat = splint_generic(&nfw_spline_reverse, unif_rnd);				
			
			//track the mean radial displacement
			mean_r_sat += r_sat;		
		
			//actual position is with respect to the i-th central
			x_sat = r_sat * sin(theta_sat) * cos(phi_sat) + particles[i][0];
			y_sat = r_sat * sin(theta_sat) * sin(phi_sat) + particles[i][1];
			z_sat = r_sat * cos(theta_sat) + particles[i][2];		
		
			//PBC
			PBC(&x_sat, &y_sat, &z_sat);
		
			//place the satellite in particle array		
			particles[satellite_index][0] = x_sat;
			particles[satellite_index][1] = y_sat;
			particles[satellite_index][2] = z_sat;	
			
			//finally, assign velocities	
			
			gsl_integration_qags(&F, r_sat/rs_global, 1000.0, 0, 1e-5, 1000, onehalo_workspace, &result_I, &error); 
	
			f_cdm = log(1.0 + cdm_global) - cdm_global/(1.0 + cdm_global);	
			sigma_sq_sat = (G*halo_mass*cdm_global*cdm_global*r_sat/(R_virial_global*R_virial_global*f_cdm))*pow(1.0 + cdm_global*r_sat/R_virial_global, 2.0)*result_I;
			
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
	
	//sanity check?
	mean_r_sat /= (double) total_sats;
	printf("mean r of sats: %lf \n", mean_r_sat); 
	
	free(no_satellites);
	return 0;
}
/*
int ZA_test_Mikes() {
	// model spline for Pk
	//MY
	char Pk_in_model[100];
	sprintf(Pk_in_model, "%s/linear_matter_pk_sig8_0.593_z_0.75.dat", home_directory);
	Pk_current = input_spline_file(Pk_in_model, Pk_model_format, false);
	
	k_bins = 300;
	k_min = 0.001;
	k_max = 11.35;
	rescaling_k_bin_info = prep_bins(k_min, k_max, k_bins, true);	//logbin
	prep_Pk_constz(GR, 0.0, &Pk_current, &rescaling_k_bin_info);
	
	rescaling_R_bin_info = prep_bins(0.001, 30.0, 500, true); //logbin
	prep_variance_splines(&Pk_spline_current, &variance_spline_current, NULL, NULL, false, z_test);
	reverse_spline(&variance_spline_current, &variance_spline_current_reversed);
	
	SplineInfo_Extended Pk_current_ext;
	prep_Pk_model_spline(&Pk_current, &Pk_current_ext);
	

	// MY
	char out_model[100];
	sprintf(out_model, "%s/output/HOD_rands/model_pk.dat", home_directory);
	Dplus_test = 0.4;	
	haloModel_out(out_model, 0.015, 0.65, 400);	
	
	n0 = n1 = n2 = cells[0];
	
	
	
	for (int i = 0; i < 100; i++) {
	
		prep_DisplacementCalc();
		printf("displ calc prepared \n");

		HaloCatalogue_NFWprofile(particle_no, -1); //total_gal? not necessary for this?	
		printf("halo catalogue ready \n");	
		
		sprintf(outFile_Pk_Mikes, "%s/output/HOD_rands/test/measured_pk_%d.dat", home_directory, i);
				
		//MY CODE (NGP)
		populateGrid();
		
		// initializes memory for overdensity
		prep_grid();
		
		//MY CODE (overdensity)
		gridIntoOverdensity();		

		prep_fftw();
		printf("starting pk calc\n");
		PkCalc();
		
		//MY
		clearGrid();
		clearParticles();
		
		
		free_pkRegression();		
	
	}

	return 0;
}

double norm, res;
double Sigma8 = 0.9;
double ShapeGamma = 0.21;   //% only needed for Efstathiou power spectrum
double UnitLength_in_cm = 3.085678e21;  // % defines length unit of output (in cm/h)
double AA, BB, CC, nu;

int init_Efstathiou_pk(){
	double norm = 1.0;

	norm = 	Sigma8*Sigma8/TopHatSigma2(8.0, norm);

	// call like this from now on. 
	// PowerSpec_Efstathiou(k, norm);

	return 0;
}

double PowerSpec_Efstathiou(double k, double norm){


	//gamma -> efstathiou.
	AA = 6.4 / ShapeGamma * (3.085678e24 / UnitLength_in_cm);
	BB = 3.0 / ShapeGamma * (3.085678e24 / UnitLength_in_cm);
	CC = 1.7 / ShapeGamma * (3.085678e24 / UnitLength_in_cm);
	nu = 1.13;

	return norm * k / pow(1. + pow(AA * k + pow(BB * k, 1.5) + CC * CC * k * k, nu), 2 / nu);
}


double TopHatSigma2(double R, norm){
	double result, abserr;
	gsl_integration_workspace *workspace;
	gsl_function F;

	workspace = gsl_integration_workspace_alloc(1000);

	F.function = &sigma2_int;

	r_tophat = R;

	gsl_integration_qag(&F, 0, 500.0 * 1 / R, 0, 1.0e-8, 1000, GSL_INTEG_GAUSS41, workspace, &result, &abserr);

	gsl_integration_workspace_free(workspace);

	return result;

	// note: 500/R is here chosen as (effectively) infinity integration boundary 
}


double sigma2_int(double k, void *param){
  double kr, kr3, kr2, w, x;
 
  norm = param[0]; //double, should be one initially. 

  kr = r_tophat * k;
  kr2 = kr * kr;
  kr3 = kr2 * kr;

  if(kr < 1e-8)
    return 0;

  w = 3 * (sin(kr) / kr3 - cos(kr) / kr2);
  x = 4 * PI * k * k * w * w * PowerSpec_Efstathiou(k, norm);

  return x;
}

*/
