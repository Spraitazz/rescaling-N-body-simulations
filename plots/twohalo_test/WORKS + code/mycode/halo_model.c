//selects halos with mass above m_min, outputs to new file.
//NOT CURRENTLY USED
/*
int select_halos_mmin(char outFile[]) {
	double log_Mmin, Mmin;
	log_Mmin = 14.06;
	Mmin = pow(10., log_Mmin);	
	//catalogue should be read in
	FILE *f = fopen(outFile, "w");	
	for (int i = 0; i < particle_no; i++) {
		if (particles[i][6] > Mmin) {
			fprintf(f, mocks_format, particles[i][0], particles[i][1], particles[i][2],
			particles[i][3], particles[i][4], particles[i][5], particles[i][6]);
			fprintf(f, "\n");
		}
	}	
	fclose(f);
	return 0;
}
*/

//sets R_nl for ZA smoothing
//VARIANCE SHOULD BE SPLINED, STORED
int set_ZA_params() {
	//R nonlinear for ZA smoothing function
	//taken to be the radius at which the overdensity linear variance is 1
	//should use sigma_nl squared, but 1 in this case!!	
	//sigma_nl = 1.0;	
	//R_nl = splint_generic(&variance_spline_current_reversed, sigma_nl);	
	printf("Rnl is: %lf\n", R_nl);	
	return 0;
}

//sets the global variables cdm_global and R_virial_global for use in all functions given a halo mass
//VARIANCE SHOULD BE SPLINED, STORED
int set_halo_params(double halo_mass) {
	//getting Mstar for cdm, using the overdensity linear variance obtained from the model power spectrum	
	//printf("here?\n");
	R_mstar = splint_generic(&variance_spline_current_reversed, delta_c * delta_c);

	//m* is the nonlinear mass scale at z = 0 defined such that sigma^2 (m*) = delta_c ^ 2
	//inverse spline returns R when sigma^2 has this value, then get m*
	Mstar = R_to_m(R_mstar, z_current_global);
	//printf("r mstar %lf, mstar %lf\n", R_mstar, Mstar);	
	cdm_global = (c0 / (1.0 + z_current_global)) * pow(halo_mass / Mstar, beta);		
	R_virial_global = pow(3.0*halo_mass/(4.0*M_PI*rho_bar_z(z_current_global)*del_nl), 1.0/3.0);	
	rs_global = R_virial_global/cdm_global;
	printf("m*: %le, cdm: %le, Rvir: %le, rs: %le\n", Mstar, cdm_global, R_virial_global, rs_global);	
	return 0;
}

//cumulative sum in rho NFW up to r, given a halo of mass m
double cumsum_nfw(double r) {
	double rs = R_virial_global / cdm_global;
	if (r == R_virial_global) {
		return 1.0;
	} else if (r == 0.0) {
		return 0.0;
	} else {
		//getting rid of outside constants, additional term comes from integral
		return (log((rs + r)/rs) + rs/(rs + r) - 1.0) * pow(log(1.0 + cdm_global) - cdm_global/(1.0 + cdm_global),-1.0);
	}
}

//NOT CURRENTLY USING THIS
/*
int prep_nfw_cumsums_Mikes(double mass, char outFile[]) {	
	double rs, r, q;
	rs = R_virial_global/cdm_global;
	
	FILE *f = fopen(outFile, "w");
	for(int j=0; j<1000; j++){  
		r = ((double) j/1000.0)*(R_virial_global/rs); // r_nfwinversion[j] is in units of r_s

		// Cumulative mass up to r
		q = (log(1.0 + r) - r*pow(1.0 + r, -1.0))*pow(log(1.0 + cdm_global) - cdm_global/(1.0 + cdm_global), -1.0);
		fprintf(f, "%le \t %le \n", r, q);
	}
	fclose(f);
	return 0;
}
*/

//outputs integral of rho_nfw * dV from 0 to r for multiple r going from 0 to r_virial
int prep_nfw_cumsums(char outFile[]) {	
	double rstart, rstop, dr, r, cumsum;
	int bins = 1000;
	
	rstart = 0.0;
	rstop = R_virial_global;
	dr = (rstop - rstart) / (double) (bins-1);
	FILE *f = fopen(outFile, "w");
	
	for (int i = 0; i < bins; i++) {
		r = rstart + (double)i * dr;
		cumsum = cumsum_nfw(r);
		fprintf(f, "%le \t %le \n", r, cumsum);
	}
	fclose(f);
	return 0;
}

//prepares nfw cumulative sums for a given catalogue, places into array of SplineInfo objects
//MUST KNOW Mmin and Mmax beforehand, stored globally
int nfw_cumsums_catalogue(int bins, const char root_dir[]) {
	double dm, this_mass;
	char out_file[100];	
	dm = (Mmax - Mmin)/(double)(bins-1);
	for (int i = 0; i < bins; i++) {
		this_mass = Mmin + (double)i * dm;
		//printf("this mass is: %le\n", this_mass);
		sprintf(out_file, "%scumsums_mass_%d.dat", root_dir, i);
		//set the r_vir and the cdm parameters!!!
		set_halo_params(this_mass);		
		prep_nfw_cumsums(out_file);	
		//nfw_cumsums_reversed[i] = input_spline(out_file, "%le \t %le \n", true);
		remove(out_file);
	}
	return 0;
}

//for model power spectrum, COPIED FROM MIKE
double ukm_nfw_profile(double k){
    // For u(k|m) = (1/M)\int_vol \rho(x|m) exp^{-ikx}
    // Spherically symmetric profile, truncated at the virial radius. 
    // u(k|m) = \int_0^{r_{vir}} dr 4pi r^2 [sin(kr)/kr] rho(r|m)/m
    
    double Interim, rs; 
    
    rs = R_virial_global/cdm_global; 
    
    Interim  = sin(k*rs)*(gsl_sf_Si((1.0 + cdm_global)*k*rs) - gsl_sf_Si(k*rs)) - sin(cdm_global*k*rs)*pow((1.0 + cdm_global)*k*rs, -1.0) + cos(k*rs)*(gsl_sf_Ci((1.0 + cdm_global)*k*rs) - gsl_sf_Ci(k*rs));

    // for the NFW profile, defined by rhos, rs and c, mass is not independent. 
    Interim /= log(1.0 + cdm_global) - cdm_global/(1.0 + cdm_global);
    return Interim;
}



double muOrderZero(double ks){
    // limits were established with signa = 2.*3./sqrt(2.), limits should scale propto velDispersion/(2.*3./sqrt(2.))
    if(ks > 0.0914289)  return 0.5*pow(pi, 0.5)*gsl_sf_erf(ks)/ks;

    else{
        return 1.0 - 0.333333*pow(ks, 2.) + 0.1*pow(ks, 4.) - 0.0238095*pow(ks, 6.);
    }
}


double muOrderTwo(double ks){
    if(ks > 0.15727469)  return pow(4.*ks*ks*ks, -1.)*(pow(pi, 0.5)*gsl_sf_erf(ks) - 2.*ks*exp(-1.*ks*ks));
    
    else{
        return 1./3. - ks*ks/5. + pow(ks, 4.)/14.;
    }
}


double muOrderFour(double ks){
    if(ks>0.2418305)  return pow(8.*pow(ks, 5.), -1.)*(3.*pow(pi, 0.5)*gsl_sf_erf(ks) -6.*ks*exp(-1.*ks*ks) - 4.*pow(ks, 3.)*exp(-ks*ks));

    else{
        return 1./5. - ks*ks/7. + pow(ks, 4.)/18.;
    }
}


double muOrderSix(double ks){
    if(ks>0.335168)  return pow(16.*pow(ks, 7.), -1.)*(15.*pow(pi, 0.5)*gsl_sf_erf(ks) - 2.*ks*exp(-ks*ks)*(15. + 10.*ks*ks + 4.*pow(ks, 4.)));

    else{
         return 1./7. - ks*ks/9. + pow(ks, 4.)/22.;
    }
}

//Mike
double kaiserGauss_Monofactor(double ks, double beta){
    return muOrderZero(ks) + 2.*beta*muOrderTwo(ks) + beta*beta*muOrderFour(ks);
}

//Mike
double kaiserGauss_Quadfactor(double ks, double beta){
    return (5./2.)*(-1.*muOrderZero(ks) + muOrderTwo(ks)*(3. - 2.*beta) + muOrderFour(ks)*(-beta*beta + 6.*beta) + 3.*beta*beta*muOrderSix(ks));
}


//FROM MIKE 
double haloModel_pk(double k, int order){
    // assumes linear bias of 1, therefor beta = f. Here 2. is the standard deviation of the gaussian variate added to z co-ordinate for non-linear RSD.
    double Interim;
  
    
    Interim = 0.0; 
    Interim = volume/particle_no;  
     
	//1halo only                                      				
	//return Interim*pow(ukm_nfw_profile(k), 2.0);

	//2halo only
	//return Interim + splint_generic(&Pk_current, k) * exp(-k*k*R_nl*R_nl);
	//return Interim + splint_generic(&Pk_current, k);
	
	double onehalo_term = 0.0;
	double twohalo_term = splint_generic(&Pk_current, k);
	//twohalo_term *= exp(-k*k*R_nl*R_nl); 
	
	//2halo + redshift space
	if (redshift_space) {
		if (order== 0) {
			
			//beta = fg because no galaxy bias assumed
			//onehalo_term += pow(ukm_nfw_profile(k), 2.0)*kaiserGauss_Monofactor(0.0, 0.0);		
			twohalo_term *= kaiserGauss_Monofactor(0.0, fg); //k dependence only with velocity dispersion	
		
		} else if (order == 2) {	
		
			//onehalo_term += pow(ukm_nfw_profile(k), 2.0)*kaiserGauss_Quadfactor(0.0, 0.0);
			twohalo_term *= kaiserGauss_Quadfactor(0.0, fg); //k dependece for lorentz/gaussian 	
		
		} else {
			printf("only monopole (order 0) and quadrupole (order 2) available \n");
			exit(0);
		}
	}
	
	return Interim + twohalo_term; 
	
	//1halo + 2halo
	//return Interim*pow(ukm_nfw_profile(k), 2.0) + pow(Dplus, 2.0)*splint_Pk_model(Pk_spline_model, k)*volume; 
}

//just outputs the halo model power spectrum to a given file
int haloModel_out(char outFile[], double kmin, double kmax, int bins) {
	double k, logdk, logk, logkmin;
	double model_mono, model_quad;	
	FILE *f = fopen(outFile, "w");	
	
	logkmin = log10(kmin);
	logdk = (log10(kmax) - logkmin) / (double)(bins-1);
	for (int i = 0; i < bins; i++) {
		logk = logkmin + (double)i * logdk;
		k = pow(10.0, logk);		
		model_mono = haloModel_pk(k, 0);
		model_quad = haloModel_pk(k, 2);
		fprintf(f, "%le \t %le \t %le \n", k, model_mono, model_quad);
	}
	fclose(f);	
	return 0;
}

int pk_out_halomodel(char outFile[], double kmin, double kmax, double kfold) {
	char out_temp[] = "/home/jonas/Testing_GR/code/Jonas/output/HOD_rands/temp_pk.dat";
	char out_temp_folded[] = "/home/jonas/Testing_GR/code/Jonas/output/HOD_rands/temp_pk_folded.dat";
	//halos_run(out_temp, 1.0);
	//halos_run(out_temp_folded, 4.0);
	//combine_folded(out_temp, out_temp_folded, outFile, kmin, kmax, kfold, true);	
	printf("combined\n");
	remove(out_temp);
    remove(out_temp_folded);
	return 0;
}

double I_integral_HOD(double t, void *params) {
	params = params;
	double f_t = log(1.0 + t) - t/(1.0 + t);
	return f_t/(pow(t,3.0)*pow(1.0 + t,2.0));
}

/*
//according to De la Torre http://arxiv.org/pdf/1202.5559v4.pdf
//which they got from VDB http://arxiv.org/pdf/astro-ph/0404033v2.pdf
int HOD_velocities(double halo_mass) {


	
	//FOR SATELLITES:
	double result_I, error, rho_nfw, f_cmd, rs, sigma_sq_sat, vx, vy, vz, sigma;
	rs = R_virial_global/cdm_global;
	
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
	gsl_function F;
  	F.function = I_integral;
  	F.params = NULL;   	
	
	for (int i = (particle_no - satellite_no); i < particle_no; i++) {
		 	
	  	gsl_integration_qags(&F, r/rs, 100000.0, 0, 1e-6, 1000, w, &result_I, &error); 
	 	
		//rho_nfw = pow(cdm_global*r/R_virial_global,-1.0)*pow(1.0 + cdm_global*r/R_virial_global,-2.0);
	
		f_cdm = log(1.0 + cdm_global) - cdm_global/(1.0 + cdm_global);	
		sigma_sq_sat = (G*halo_mass*cdm_global*cdm_global*r/(R_virial_global*R_virial_global*f_cdm))*pow(1.0 + cdm_global*r/R_virial_global, 2.0)*result_I;
	
		//this is the relative velocity to the halo, so add together with the host halo velocity
		sigma = sqrt(sigma_sq_sat);
		vx = gsl_ran_gaussian(gsl_ran_r, sigma_sq_sat);
		vy = gsl_ran_gaussian(gsl_ran_r, sigma_sq_sat);
		vz = gsl_ran_gaussian(gsl_ran_r, sigma_sq_sat);
		
		particles[i][4] = 
	
	}

	return 0;
}
*/

//omega_m in current cosmology is declared in main.c
int scale_masses(double s, double omega_m_primed) {
	for (int i = 0; i < particle_no; i++) {
		//OMEGA M DEPENDS ON Z!!!!!!!!!!!
		particles[i][6] = pow(s, 3.0)*omega_m_primed*particles[i][6]/omega_m;
	}
	return 0;
}

//a is declared in main.c
int scale_velocities(double s, double H, double H_primed, double a_primed) {
	double fg, fg_primed;
	//fg = d lng/ d lna
	for (int i = 0; i < particle_no; i++) {
		particles[i][3] = particles[i][3]*s*(H_primed*a_primed*fg_primed)/(H*a*fg);
		particles[i][4] = particles[i][4]*s*(H_primed*a_primed*fg_primed)/(H*a*fg);
		particles[i][5] = particles[i][5]*s*(H_primed*a_primed*fg_primed)/(H*a*fg);	
	}
	return 0;
}

//actual "running" code for randoms, generates centrals and then adds satellites assuming given halo mass
int randoms_halo_model(char root_path[], int run_no) {
	double box_mass, cos_theta_sat, theta_sat, phi_sat, r_sat, unif_rnd, mean_r_sat;
	double x_sat, y_sat, z_sat, erf_central, mean_N_central, mean_N_satellite;
	//velocity assignment for satellites variables
	double result_I, error, rho_nfw, f_cdm, sigma_sq_sat, sigma_sq_sat_tcdm, vx, vy, vz, sigma;
	int satellite_number, has_central, particle_no_init, total_sats;
	
	char final_out[100];
	sprintf(final_out, "%s%s%d%s", root_path, "finaltest/randoms_final_", run_no, ".dat");
	char onehalo_out[100];	
	sprintf(onehalo_out, "%s%s%d%s", root_path, "onehalotest/randoms_onehalo_", run_no, ".dat");
	char twohalo_out[100];	
	sprintf(twohalo_out, "%s%s%d%s", root_path, "twohalotest/randoms_twohalo_", run_no, ".dat");
	
	//box_mass = rho_bar * volume;
	//halo_mass = box_mass / particle_no;	
	//printf("box mass / particle number: %le\n", box_mass / (double) particle_no);	
		
	//save initial number for loops, as particle number will change with adding satellites
	central_no = particle_no;
	particle_no_init = particle_no;

	//centrals	
	populateParticles();	
	

	
	//CENTRALS GET THE VELOCITIES OF THEIR HOST HALOES
	/*
	for (int i = 0; i < particle_no; i++) {
		particles[i][0] = halo_vel...
	
	}
	*/
	
	
	//2halo term, restoring displacements from Pk model, wont need this memory so free it(clogs up laptop)
	//initialise ZA fft array
	printf("starting ZA\n");
	//fftw_allocateMemory_ZA();	
	//restore_disp_ZA();	
	//fftw_free(displacements);
	//fftw_free(displacements_fourier);	
	printf("ZA done\n");
	
	/*
	//1halo term
	//satellites, each central is assigned randomly from distribution
	//stored in the no_satellites array - how many satellites for each central
	//find the mean number of satellites at this halo mass
	mean_N_satellite = pow((halo_mass - M0)/M1, alpha);	
	mean_N_sats = mean_N_satellite;
	mean_r_sat = 0.0;
	total_sats = 0;
	
	int* no_satellites = calloc(particle_no, sizeof(*no_satellites));
	double mean_N_minus = 0.0;
		
	for (int i = 0; i < central_no; i++) {
		satellite_number = gsl_ran_poisson(gsl_ran_r, mean_N_satellite);
		mean_N_minus += (double) (satellite_number * (satellite_number - 1));
		//satellite_number = 16;
		no_satellites[i] = satellite_number;
		total_sats += satellite_number;
	}
	printf("adding %d satellites\n", total_sats);
	
	satellite_no = total_sats;
	particle_no += total_sats;
	mean_N_minus /= (double) central_no;
	
	//realloc to add memory for satellites
	particles = realloc(particles, particle_no * sizeof(*particles));
	//calloc zeros out the values of satellites
	for (int i = central_no; i < particle_no; i++) {
		particles[i] = (double*) calloc(7, sizeof(**particles));
	}		
	
	//prepare integral for velocity assignment of satellites		
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
	gsl_function F;
  	F.function = &I_integral_HOD;
  	F.params = 0;   
  	//end velocity assignment stuff	
				
	//add satellites, this counter is to place satellites in the now expanded particles array
	int satellite_index = particle_no_init;
	for (int i = 0; i < central_no; i++) {
		//each central i will have no_satellites[i] as found before
		for (int j = 0; j < no_satellites[i]; j++) {
			//cos theta is uniform random in [-1,1]
			cos_theta_sat = randomDouble(-1.0, 1.0);
			theta_sat = acos(cos_theta_sat);
			
			//phi is uniform random in [0,2pi]
			phi_sat = randomDouble(0.0, 2.0*M_PI);
			
			//selects a cumulative sum in range [cumsum_min, cumsum_max]
			//splint (inverse) then relates that cumsum to some r value between 0 and r_virial
			unif_rnd = randomDouble(nfw_spline_reverse.xmin, nfw_spline_reverse.xmax);

			r_sat = splint_generic(nfw_spline_reverse, unif_rnd);			
			//IF USING MIKES NFW CUMSUMS, UNCOMMENT BELOW, as r will be in units of rs
		    //r_sat *= R_virial_global;
			//r_sat /= cdm_global;
			
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
			
			gsl_integration_qags(&F, r_sat/rs_global, 1000.0, 0, 1e-5, 1000, w, &result_I, &error); 
	
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
	printf("mean r of sats: %lf, <N(N-1)>: %lf\n", mean_r_sat, mean_N_minus); 	
	*/
	
	//outputs
	printf("starting ffts\n");
	
	//allocate fftw memory for P(k), set the plan
	//fftw_allocateMemory();	
	//pk_out_halomodel(onehalo_out);
	pk_out_halomodel(twohalo_out, 0.03, 0.6, 0.15);
	//pk_out_halomodel(final_out);
	//fftw_free(overdensity);
	//fftw_free(overdensity_fourier);
	printf("ffts done\n");
	
	//for running multiple times, must clear memory		
	clearMemory();
	particle_no = particle_no_init;
	return 0;
}








//FOR OUTPUT OF A NEW CATALOGUE OF GALAXIES BASED ON AN INPUT HALO CATALOGUE
/*
int populate_halos(char out_centrals[], char out_satellites[]) {

	double erf_central;
	unsigned int has_central;
	int satellite_no;
	
	
	//R_virial = 10.;
	
	char ovl_in[] = "/home/jonas/Testing_GR/code/Jonas/output/variances/var0_logstep.dat";	
	char format_variances[] = "%le \t %le \n";
	variance_spline = input_spline(ovl_in, format_variances, true); //reversed
	R_mstar = splint_generic(variance_spline, delta_c * delta_c);
	Mstar = 4.0 * M_PI * pow(R_mstar, 3.0) * rho_bar / 3.0;	
	//printf("mstar is: %lf\n", Mstar);
	//cdm = (c0 / (1. + z_global)) * pow(m / Mstar, beta);
	//R_virial = pow(3.*m/(4.*M_PI*rho_bar*del_nl), 1./3.);
	
	
	
	char sim_in[100];
	sprintf(sim_in, "%s%s", mock_directory, "AllHalo_GR_1500Mpc_a0.7_1.txt");	
	//read_full_xyz(sim_in, 1, mocks_format);	
	//particles read in, positioned	
	
	//calculate cumulative sums for NFW distribution
	//rho_nfw_out();
	
	char cumsum_out[] = "/home/jonas/Testing_GR/code/Jonas/output/nfw_cumsums/mass_1E14";
	SplineInfo nfw = input_spline(cumsum_out, "%le \t %le \n", true);
	double norm = nfw.xmax;
	double rnd = randomDouble(0.0, 1.0);
	double rnd2 = randomDouble(0.0, 1.0);
	printf("for a normalised random: %le, the r is: %lf\n", rnd*norm, splint_generic(nfw, rnd*norm));
	printf("for a normalised random: %le, the r is: %lf\n", rnd2*norm, splint_generic(nfw, rnd2*norm));
	
	
	
	int massBins = 50;
	char cumsums_dir[] = "/home/jonas/Testing_GR/code/Jonas/output/nfw_cumsums/";
	char cur_file[100];
	SplineInfo* nfw_splines = malloc(massBins*sizeof(*nfw_splines));

	for (int i = 0; i < massBins; i++) {
		sprintf(cur_file, "%s%s%d", cumsums_dir, "cumsum_", i);
		nfw_splines[i] = input_spline(cur_file, "%le \t %le \n", true);		
	}
	//printf("splint: %lf\n", splint_generic(nfw_splines[0], 0.25));	
	//printf("splint result: %lf\n", splint_generic(nfw_splines[0], 0.26));
	
	
	
	
	//char new_catalogue[] = "/home/jonas/Testing_GR/code/Jonas/output/HOD/centrals_only.dat";
	FILE *centrals = fopen(out_centrals, "w");
	FILE *satellites = fopen(out_satellites, "w");
	
	init_rng();
	for (int i = 0; i < particle_no; i++) {
		halo_mass = particles[i][6];
		erf_central = gsl_sf_erf((log10(halo_mass) - log_Mmin)/sigma_logm);
		mean_N_central = 0.5 * (1. + erf_central);
		mean_N_satellite = pow((halo_mass - M0)/M1, alpha);
		has_central = gsl_ran_bernoulli(gsl_ran_r, mean_N_central);
		if (has_central == 1) {
			fprintf(centrals, galaxies_format, particles[i][0], particles[i][1],
			particles[i][2], particles[i][3], particles[i][4], particles[i][5]);
			fprintf(centrals, "\n");
		}
		satellite_no = gsl_ran_poisson(gsl_ran_r, mean_N_satellite);
				
		
		for (int j = 0; j < satellite_no; j++) {
			cos_theta_sat = randomDouble(-1., 1.);
			theta_sat = acos(cos_theta_sat);
			phi_sat = randomDouble(0., 2.*M_PI);
			//r_sat = pow(randomDouble(0.000001, R_virial), 1./3.);
			cdm_global = (c0 / (1. + z_global)) * pow(halo_mass / Mstar, beta);
			R_virial_global = pow(3.*halo_mass/(4.*M_PI*rho_bar*del_nl), 1./3.);
			
			x_sat = r_sat * sin(theta_sat) * cos(phi_sat) + particles[i][0];
			y_sat = r_sat * sin(theta_sat) * sin(phi_sat) + particles[i][1];
			z_sat = r_sat * cos(theta_sat) + particles[i][2];	
			fprintf(satellites, galaxies_format, x_sat, y_sat, z_sat, particles[i][3],
			particles[i][4], particles[i][5]);
			fprintf(satellites, "\n");		
		}
		
	}

	fclose(centrals);
	fclose(satellites);

	return 0;
}
*/
