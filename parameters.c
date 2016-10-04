int set_params() {

	//cosmology-specific parameters
	rho_cosm = 2.78E+11; //[h^2 M_sol Mpc^-3]
	//current cosmology
	h = 0.678; //dimensionless hubble parameter
	H0 = 100.0*h; //[kms^-1 Mpc^-1]
	omega_m = 0.24;
	omega_v = 1.0 - omega_m;
	omega_r = 0.0;
	omega = 1.0;	
	G = 6.673e-11; //[N m^2 kg^-2]
	fg = pow(omega_m, 0.545);
	
	//HOD parameters at Mmax -19.5
	//http://arxiv.org/pdf/1005.2413v2.pdf
	log_Mmin = 11.57;
	sigma_logm = 0.17;
	log_M0 = 12.23;
	log_M1 = 12.87;
	M0 = pow(10.0, log_M0); //[M_sol]
	M1 = pow(10.0, log_M1); //[M_sol]
	alpha = 0.99;
	
	//critical linear overdensity for top-hat collapse
	delta_c = 1.686;
	
	//mass-concentration-redshift relation for CDM haloes
	c0 = 11.0;
	beta = -0.13;
	
	//virialisation overdensity treshold
	del_nl = 178.0;		

	//for D+ calculation
	lambda_F4 = 1.0/0.042; //[h Mpc^-1]	
	lambda_F5 = 1.0/0.133;
	lambda_F6 = 1.0/0.419;
	lambda_GR = 0.0;	
	
	
	
	return 0;	
}
