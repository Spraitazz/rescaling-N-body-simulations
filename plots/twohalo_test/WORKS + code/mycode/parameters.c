int set_params() {

	//cosmology-specific variables	
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
	//rho_bar = rho_cosm*a*omega_m;
	
	
	//HOD variables at Mmax of -20.0
	log_Mmin = 11.83;
	sigma_logm = 0.25;
	log_M0 = 12.35;
	log_M1 = 13.08;
	M0 = pow(10.0, log_M0); //[M_sol]
	M1 = pow(10.0, log_M1); //[M_sol]
	alpha = 1.00;
	delta_c = 1.686;
	c0 = 11.0;
	beta = -0.13;
	del_nl = 200.0;			
	
	//target cosmology
	//omega_m_primed = 0.24;		

	//for D+ calculation
	lambda_F4 = 1.0/0.042; //[h Mpc^-1]	
	lambda_F5 = 1.0/0.133;
	lambda_F6 = 1.0/0.419;
	lambda_GR = 0.0;
	fR0_F4 = -1e-4;
	fR0_F5 = -1e-5;
	fR0_F6 = -1e-6;
	fR0_GR = 0.0;
	
	
	
	
	
	return 0;	
}
