//size limits for box
double volume_limits[3];
double volume_limits_save[3];
double volume;

//jenkins folding
double Jenkins_foldfactor = 1.0;

//particle number in "universe"
int particle_no, halo_no, central_no, satellite_no;
int max_sats;
Particle_Catalogue *all_haloes;

//grid size, number of cells 
int cells[3];
int cells_displ[3];

//particles and overlay grid
//Particle* particles;
double** satellites;
double** centrals;
Particle* particles_save;
double* grid;

//data input formats
const char mocks_format[] = "%lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %le";
const char galaxies_format[] = "%le \t %le \t %le \t %le \t %le \t %le";
const char cube_format[] = "%lf \t %lf \t %lf \t %*f \t %*f \t %lf";
const char spline_format_Pk[] = "%lf \t %lf \t %*f \t %*d \n";
const char monopoles_rescaling_format[] = "%lf \t %lf \t %*f \n";
const char Pk_model_format[] = "%le \t %le \n";
const char variance_format[] = "%le \t %le \n";

//fftw data
fftw_plan overdensity_plan, ZA_displacements_plan;
fftw_complex *overdensity, *overdensity_fourier;
fftw_complex *displacements_fourier, *displacements;

int ph;

//Pk calculation data
int spectrum_size, spectrum_size_final;
double* fk_abs_squares;
double* monopole;
double* quadrupole;
double* hexadecapole;
int* k_number;
bool* toRemove_bins;
double* bin_k_average;
double* bin_k_min;
double* bin_k_max;
//1st element is sum of L2 for that bin, 2nd is sum for [L2 squared], 3rd is sum of [L2 * fourier coefficient abs squared]
double** binsums_L2;
double** binsums_L4;
double WindowFunc_pow;


char pk_format[] = "%lf \t %le \t %le \t %le \t %d \t %lf \t %lf \n";
char k_Pkmono_format[] = "%lf \t %le \t %*e \t %*e \t %*d \t %*f \t %*f \n";

char pk_format_folding[] = "%lf \t %le \t %le \t %*d \t %*f \t %*f \n";


//catalogue input files
char mock_directory[] = "/home/jonas/Testing_GR/Mocks/L1500_np1024_a0.7_halo";
const char* home_directory;
char temp_directory[100];
char cube_file[] = "/home/jonas/Testing_GR/Mocks/cube_gal_-20.0_vel.dat";
//char cube_directory[];

//HOD population stuff
double cdm_global, R_virial_global, rs_global, mean_N_sats;
double alpha, log_M1, M0, M1, delta_c, R_mstar, Mstar, c0;
double beta, z_global, z_current, z_target, del_nl, b_eff;

//cosmological parameters
double rho_cosm, omega_m_0, omega_m_primed, rho_bar, H0;
double omega_v_0, omega_r_0, omega, H, fg, fg_primed, G, h;
double log_Mmin, sigma_logm, log_M0;




//int k_bins, z_bins, R_bins, Dplus_bins, current_gravity, target_gravity;
int Dplus_bins;

double weighted_shotnoise_sats;



//rescaling
double lambda_F4, lambda_F5, lambda_F6, lambda_GR;
double fR0_F4, fR0_F5, fR0_F6, fR0_GR;
double R1_primed, R2_primed, Mmin, Mmax;
double z_guess;
double dsq_minimized_value, s, z_rescaled;
bool vary_z_current;
//double k_min, k_max, z_min, z_max;


//for ZA
double R_nl;
double sigma_exp_smoothed;
bool ZA_mem = false;


int prep_oneHalo(void);
int oneHalo(void);





/*
SplineInfo *variance_spline_current, *variance_spline_current_reversed;
SplineInfo *variance_spline_target, *variance_spline_const;
SplineInfo *nfw_spline_reverse;
SplineInfo *nfw_cumsums_reversed;
SplineInfo *redshift_spline_reversed;
SplineInfo *Pk_current, *Pk_target;
SplineInfo **Pk_splines_zBins;
SplineInfo **variance_splines_zBins;
*/
//SplineInfo_Extended *Pk_target_extended, *Pk_current_extended;
//SplineInfo_Extended **Pk_splines_zBins_extended; 

//BinInfo *rescaling_k_bin_info, *rescaling_z_bin_info, *rescaling_R_bin_info, *HOD_mass_bin_info;
//BinInfo *Pk_calc_bins;

//functions.c
double randomDouble(double min, double max);
int PBC(double *x, double *y, double *z);
int countLines(FILE *f);
int arr_ind(int x, int y, int z);
int arr_ind_displ(int x, int y, int z);
int mass_min_max(Particle_Catalogue *catalogue);
double m_to_R(double mass, Parameters *params);
double R_to_m(double R, Parameters *params);
Spline* input_spline_values(int lines, double* x_vals, double* y_vals, int model);
Spline* input_spline_file(char inFile[], const char format[], int reverse, int model);
Spline* prep_spline_generic(BinInfo *binInfo, int reverse, double (*func)(double));
double splint_generic(Spline *spline, double x);


//memory.c
Particle* initParticles(int particle_no);
int init_rng(void);
int clearParticles(Particle** particles);
int clearGrid(void);
int overdensity_clearMemory(void);
int overdensity_allocateMemory(void);
int overdensity_freeMemory(void);
int ZA_clearMemory(void);
int ZA_freeMemory(void);
int ZA_allocateMemory(void);
int clearPk(void);
int clearMemory(void);
int memCheck(void);
int init_spectrum_storage(int bins);
void calloc2D_double(double*** arr, int dim1, int dim2);
void free2D_double(double*** arr, int dim2);

//power_spectrum.c
int populateParticles(Particle** particles, int particle_no);
int toRedshift(Particle_Catalogue *cat, Parameters *params);
int initGrid(void);
int populateGrid(Particle_Catalogue *catalogue, int grid_func);
int populateGridNGP(Particle_Catalogue *catalogue);
int populateGridCIC(Particle_Catalogue *catalogue);
int gridIntoOverdensity(void);
int fftw_overdensity_calculate(void);
int overdensity_pk();
int pk_logbin(int overdensity_index, double k, double mu, BinInfo* Pk_binInfo);
int multipole_calc(void);
int pk_to_file_logbin(const char filename[]);
int delsq_to_file_logbin(const char filename[]);

//Cubic spline
void spline(double  x[], double  y[], int n, double yp1, double ypn, double y2[]);
void splint(double xa[], double ya[], double y2a[], int n, double x, double *y);
int powerlaw_regression(int N, double rmin, double rmax, double sign, double ri[], double Di[], double* norm, double* index);

//halo_model.c
int scale_velocities(double s, Particle_Catalogue *cat, Parameters *cur_params, Parameters *targ_params);
int scale_masses(double s, double omega_m_primed);
int nfw_cumsums_catalogue(int bins, const char root_dir[]);
double halo_model_Pk(double k, double shot_noise, double twohalo_prefactor, int redshift_space, int order, Spline *Pk_spline, double fg);
Spline* prep_nfw_cumsums(double halo_mass, int bins);
int haloModel_out(char out[], BinInfo *k_bins, double shot_noise, double twohalo_prefactor, int redshift_space, Spline *Pk_spline, Parameters *par);
double cumsum_nfw(double r);
int set_halo_params(double halo_mass, Spline *variance_spline, Parameters *params);
int set_ZA_params(void);
double ukm_nfw_profile(double k);
double I_integral_HOD(double t, void *params);

//catalogue_input.c
int read_full_xyz(const char filename[], int skipLines, const char format[]);

//HOD.c
int add_centrals(void);
int prep_HOD_Pk(bool sats);
int add_satellites(void);
double weighted_shotnoise(void);
int prep_bin_info(void);
int HOD_main(void);

//rescaling.c
Spline* extend_spline_model(SplineInfo *spline);
Spline* extend_spline(SplineInfo *spline, int model);
//double splint_Pk_model(SplineInfo_Extended *pk_model, double k);
double splint_Pk(Spline *Pk_spline, double k);
int halo_mass_func(void);
int rescale(void);
int rescale_testRun(void);
int rescale_catalogue(void);
int rescale_model(void);

//rescaling_functions.c
double delta_sq_int_func(double R, void *params);
double OLV(double k, void *params);
double integrate_OLV(double R, Spline *pk_spline);
double delta_sq_rms(const gsl_vector *inputs, void *params);
int dsq_multimin(bool vary_z_cur, double z_init, Spline* variance_const, Spline** variances_varz, BinInfo *z_binInfo);
int dsq_multimin_test(void);
Spline* prep_variances_test(BinInfo *R_bin_info, Spline *Pk_spline, double s_test);
double delsq_lin(Spline *Pk_spline, double k);
double expected_variance_smoothed(double R_nl_this, Spline *Pk_spline);
double OLV_smoothed(double k, void *params);

Spline* prep_Pk_constz(int gravity, double z_from, double z_to, Spline *Pk_spline_in, BinInfo *k_bin_info, BinInfo *z_binInfo, Parameters *params);
Spline* Dplus_spline(int gravity, double k, BinInfo *z_binInfo, Parameters *params);
Spline* prep_variances(BinInfo *R_bin_info, Spline *Pk_spline);
int generate_sigma_plots(void);


//Dplus.c
double z_to_t(double z, Parameters *params);
double z_to_a(double z);
double a_to_z(double a);
double a_to_t(double a, Parameters *params);
double H_a(double a, Parameters *params);
double z_to_t_int(double z, void * params);
Spline* prep_redshift_spline(BinInfo *z_binInfo, Parameters *params);
double t_to_z(double t, Spline *redshift_reverse);
double H_z(double z);
int growth_function(double t, const double y[], double f[], void *params);
int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
int Dplus_calc(int gravity, double k, double** zs, double** Dpluses, BinInfo *z_binInfo, Parameters *params);

//bias
double m_to_v(double mass, Parameters *params, Spline* variance_spline);
double f_v(double v);
double ln_fv_deriv(double v);
double v_to_m(double v, Bias_Params *bias_params);
double b_v(double v);
double b_m(double m, Parameters *params, Spline* variance_spline);
double bf_over_m(double v, void * params);
double f_over_m(double v, void * params);
int bias_plot(double z);
double calc_b_eff(Parameters *params, Spline* variance_spline, double Mmin, double Mmax);

//fold.c
Particle* Jenkins_fold_volume(Particle_Catalogue *catalogue);
int Jenkins_restore(Particle_Catalogue *catalogue, Particle** particles_saved);
FoldingInformation* combine_folded(const char unfolded_path[], const char folded_path[], char pk_format[], char outFile[], double kmin, double kmax, double k_fold, int* foldPoints, bool print);

//covariance_matrices
int generate_covariance(char inPath[], char outPath[], int sim_number, int startIndex);

//runs.c
int pre_covariance_run(const char inputData[], const char input_format[], const char outFile[], bool save_files, bool fold, int skipLines);
int randoms_fullRun(int particle_number, const char outFile[], bool fold, bool test_mode);
int haloes_measure_Pk(Particle_Catalogue *catalogue, char outFile[], double foldfactor, int grid_func);
	
//tests.c
int test_run_randoms(void);
int full_test(bool verbose);
int memory_test(bool verbose);
int overdensity_pk_test(bool verbose);
int randoms_populate_test(bool verbose);
int test_file_io(bool verbose);
int test_random(bool verbose);


//#JustGSLThings
const gsl_rng_type* gsl_ran_T;
struct gls_rng* gsl_ran_r;
gsl_integration_workspace* z_to_t_workspace;
size_t z_to_t_workspace_size = 1000;
gsl_integration_workspace* dsq_workspace; 
gsl_integration_workspace* OLV_workspace;
size_t OLV_workspace_size = 2000;
gsl_integration_workspace* onehalo_workspace;
gsl_integration_cquad_workspace* dsq_workspace_cquad;






