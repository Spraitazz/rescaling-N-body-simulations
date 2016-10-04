//size limits for box
double volume_limits[3];
double volume_limits_save[3];
double volume;

//jenkins folding
double Jenkins_foldfactor = 1.0;

//particle number in "universe"
int particle_no, halo_no, central_no, satellite_no;

//grid size, number of cells 
int cells[3];
int cells_displ[3];

//particles and overlay grid
double** particles;
double** particles_save;
double* grid;

bool redshift_space;

/*

//Mike's halomodel.c
double* xDisplacement;
double* yDisplacement;
double* zDisplacement;
fftw_complex *outx, *outy, *outz, *H_k;
fftw_plan iplan_x, iplan_y, iplan_z;
int n0, n1, n2;
int DisplacementCalc();
int HaloCatalogue_NFWprofile(int halo_Number, int total_gal);
int prep_DisplacementCalc();
int prep_grid();
int free_pkRegression();


//Mikes power spectrum stuff
fftw_complex *smooth_overdensity, *H2_k;
fftw_plan p;

double*  kHexadecapole; 
double*  kQuadrupole;
double*  kMonopole;

double** polar_pk;

double* mean_modk;
double*  logk_limits        = NULL;
int* modes_perbin;
double* binnedPk;
double** d2_binnedpk;
double** twodim_pk;

int i, j, k;

FILE *output;
int polarPk_modeCount = 0;
int kbin_no;
int modeCount;
int fieldFlag;

char outFile_Pk_Mikes[100];

int polar_pkcount;

double nz_smoothRadius;

char root_dir[] = "/shome/jonasp/Work/Testing_GR/code/Jonas/output";



int MonopoleCalc(int modBinNumb, double mean_modBin[], double Monopole[], double** Array, int modeCount, char filepath[], double mu_lolimit, double mu_hilimit, int fileOutput);
int PkCorrections();
int PkCalc();
int prep_fftw();

*/


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
double* power_spectrum_monopole;
double* power_spectrum_quadrupole;
int* k_number;
bool* toRemove_bins;
double* bin_k_average;
double* bin_k_min;
double* bin_k_max;
//1st element is sum of L2 for that bin, 2nd is sum for [L2 squared], 3rd is sum of [L2 * fourier coefficient abs squared]
double** binsums_L2;
double WindowFunc_pow;


const char pk_format[] = "%lf \t %le \t %le \t %d \t %lf \t %lf \n";
char pk_format_folding[] = "%lf \t %le \t %le \t %*d \t %*f \t %*f \n";


//catalogue input files
char mock_directory[] = "/home/jonas/Testing_GR/Mocks/L1500_np1024_a0.7_halo/";
char* home_directory;
char cube_file[] = "/home/jonas/Testing_GR/Mocks/cube_gal_-20.0_vel.dat";
//char cube_directory[];

//HOD population stuff
double cdm_global, R_virial_global, rs_global;
//rescaling variables
double Mmin, Mmax, mean_N_sats, H0, omega_v, omega_r, omega, H, fg, fg_primed;
double rho_cosm, a, omega_m, omega_m_primed, rho_bar, log_Mmin, sigma_logm, log_M0, G, h;
double alpha, log_M1, M0, M1, delta_c, R_mstar, Mstar, c0, beta, z_global, z_current, z_target, del_nl, b_eff;

double lambda_F4, lambda_F5, lambda_F6, lambda_GR;
double fR0_F4, fR0_F5, fR0_F6, fR0_GR;
double R1_primed, R2_primed;

double sigma_exp_smoothed;

double dsq_minimized_value, s, z_rescaled;
bool vary_z_current, testmode;
double k_min, k_max, z_min, z_max, dk, dz;
int k_bins, z_bins, Dplus_bins, current_gravity, target_gravity;



//for ZA

double z_current_global, R_nl;

double* phases;

int prep_oneHalo(void);
int oneHalo(void);


//structs
typedef struct {
	double** folded_values;
	double** unfolded_values;
	int unfold_last_ind;
	int fold_first_ind;
	int min_k_ind;
	int foldover_index;
} FoldingInformation;

typedef struct {
	double* x_vals;
	double* y_vals;
	double* coeffs;
	int lines;
	double xmin;
	double xmax;
} SplineInfo;

typedef struct {
	SplineInfo *spline;
	double pk_loA;
	double pk_hiA;
	double pk_lon;
	double pk_hin;
	bool model;
} SplineInfo_Extended;

typedef struct {
	SplineInfo_Extended *spline;
	double R;	
} OLV_parameters;

typedef struct {
	double xmin;
	double xmax;
	int bins;
	double dx;
	bool logbin;
} BinInfo;


//splines
SplineInfo variance_spline_current, variance_spline_target, variance_spline_current_reversed;
SplineInfo nfw_spline_reverse;

SplineInfo* nfw_cumsums_reversed;
SplineInfo* Pk_splines_zBins; 
SplineInfo* variance_splines;
SplineInfo redshift_spline_reversed;

SplineInfo Pk_current;

SplineInfo_Extended Pk_spline_current, Pk_spline_target, Pk_current_extended;

//bins
BinInfo rescaling_k_bin_info, rescaling_z_bin_info, rescaling_R_bin_info;
BinInfo Pk_calc_bins;

//functions.c
double randomDouble(double min, double max);
int PBC(double *x, double *y, double *z);
int countLines(FILE *f);
int arr_ind(int x, int y, int z);
int arr_ind_displ(int x, int y, int z);
int mass_min_max(void);
double m_to_R(double m, double z);
double R_to_m(double R, double z);
SplineInfo input_spline_values(int lines, double* x_vals, double* y_vals);
SplineInfo input_spline_file(char inFile[], const char format[], bool reverse);
SplineInfo prep_spline_generic(double xmin, double xmax, int bins, bool logbin, bool reverse, double (*func)(double));
double splint_generic(SplineInfo *info, double x);

//memory.c
int initParticles(void);
int init_rng(void);
int clearParticles(void);
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
double init_spectrum_storage_len(int bins, double kmin, double kmax);
void calloc2D_double(double*** arr, int dim1, int dim2);
void free2D_double(double*** arr, int dim2);

//power_spectrum.c
int populateParticles(void);
int toRedshift(void);
int initGrid(void);
int populateGrid(void);
int populateGridCIC(void);
int gridIntoOverdensity(void);
int fftw_overdensity_calculate(void);
int overdensity_pk(bool randoms);
int pk_logbin(int index, double k, double logkmin, double logkmax, double logdk, double mu);
int multipole_calc(void);
int remove_bins(void);
int pk_to_file_logbin(const char filename[]);
int delsq_to_file_logbin(const char filename[]);

//Cubic spline
void spline(double  x[], double  y[], int n, double yp1, double ypn, double y2[]);
void splint(double xa[], double ya[], double y2a[], int n, double x, double *y);
int powerlaw_regression(int N, double rmin, double rmax, double sign, double ri[], double Di[], double* norm, double* index);

//halo_model.c
int scale_velocities(double s, double H, double H_primed, double a_primed);
int scale_masses(double s, double omega_m_primed);
int randoms_halo_model(char root_path[], int run_no);
int nfw_cumsums_catalogue(int bins, const char root_dir[]);
int pk_out_halomodel(char outFile[], double kmin, double kmax, double kfold);
double haloModel_pk(double k, int Order);
int prep_nfw_cumsums(char outFile[]);
int haloModel_out(char outFile[], double kmin, double kmax, int bins);
double cumsum_nfw(double r);
int set_halo_params(double halo_mass);
int set_ZA_params(void);
double ukm_nfw_profile(double k);
double I_integral_HOD(double t, void *params);

//catalogue_input.c
int read_full_xyz(const char filename[], int skipLines, const char format[]);

//rescaling.c
int prep_Pk_model_regression(SplineInfo *spline, SplineInfo_Extended *model_spline);
double splint_Pk_model(SplineInfo_Extended *pk_model, double k);
int displacements_ZA(int xyz);
int restore_disp_ZA(void);
int Pk_out(char outFile[]);
int get_Pk(char sim_path[], SplineInfo *spline);
int halo_mass_func(void);
int rescale(void);

//rescaling_functions.c
double delta_sq_int_func(double R, void *params);
double OLV(double k, void *params);
double integrate_OLV(double R, SplineInfo_Extended *pk_spline);
double delta_sq_rms(const gsl_vector *inputs, void *params);
int dsq_multimin(void);
double Dplus_k_z(double k, double z, bool z_const);
int prep_Pk_splines(int gravity, SplineInfo *Pk_spline_in, bool vary_z, double z_fixed, SplineInfo *Pk_spline_out, SplineInfo* *Pk_splines_zBins);
SplineInfo Dplus_spline(int gravity, double k);
int prep_variance_splines(SplineInfo_Extended *Pk_spline_in, SplineInfo *variance_spline_out, SplineInfo_Extended* *Pk_splines_zBins, SplineInfo* *variance_splines_zBins, bool vary_z, double z_fixed);
int prep_variances(BinInfo *R_bin_info, char OLV_out[], SplineInfo_Extended *pk_spline);
int generate_sigma_plots(void);

//Dplus.c
double z_to_t(double z);
double z_to_a(double z);
double a_to_z(double a);
double a_to_t(double a);
int prep_redshift_spline(double zmin, double zmax, int bins);
double t_to_z(double t);
double H_z(double z);
double a_to_t_int(double z, void * params);
int growth_function(double t, const double y[], double f[], void *params);
int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
int Dplus_calc(int gravity, double k, char outFile[], double** zs, double** Dpluses);

//bias
double m_to_v(double mass);
double f_v(double v);
double ln_fv_deriv(double v);
double v_to_m(double v);
double b_v(double v);
double bf_over_m(double v, void * params);
double f_over_m(double v, void * params);
int calc_b_eff(void);

//fold.c
int Jenkins_fold_volume(void);
int Jenkins_restore(void);
FoldingInformation combine_folded(const char unfolded_path[], const char folded_path[], char pk_format[], char outFile[], double kmin, double kmax, double k_fold, bool print);

//covariance_matrices
int generate_covariance(char inPath[], char outPath[], int sim_number, int startIndex);

//runs.c
int pre_covariance_run(const char inputData[], const char input_format[], const char outFile[], bool save_files, bool fold, int skipLines);
int Gaussian_randoms(int particle_number, const char outFile[], bool fold);
int noise_randoms(int particle_number, const char outFile[], bool fold, bool test_mode);
int randoms_fullRun(int particle_number, const char outFile[], bool fold, bool test_mode, bool gaussian);
int haloes_measure_Pk(char outFile[], double foldfactor, bool CIC);
	
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
gsl_integration_workspace* dsq_workspace; 
gsl_integration_workspace* OLV_workspace;
gsl_integration_workspace* onehalo_workspace;
gsl_integration_cquad_workspace* dsq_workspace_cquad;






