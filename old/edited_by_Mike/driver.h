//size limits for Pk
double limit_x;
double limit_y;
double limit_z;

//particle number in "universe"
int particle_no;

//grid size
int cells_x;
int cells_y;
int cells_z;

//particles and overlay grid
double** particles;
double* grid;

//ft of overdensity (without extern throws compiler warning??)
double** overdensities_ft;  //[][2];

double power_spectrum[100];

//functions
double randomDouble();
int populateParticles();
int initGrid();
int populateGrid();
int gridIntoOverdensity();
int fourier_coefficient();
int overdensity_pk();

const gsl_rng_type* gsl_ran_T;
gsl_rng*            gsl_ran_r;
