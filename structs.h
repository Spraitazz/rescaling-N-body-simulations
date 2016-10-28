
/*
typedef struct {
	double** folded_values;
	double** unfolded_values;
	double kmin;
	double kmax;
	int unfold_last_ind;
	int fold_first_ind;
	int min_k_ind;
	int aliasing_index;
	int lines;
	bool max_reached;
} FoldingInformation;
*/

typedef struct {
	double x;
	double y;
	double z;
	double vx;
	double vy;
	double vz;
	double mass;
} Particle;

typedef struct {
	Particle* particles;
	int particle_no;
} Particle_Catalogue;

typedef struct {	
	int uf_first_ind;
	int uf_last_ind;
	int f_first_ind;
	int f_last_ind;
	bool max_reached;
	bool success;	
} FoldingInformation;

typedef struct {
	int folds;
	int** foldPoints;
} Foldings;

typedef struct {
	double* x_vals;
	double* y_vals;
	double* coeffs;
	int lines;
	double xmin;
	double xmax;
} SplineInfo;

typedef struct {
	SplineInfo *splineInfo;
	double pk_loA;
	double pk_hiA;
	double pk_lon;
	double pk_hin;
	int model;
} Spline;



/*
typedef union {
	SplineInfo *sp;
	SplineInfo_Extended *sp_ext;
} Spline;
*/



typedef struct {
	double xmin;
	double xmax;
	int bins;
	double dx;
	int bin_type;
} BinInfo;

typedef struct {
	double H0;
	double omega_m_0;
	double omega_v_0;// = 1.0 - omega_m;
	double omega_r_0;
	double omega;
	double gamma; //f = omega_m ^ gamma
	double z;	
} Parameters;

typedef struct {
	Parameters *parameters;
	Spline *variance;
} Bias_Params;

typedef struct {
	Spline *Pk_spline;
	double R;	
} OLV_parameters;

typedef struct {
	bool vary_z_current;
	double R1_primed;
	double R2_primed;
	Spline *variance_const;
	Spline **variances_varz;
	BinInfo *z_binInfo;
} Multimin_Params;

typedef struct {	
	double s;
	double z_var;
	bool vary_z_current;
	Spline *variance_const;
	Spline **variances_varz;
	BinInfo *z_binInfo;
} Dsq_Params;

typedef struct {
	double lambda0;
	double k;
	double t0;
	Spline *redshift_reverse;
	Parameters *cosmo_params;
} Minim_Params;
