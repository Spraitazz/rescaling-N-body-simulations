
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
	bool vary_z_current;
	double R1_primed;
	double R2_primed;
	SplineInfo* variance_const;
	SplineInfo** variances_varz;
} Multimin_Params;

typedef struct {	
	double s;
	double z_var;
	bool vary_z_current;
	SplineInfo* variance_const;
	SplineInfo** variances_varz;
} Dsq_Params;

typedef struct {
	double xmin;
	double xmax;
	int bins;
	double dx;
	int bin_type;
} BinInfo;

typedef struct {
	double h;
	double H0;
	double omega_m;
	double omega_v;// = 1.0 - omega_m;
	double omega_r;
	double omega;
	double fg;	
} Parameters;
