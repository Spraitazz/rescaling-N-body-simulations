// Available functions. 
int          comovDistReshiftCalc();

int          DisplacementCalc();

int          JenkinsCoordinates();
int          JenkinsFold(double original[], int lenArray, int axis);
int          ApplyJenkins();

double       invert_StefanoBasis(double centreRA, double centreDec, double* xval, double* yval, double* zval);

int          EvaluateGridParameters();
double       SolidAngleCalc(double decLowerBound, double decUpperBound, double raInterval);

int          ApplyFKPweights();

int          CalcCorrections();

int          randNGP();

// double        zChi_Integrand(double x);
int          CoordinateCalc();
int          NGPCalc();
int          PkCalc();

double       splintMatterPk(double EvalPoint);

double       (*pt2Pk)(double k)                 	= NULL;
double       (*pt2Pk_z)(double a, double k)       	= NULL;
double       (*pt2Xi)(double)       				= NULL;

double       (*pt2RSD_k)(double, double, int)       = NULL;
double       (*pt2RSD_r)(double, double, int)       = NULL;

double       (*pt2nz)(double)                       = NULL;

double       volavg_invnbar;
double       volavg_redshift;

int          VIPERSbasis(double centerRA, double centerDec, double xCoors[], double yCoors[], double zCoors[], int len);

int          Celestialbasis(double centerRA, double centerDec, double xCoors[], double yCoors[], double zCoors[], int len);

int          projectVIPERSsystem();


// Artificial window fn's for cubic run
int          FullCube();
int          PencilBeamSurvey(int xlow, int xhi, int ylow, int yhi);

double       interp_nz(double x);

int          FirstColumnCompare(const void* a, const void* b);
int          SecondColumnCompare(const void* a, const void* b);

double       ArrayMax(double a[], int n);
double       ArrayMin(double a[], int n);

double       SumDoubleArray(double array[], int len);

int          printInterpPk();

int          assignbinninginterval(double interval);

int          zCubeCreate();

// padded window fn. calculation.

int          assignNGPMemory();
int          assignFFTMemory();
int          assignNGPMemory();
int          assignFFTMemory();
int          assign2DPkMemory();
int          assignConvolutionMemory();

int          freeRand();
int          freeHOD();
int          freeNGP();
int          freeFFTw();
int          freePk();
int          free_sdltInterp();
int          free_linear();
int          freeClipped();
int          free2dPk();


int          setInputPk();


double       interp_comovingDistance(double z);
double       interp_inverseComovingDistance(double r);

// Pointers to interpolation functions. 
float        (*pt2zChiIntegrand)(float);
double       (*pt2interp_comovingDistance)(double)         = &interp_comovingDistance;
double       (*pt2interp_inverseComovingDistance)(double)  = &interp_inverseComovingDistance;

char         vipersHOD_dir[200];

/*
// VIPERS HOD mock parameters.              
// 500s mocks. Big multidark, http://www.multidark.org/MultiDark/Help?page=simulations corresponding to planck 1: http:arxiv.org/abs/1303.5076
double  Om_v      =    0.69;          // 0.73; value added mocks + HOD Cube. outdated unless using one of these. 
double  Om_r      =    0.00;          // 0.0;

// total matter, cdm+baryons. 
double  Om_m      =    0.31;          // 0.27;
double  Om_b      =   0.048;          // 0.0469;
double  Om_tot    =   1.000;          // 1.0;
double  h         =   0.673;          // 0.7;
double  sigma_8   =    0.82;          // 0.82;
double  ns        =    0.96;          // 0.95;

double  p_index   =     1.0;
double  pk_gamma  =     0.0;          // Bond Efstathiou shape parameter.
double  hz        =     0.0;
*/

fftw_complex *overdensity, *smooth_overdensity, *H_k, *H2_k;

double        H_0;
double        H_0inPerSec;
double        HubbleTime;

double        HubbleCnst(double z);

// Selection parameters. Volume limited sample between redshift 0.7 and 0.9
double        lo_zlim;
double        hi_zlim;
double          z_eff; 
double           aexp;

double        lo_MBlim;
double        hi_MBlim;  

// Array to hold the coordinate limits of the VIPERS survey. 
double        AxisLimsArray[2][3];
         
double        xCellSize;                                                     // Cell size, comoving distance, h^-1 Mpc
double        yCellSize;                                                     // Cell size, comoving distance, h^-1 Mpc
double        zCellSize;                                                     // Cell size, comoving distance, h^-1 Mpc

double        MinChi3;                                                      // h^-1 Mpc, Approximately, fig. 14, Guzzo et al.  2013
double        MaxChi3;                                                      // Redshift limited sample, 0.7 < z < 0.9
double        IntervalChi3;

double        hiChi;
double        loChi;

double        UpperRAlimit;
double        LowerRAlimit;

double        UpperDecLimit;
double        LowerDecLimit;

double        CellVolume;

// double*      densityArray        = NULL;
// double*      FKPweights          = NULL;
// double*      booldensity         = NULL;
double*      meanCellRedshift    = NULL;

double       TotalVolume         = 0.0;
double       TotalSurveyedVolume = 0.0;

double               fkpWeightedVolume   = 0.0;
double               fkpSqWeightsVolume  = 0.0;
double               fkpShotNoiseCorr    = 0.0;
double       analyticfkpWeightedVolume   = 0.0;

// Dimensions of the padded volume (TotalVolume) in units of CellSize. 
int          n0, n1, n2;      
int          loopCount;

// z - Comoving distance data points for interpolation. 
int          nPoints = 1000;
double       z_Array[1001];
double       ComovingDistance_z[1001];
double       redshiftInterval  = 0.015;

// Spline and splint holders. 
double       z_ComovingDistance_2derivatives[1001];
double       ComovingDistance_z_2derivatives[1001];

// Memory assignment for ZADE catalogue input. 

int                ch;
int        Vipers_Num;
int         fieldFlag;

int*                id  = NULL;
double*             ra  = NULL;
double*            dec  = NULL;
double*           zobs  = NULL;
double*           zcos  = NULL;
double*            M_B  = NULL;
double*          zflag  = NULL;
int*              type  = NULL;
int*          photoMask = NULL;


// Value added catalogue parameters. 
double*           zpec   = NULL;
double*          zphot   = NULL;

double*	         gal_z   = NULL;
double*            csr   = NULL;
double*       sampling   = NULL;  
double*     sampling35   = NULL;
char**       pointing    = NULL;
char**       quadrant    = NULL;
int*      flag_Nagoya    = NULL;
int*       flag_SSPOC    = NULL;
int*     flag_SSPOC35    = NULL;
double*       rand_sel   = NULL;
double*   fkp_galweight  = NULL;
double*   clip_galweight = NULL;


// derived parameters
bool*    Acceptanceflag = NULL;
double*     polarAngle  = NULL;
double*          rDist  = NULL;
double*          xCoor  = NULL;
double*          yCoor  = NULL;
double*          zCoor  = NULL;
double*          xVel   = NULL;
double*          yVel   = NULL;
double*          zVel   = NULL; 

// randoms.
double       alpha;

int          accepted_gals;
int          rand_number   = 0;
int          accepted_rand = 0;

double*      rand_ra     = NULL;
double*      rand_dec    = NULL;
double*      rand_chi    = NULL; 
double*      rand_x      = NULL;
double*      rand_y      = NULL;
double*      rand_z      = NULL;
double*      rand_weight = NULL;
bool*        rand_accept = NULL;

// TotalWeight is the sum of ZADE weight for the ZADE catalogue = Number of spec z + Number of used Photometric galaxies used (including compensation for sampling).
double       TotalZADEWeight   = 0.0;
double       MeanNumberDensity = 0.0;

int          boxlabel;
int          xlabel, ylabel, zlabel;

// FFTw calc parameters.
int          Index;
int          local_Index;
int          xNyquistIndex; 
int          yNyquistIndex;
int          zNyquistIndex;

double        kSq;
double        kIntervalx;
double        kIntervaly;
double        kIntervalz;

double        kmodulus;
double        mu;

double        kbinInterval;

double        k_x, k_y, k_z;

double        H_kReal;          
double        H_kImag;           

double        xNyquistWaveNumber;
double        yNyquistWaveNumber;
double        zNyquistWaveNumber;

// First column is mod k, second Pk.
// double**     PkArray            = NULL;
double**     twodim_pk        = NULL;

double**     muIntervalPk       = NULL;

double*      legendre2weights   = NULL;
// Binned Pk parameters.

// int          kBinNumb;
int          kbin_no;
int          Num_ModesInMuInterval = 0;
int*         modes_perbin        = NULL;

double       modkMax;

double*      del2               = NULL;
// double*      midKBin            = NULL;
// double*      meanKBin           = NULL;        
double*      mean_modk          = NULL;

double*      binnedPk           = NULL;
// double*      kBinLimits         = NULL;
double*      logk_limits        = NULL;
// double*      linearErrors       = NULL;

// Randoms generation
int          lineNum      = 0;
int          NuRandoms    = 0; 
int          NuQuadrants  = 0;

// Jenkins scaling trick. 
double       Jenkins_foldfactor;

double*      sdltk               = NULL;
double*      sdltPk              = NULL;
double*      sdlt2d              = NULL;

// Camb linear P(k)
double*      lineark             = NULL;
double*      linearPk            = NULL;
double*      linear2d            = NULL;

// Factors for theoretical prediction of redshift space P(k).
double       kaiserFactor;
double       velDispersion;
double       f;
double       beta;
double       y;

double       fsigma8;
double       bsigma8;
double        alpha_pad;
double      epsilon_pad;

const double gamma_GR           =    0.545;
const double gamma_DGP          =  11./16.;

//const double gamma_fR         = 0.41-0.21*z;          

double       Om_mz;

// Clipping
double*     PkCube              = NULL;

// double*     Corrfn              = NULL;
// double*     suppressedCorrfn    = NULL;
// double*     distortedCorrfn     = NULL;
// double*     clippedPk           = NULL;

// double*     W2_veck             = NULL;
// double*     W2_vecr             = NULL;
// double*     FFTW2_vecr_re       = NULL;
// double*     FFTW2_vecr_im       = NULL;

// Binning a 2D, redshift space P(k).
int**        zSpacemodesPerBin  = NULL;
double**     mean_perpk         = NULL;
double**     mean_losk          = NULL;
double**     d2_binnedpk     = NULL;

// Comoving number density calculation.
// double       zBinWidth;

double*      zbins    = NULL;
double*      chibins  = NULL;
double*      Nchi     = NULL;
double*      nbar     = NULL;
double*      nbar_2d  = NULL;

double*      comovVol = NULL;

double chi_interval;
int    chibin_no;

// double*        filteredComovingNumberDensity = NULL;
// double*       MeanSliceRedshift             = NULL;
// double*       LuminosityDistance            = NULL;
// double*       Schechter_fn                  = NULL;

// spline and splint comoving number density, n(z).

// FKP weights. 
double       fkpPk;
double       TotalFKPweight;
double       Interim;

char         surveyType[200];

// clipping threshold
double       appliedClippingThreshold;

// rolling periodic cube. 
double       xcentre;
double       ycentre;
double       zcentre;

double       xroll;
double       yroll;
double       zroll;

double*      inputPk;

// VIPERS ra and dec of cell co-ordinates.

double*   surveyMask;


int CatalogueInput(char filepath[]);

int   CalcCellraDec();

double sdltNz(double z);

int   CatalogNumber;

double steradians2sqdegs(double inSteradians);

int    len;


double*    xdata = NULL;
double*    ydata = NULL;
double*    edata = NULL;
double*  weights = NULL;

double*    zdata = NULL;
int*      Nzdata = NULL;

int   maxiter    = 20000000;


int lineNo;

// gsl random number generation. 
const gsl_rng_type* gsl_ran_T;
gsl_rng*            gsl_ran_r;

int dof;

double sqdegs2steradians(double inSqdegs);

// float* filteredNumberAtRedshift;
// float* filtered_divfuncln_Atz;

// float nzSigma;

//  Apodise the window fn. to supress the Gibb's phenomenon. 
double  GibbsSkinDepth;


double* Cell_SurveyEdge;
double* Cell_ApodiseWeights; 
double* Cell_ShortDist2edge;

double  apodisedVolume;

double  lightconeShot();
double  CubeShot();

int     muBinNumb;
int     modkBinNumb;

double* muBinLimits;

int     fft_size;
int     polar_pkcount;
// int**    polar_modesPerBin;

// double** mean_mu;
// double** mean_modk;
double** polar_pk;
// double** polar2DBinnedPk;


double*  kHexadecapole; 
double*  kQuadrupole;
double*  kMonopole;

double*  kMonopole_expError;
double*  kQuadrupole_expError;

double*  wfMonopole;
double*  wfQuadrupole;

double   TotalW1W4area;                                                                          

double  CentreRA;
double  CentreDec;

double  W1area;
double  W4area;


double  TotalObservedGalaxies    =    0.0;
double  dimmestAcceptedMagnitude =  -99.0;

double xtranslateDist;
double ytranslateDist;

double kaiserGauss_Monofactor(double ks, double beta);
double kaiserGauss_Quadfactor(double ks, double beta);
double kaiserGauss_Hexfactor(double ks, double beta);

double splintLinearPk(double k);

char      theoryPk_flag[200];
char      theoryRSD_flag[200];

double    shotNoiseCorrection_clipfactor;

// translation parameters for Stefano's co-ordinates. 
double    stefano_trans_x;
double    stefano_trans_y;
double    stefano_trans_z;

// mean sampling rate, e.g for FKP calc. 
double    meanSampling;

int min_zShift = 99, max_zShift = 0;
int min_yShift = 99, max_yShift = 0;
int min_xShift = 99, max_xShift = 0;

int polarPk_modeCount = 0;


// Calculation fof the age of the universe, leading to linear growth rate.  Change to double precision.
double*  AgeOftheUniverse;  // Units of H_0
double*  HubbleCnstWithTime;

// Calculation of the linear growth rate. 
void  (*pt2derivs)(float, float[], float[]) = NULL;
void  (*pt2rkqs)(float [], float [], int, float *, float, float, float [], float *, float *, void (*)(float, float [], float[])) = NULL;

int    nVar    = 2;

float* yStartArray;
float* yFinalArray;
float* y2derivsStartArray;

int   linearGrowth_nPoints;

float eps                = 0.00001;
float defaultStepSize    = 0.001;
float MinAllowedStepSize = 0.000001;

float InitialStartTime;
float FinalTime;
float Initial_lnScalefactor;
float Final_lnScalefactor;

    
// Defining declarations for these variables, with memory allocation xp[1..kmax] and yp[1..nvar][1..kmax]
// for the arrays, should be in the calling program.

// float* TimeArray;
double* Age2derivatives;
double* Om_mOfa;
double* f_Om_mOfa545;
double* f_Om_mOfa545_2derivs;
double* approx2linear_growthfactor;

double* lnAarray;
double* lnA2derivatives;
double* linear_growthfactor;
double* SplineParams_ofgrowthfactor;
double* SplineParams_ofgdot;

float* xVals;
float* length_scales;
float* gdot;
float* derivs_error;

float* Second_derivs;
float* Second_derivs_error;
float* Second_deriv_lengthscales;

float* SplineParams_ofgdotdot;

int    nDerivs;
double growthfactor_today;
double approx_growthfactor_today;

// float* HubbleConstantArray;
double* HubbleConstant2derivatives;

// Likelihood calculation.
double**     Multipoles;
double**    dMultipoles;

double*       sigmaNorm;

double  app_sigma8; 

double  ChiSq_kmax;
double  ChiSq_kmin;

// Inclusion of Alcock-Paczynski in Likelihood analysis.
double  min_alpha_pad;
double  max_alpha_pad;

double  min_epsilon_pad;
double  max_epsilon_pad;

double  min_fsigma8;
double  max_fsigma8;

double  min_bsigma8;
double  max_bsigma8; 

double* bsigma8Posterior;
double* fsigma8Posterior;
double*   sigmaPosterior;

double  PosteriorNorm = 0.0;

double    detCovariance;

// int*      ModeNumber;
gsl_matrix* Covariance;
gsl_matrix* sigma_norm;

gsl_vector* eval;
gsl_matrix* evec;

gsl_eigen_symmv_workspace* w;

gsl_vector* col;

double*   MeanMultipoles;

double    A11Sq;

double    linearBias;
double    ChiSqEval();

int       Res;
double    dRes;

double***** ChiSqGrid;
double***** lnLikelihoodGrid;
double    minChiSq;

double    minChiSq_fsigma8;
double    minChiSq_A11Sq;
double    minChiSq_sigma;
double    minChiSq_bsigma8;
double    minChiSq_alpha_pad;
double    minChiSq_epsilon_pad;

double    paramNumber;

double    min_velDisperse;
double    max_velDisperse;

double    min_A11Sq;
double    max_A11Sq;

double    clippedVolume;

// Boot strap catalogue generation.
int*      BootStrap_flag;
double*   BootStrap_Wght;

double*   xBoot;
double*   yBoot;
double*   zBoot;
double*   deltaBoot;

double*  eigenVals;
double** eigenVecs;
    
// double** sigmaNorm;
int order, nrotations, mono_order;

double* xdata;
double* xtheory;
    
double* ydata;
double* ytheory;

double  smallestEigenvalue;

int     chiSq_kmaxIndex;
int     chiSq_kminIndex;

// Index at which 
int     jenkins_foldIndex_foldedfile;
int     jenkins_foldIndex_unfoldedfile;
double  jenkins_fold_kjoin;

int     hiMultipoleOrder;
int     lineNo;

float   (*pt2_pkIntegrand)(float);

int     besseltransform_order;
float   q0;

double Pk_powerlaw(double k, double r0, double gamma);

double* mono;
double* quad;
double*  hex;

double* monop;
double* quadp;
double*  hexp;

double*  cmono;
double* cmonop;

float* rVals;
float* fcmono;
float* fcmono2D;

int    FFTlogRes;

double Pk_powerlaw_truncated_xi(double r);

double HODPk_Gaussian(double k);

// double*  theta;
// double** mixingmatrix;

// double**  rmodulus_vec;
// double*   kmodulus_vec;
// double*         mu_vec;

double LegendrePolynomials(double x, int n);

// double unitTheory(double r, double k, double d, int order);

double u0, variance;

double splintConvQuad(double k);
double splintConvMono(double k);

double app_mean;

double* windfn_rvals;
double* windfn_rMonopole;
double* windfn_rQuadrupole;
double* windfn_rHexadecapole;

double* windfn_rMonopole2d;
double* windfn_rQuadrupole2d;
double* windfn_rHexadecapole2d;

int     rbinNumb;

double* monoCorr;
double* monoCorr2d;
double* rmonoCorr;

double  AbelIntegrand(double);

double (*pt2funcforAbel)(double) = NULL;

double Abel_rp;

double volAvgOfXi;

double volAverage_Corrfn(double rmax);
double volNormedCorrfn(double r);

double  dra;
double  ddec;

int rabin;
int decbin;

int radecbinNumb = 300;

double rabinInterval;
double decbinInterval;

double*  ra_bins;
double* dec_bins;

double** binned_pairs;

double* fiber_flag;
double remaining;

double* rvals;
double*    xi;
double*  xi2D;

int splinexi_N =0;

int*    shuffle_rows;
double* rand_redshift;

double* xDisplacement;
double* yDisplacement;
double* zDisplacement;

int startint = -99;
int endint =-99;

double* sigma2rvals;

double* sigma2perp;
double* sigma2para;

double* sigma2perp_2d;
double* sigma2para_2d;

double zel_vecr[3];
double zel_vecx[3];

double zel_qmin[3];
double zel_qmax[3];

double dummy;
double* kVals;

double** invA;

double zel_r;
double zel_xi;
double zel_err;
double zel_scale;

double splint_sigma2perp(double r);
double splint_sigma2para(double r);

double ukm_nfw_profile(double k);

double haloModel_pk(double k, double nbar, double beta, int Order);
double ukm2(double k);
double linearPk_Gaussian(double k);

double nfw_rvir;
double nfw_conc;

double* r_nfwinversion;
double* q_nfwinversion;
double* r2D_nfwinversion;

double rhobox, Mbox, Mhalo, Mpart;
double Delta_crit = 200.;
double Dplus      = 0.5;
double fGR        = 0.7;

double* hod_disp;

fftw_complex *outx, *outy, *outz;

fftw_plan    iplan_x, iplan_y, iplan_z;

int*     shuffle_rows;

double anisoGauss_delta;
double anisoGauss_asigma;
double anisoGauss_bsigma;


int     VIPERS_mask_lineNo_lo;

double* VIPERS_maskr_lo;

double* VIPERS_maskMono_lo;
double* VIPERS_maskQuad_lo;
double* VIPERS_maskHex_lo;
double* VIPERS_maskOct_lo;
double* VIPERS_maskDec_lo;

double* VIPERS_maskMono2D_lo;
double* VIPERS_maskQuad2D_lo;
double* VIPERS_maskHex2D_lo;
double* VIPERS_maskOct2D_lo;
double* VIPERS_maskDec2D_lo;

double mask_monopolenorm_lo;

int     VIPERS_mask_lineNo_hi;
double* VIPERS_maskr_hi;
double* VIPERS_maskMono_hi;
double* VIPERS_maskQuad_hi;
double* VIPERS_maskHex_hi;
double* VIPERS_maskOct_hi;
double* VIPERS_maskDec_hi;

double* VIPERS_maskMono2D_hi;
double* VIPERS_maskQuad2D_hi;
double* VIPERS_maskHex2D_hi;
double* VIPERS_maskOct2D_hi;
double* VIPERS_maskDec2D_hi;

double mask_monopolenorm_hi;

double loRes_highRes_join;

int     VIPERS_mask_lineNo_hihi;
double* VIPERS_maskr_hihi;
double* VIPERS_maskMono_hihi;
double* VIPERS_maskQuad_hihi;
double* VIPERS_maskHex_hihi;
double* VIPERS_maskOct_hihi;
double* VIPERS_maskDec_hihi;

double* VIPERS_maskMono2D_hihi;
double* VIPERS_maskQuad2D_hihi;
double* VIPERS_maskHex2D_hihi;
double* VIPERS_maskOct2D_hihi;
double* VIPERS_maskDec2D_hihi;

double mask_monopolenorm_hihi;

double loRes_highRes_join;
double hiRes_hihiRes_join;



double splint_VIPERS_maskMono(double r);
double splint_VIPERS_maskQuad(double r);
double splint_VIPERS_maskHex(double r);
double splint_VIPERS_maskOct(double r);
double splint_VIPERS_maskDec(double r);

double splint_VIPERS_maskMultipoles(double r, int transformOrder);

int VIPERS_kSpace_multipoles_lineNo;

double* VIPERS_k;
double* VIPERS_kMono;
double* VIPERS_kMono2D;
double* VIPERS_kQuad;
double* VIPERS_kQuad2D;

double splint_VIPERS_kSpaceMono(double k);
double splint_VIPERS_kSpaceQuad(double k);

int*   fftlog_indices;


// mass functions etc.
double sig2_R;

int   sigres;

// for splint of sigma^2 of R.
double*     sig2;
double*    radii;
double*  sig2_2D;

double*  lnm_jenkins;
double*  ln_invSig_jenkins;
double*  ln_invSig_jenkins_2D;

int      pk_lineNo;

double dummy;

double pk_hin;
double pk_lon;
double pk_hiA;
double pk_loA;

double modeCorrections(int k, int j, int i);

double* Rr;
double* Ir;
double* PkCube;

double AcceptedMax(double a[], bool b[], int n);
double AcceptedMin(double a[], bool b[], int n);

double chi_invertedStefanoBasis(double* xval, double* yval, double* zval);

double* expected_foldedRandCounts;
double calc_volavg_fkpweights();

double calc_fsigma8Posterior();

double* gal_clippingweights;

double clipping_fSq;

double calc_bsigma8Posterior();
double calc_velDispPosterior();

double spherical_tophat(double k, double R);

double fraction_clipped;

double underlyingGaussian_sigma;

double GaussianFilter_radius;
double depletion_factor;

double clipping_smoothing_radius;

int    data_mock_flag;

double* spline_lnk;
double* dlnPR_dlnk;
double* dlnPR_dlnk_2D;

double        cumulative_nbar[100];
double     cumulative_nbar_2d[100];
double    chi_cumulative_nbar[100];

double inverse_cumulative_nbar(double arg);
double nz_smoothRadius;

int     Res_ap;
double dRes_ap;

double (*pt2maskMultipoles)(double r, int transformOrder) = NULL;
double splint_unit_maskMultipoles(double r, int transformOrder);
