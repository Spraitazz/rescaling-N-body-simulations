int FirstColumnCompare(const void *pa, const void *pb){
    // pointer to const. int      // This is a derefernce (indirection) operator. 
    const double *a                = *(const double **) pa;
    const double *b                = *(const double **) pb;

    if(a[0] == b[0]){
        if(a[1] - b[1] < 0)   return -1;
        else                  return  1;
    }
    
    else{ 
        if(a[0] - b[0] < 0)   return -1;
        
        else                  return  1;
    }
}

int SecondColumnCompare(const void *pa, const void *pb){
    // pointer to const. int      // This is a derefernce (indirection) operator. 
    const double *a                = *(const double **) pa;
    const double *b                = *(const double **) pb;

    if(a[1] == b[1]){
        if(a[2] - b[2] < 0)   return -1;
        else                  return  1;
    }
    
    else{ 
        if(a[1] - b[1] < 0)   return -1;
        else                  return  1;
    }
}


int TwoColumnCompareTest(){
    double** testArray;

    int length = 40;

    testArray = malloc(length*sizeof(double*));                         // rows

    for(j=0; j<length; j++){ 
        testArray[j]     = malloc(3*sizeof(double));
        testArray[j][0]  = (double) 1./(rand() % 7 + 1.);
        testArray[j][1]  = (double) 1./(rand() % 7 + 1.);
        testArray[j][2]  = (double) 1./(rand() % 7 + 1.);
    
        printf("\n%.6f \t \t %.6f \t \t %.6f", testArray[j][0], testArray[j][1], testArray[j][2]);
    }

    printf("\n\n");

    qsort(testArray, length, sizeof(testArray[0]), FirstColumnCompare);

    printf("\n\n");

    for(j=0; j<length; j++) printf("\n%.6f \t \t %.6f \t \t %.6f", testArray[j][0], testArray[j][1], testArray[j][2]);

    qsort(testArray, length, sizeof(testArray[0]), SecondColumnCompare);

    printf("\n\n");

    for(j=0; j<length; j++) printf("\n%.6f \t \t %.6f \t \t %.6f", testArray[j][0], testArray[j][1], testArray[j][2]);

    for(j=0; j<length; j++) free(testArray[j]);
    free(testArray);
    
    return 0;
}


int prep_grid(){    
    overdensity                = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n0*n1*n2);
  
    // initialise overdensity[0]
    for(j=0; j<n0*n1*n2;   j++) overdensity[j][0] = 0.;
    for(j=0; j<n0*n1*n2;   j++) overdensity[j][1] = 0.;
    
    return 0;
}

int prep_fftw(){    
    H_k                = (fftw_complex*) fftw_malloc(n0*n1*n2*sizeof(fftw_complex)); 
    
    p                  = fftw_plan_dft_3d(n0, n1, n2, overdensity, H_k, FFTW_FORWARD,  FFTW_ESTIMATE);

    return 0;
}

int assign2DPkMemory(){
    twodim_pk           = (double **) realloc(twodim_pk, n0*n1*n2*sizeof(double*));                         
 
    for(j=0; j<n0*n1*n2; j++)  twodim_pk[j]    = (double *)  malloc(3*sizeof(double));
    
    d2_binnedpk         = (double **) realloc(d2_binnedpk, 50*sizeof(double* ));
    
    for(j=0; j<50; j++) d2_binnedpk[j] = (double *) malloc(50*sizeof(double));
    
    for(k=0; k<50; k++) for(j=0; j<50; j++)  d2_binnedpk[k][j] = 0.0;
    
    return 0;
}

int prep_pkRegression(double logk_min, double logk_max, double bin_no){    
    double logk_interval;

    polar_pk             = (double **) malloc(n0*n1*n2*sizeof(double*));
    
    // cols of mod k, mu, pk.
    for(j=0; j<n0*n1*n2; j++)  polar_pk[j]    = (double *)  malloc(3*sizeof(double)); 

    logk_interval        = (logk_max - logk_min)/bin_no;
 
    logk_limits          = (double *) malloc(bin_no*sizeof(*logk_limits));
    
    for(j=0; j<bin_no; j++) {    
    	logk_limits[j] = pow(10., logk_min + j*logk_interval);
    }

    mean_modk            = (double *) malloc((bin_no-1)*sizeof(*mean_modk));
    
    modes_perbin         = (int *)    malloc((bin_no-1)*sizeof(*modes_perbin));
    
    binnedPk             = (double *) malloc((bin_no-1)*sizeof(*binnedPk));
    
    // linearErrors      = (double *) malloc((kBinNumb-1)*sizeof(*linearErrors));
    
    kMonopole            = (double  *) malloc((bin_no-1)*sizeof(*kMonopole));
    kQuadrupole          = (double  *) malloc((bin_no-1)*sizeof(*kQuadrupole));
    kHexadecapole        = (double  *) malloc((bin_no-1)*sizeof(*kHexadecapole));
    
    return 0;
}

int free_pkRegression(){
  fftw_free(H_k);
  
  fftw_destroy_plan(p);
    
  for(j=0; j<n0*n1*n2; j++)  free(polar_pk[j]);
    
  free(polar_pk);
    
  // binned Pk arrays.                                                                                                   
  free(logk_limits); 
  free(mean_modk);
  free(binnedPk);
  free(modes_perbin);
  
  free(kMonopole);
  free(kQuadrupole);
  free(kHexadecapole);
  
  return 0;
}

int PkCalc(){  
    fftw_execute(p);
	printf("generated fourier overdensity \n");
    assign2DPkMemory();
    printf("assigned 2dpk memory \n");
    
    PkCorrections();
    
    //observedQuadrupole(polar_pkcount);
    
    return 0;
}

int PkCorrections(){
    double pk, WindowFunc, k_x, k_y, k_z;
    double mu, kmodulus, kSq;
    int Index;

    double rand_shot = 0.0, gal_shot = 0.0;

	// GRID   -> no shot noise
	// RANDOM -> 1/<n>
	
	double cellSizes[3], fundmodes[3], nyquist[3];
	double xNyquistWaveNumber, yNyquistWaveNumber, zNyquistWaveNumber;
	double kIntervalx, kIntervaly, kIntervalz;
	
	cellSizes[0] = volume_limits[0] / (double) cells[0];
	cellSizes[1] = volume_limits[1] / (double) cells[1];		
	cellSizes[2] = volume_limits[2] / (double) cells[2];

    fundmodes[0] = 2.0*pi/volume_limits[0];
    fundmodes[1] = 2.0*pi/volume_limits[1];
    fundmodes[2] = 2.0*pi/volume_limits[2];
    
    kIntervalx = fundmodes[0];
    kIntervaly = fundmodes[1];
    kIntervalz = fundmodes[2];
    
    nyquist[0] = pi / cellSizes[0];
    nyquist[1] = pi / cellSizes[1];
    nyquist[2] = pi / cellSizes[2];
    
    xNyquistWaveNumber = nyquist[0];
    yNyquistWaveNumber = nyquist[1];
    zNyquistWaveNumber = nyquist[2];
    
    double logkmin = log10(kIntervalx * sqrt(3.0));
    double logkmax = log10(0.8);
    
    prep_pkRegression(logkmin, logkmax, spectrum_size);

    polar_pkcount = 0;
    printf("starting P(k) loops \n");
    for(k=0; k<n0; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
                k_x = kIntervalx*i;
                k_y = kIntervaly*j;
                k_z = kIntervalz*k;
                
                if(k_x>xNyquistWaveNumber)  k_x   -= n2*kIntervalx;
                if(k_y>yNyquistWaveNumber)  k_y   -= n1*kIntervaly;
                if(k_z>zNyquistWaveNumber)  k_z   -= n0*kIntervalz;

                Index                              = k*n1*n2 + j*n2 + i;

                kSq                                = pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);
                
                kmodulus                           = pow(kSq, 0.5);
                
                mu                                 = k_z/kmodulus;
                if(kmodulus < 0.000001)       mu   = 0.0;     

                WindowFunc                         = 1.;

                if(k_x != 0.)  WindowFunc         *= sin(pi*k_x*0.5/xNyquistWaveNumber)/(pi*k_x*0.5/xNyquistWaveNumber);

                if(k_y != 0.)  WindowFunc         *= sin(pi*k_y*0.5/yNyquistWaveNumber)/(pi*k_y*0.5/yNyquistWaveNumber);
    
                if(k_z != 0.)  WindowFunc         *= sin(pi*k_z*0.5/zNyquistWaveNumber)/(pi*k_z*0.5/zNyquistWaveNumber);
                
                // Correct for mass assignment of randoms.
                // If smoothing the field, do not correct for mass assignment.  Obvious in hindsight.    
                // Cloud in cell. NGP corresponds to WindowFunc rather than WindowFunc^2. 
                
                
				H_k[Index][0] /= ((double)cells[0]*cells[1]*cells[2]);
				H_k[Index][1] /= ((double)cells[0]*cells[1]*cells[2]);
				
				//H_k[Index][0] /= volume;
				//H_k[Index][1] /= volume;
                
                //WHAT IS THIS????????????????
               // H2_k[Index][0] /= pow(WindowFunc, 1.);
               // H2_k[Index][1] /= pow(WindowFunc, 1.);
                
                                                     
               // H_k[Index][0] = H2_k[Index][0]; 
		       // H_k[Index][1] = H2_k[Index][1]; 
		       
		       //printf("b4 hk, [%d][0]: %lf, [%d][1]: %lf \n", Index, H_k[Index][0], Index, H_k[Index][1]);
		       H_k[Index][0] /= pow(WindowFunc, 1.0);
		       H_k[Index][1] /= pow(WindowFunc, 1.0);
            	//printf("after hk, pk count: %d \n", polar_pkcount);
                pk = pow(H_k[Index][0], 2.) + pow(H_k[Index][1], 2.);

				// shot noise correction.
                //pk -= volume / particle_no;
                pk *= volume;
                
	            if(kmodulus > 0.000001){
	                //// Only half the modes are independent. ////
	            	if(k_z>0.){
	            	    // One hemi-sphere is independent, e.g. k_z >= 0.
		                polar_pk[polar_pkcount][0] = kmodulus;
		                polar_pk[polar_pkcount][1] = fabs(mu);
		                polar_pk[polar_pkcount][2] = pk;		            	   		            	            		            		           
		            
		                polar_pkcount               += 1;
		            }
		            
		            else if((k_z == 0.0) && (k_y > 0.0)){
                        // in the k_z=0 plane one semi-circle is independent, k_y>0.
		                polar_pk[polar_pkcount][0]   = kmodulus;
		                polar_pk[polar_pkcount][1]   = fabs(mu);
		                polar_pk[polar_pkcount][2]   = pk;
		            
		                polar_pkcount                += 1;
		            }
		            
		            else if((k_z == 0.0) && (k_y == 0.0) && (k_x > 0.0)){
		                // on the line k_z=k_y=0, one half is independent, k_x>=0.
		                                        // in the k_z=0 plane one semi-circle is independent, k_y>0.
		                polar_pk[polar_pkcount][0]    = kmodulus;
		                polar_pk[polar_pkcount][1]    = fabs(mu);
		                polar_pk[polar_pkcount][2]    = pk;

		                polar_pkcount                += 1;
		            }
	            }
	        }
        }
    }

	// output filepath.  
	char filepath[100];  
    sprintf(filepath,"%s/output/Pk_Mikes_test.dat", home_directory);
    
    //MultipoleCalc(kbin_no, mean_modk, kMonopole, kQuadrupole, polar_pk, modeCount, filepath, 0.0, 1.0, 1);
    MonopoleCalc(spectrum_size, mean_modk, kMonopole, polar_pk, polar_pkcount, outFile_Pk_Mikes, 0.0, 1.0, 1);

    return 0;
}

int MonopoleCalc(int modBinNumb, double mean_modBin[], double Monopole[], double** Array, int modeCount, char filepath[], double mu_lolimit, double mu_hilimit, int fileOutput){
    printf("\n\nPerforming multipole calculation. (to Monopole order)\n");
    
    int          loIndex;
    int          hiIndex;
    int     modes_perbin;

    qsort(Array, modeCount, sizeof(Array[0]), FirstColumnCompare);

    loIndex = 0;
    hiIndex = 0;
    
    for(j=0; j<modBinNumb-1; j++){
        mean_modBin[j]    = 0.0;
        
        Monopole[j]       = 0.0;
    }
    
    for(j=0; j<modeCount; j++){
        if(Array[j][0]   >= logk_limits[0]){
            loIndex = j; 
            
            break;
        }
    }
   
    if(fileOutput==1)  output = fopen(filepath, "w");
  
    
    for(j=0; j<modBinNumb-1; j++){
        modes_perbin = 0;
    
        for(i=loIndex; i<modeCount; i++){
        	//printf("array[i][0]: %lf \n", Array[i][0]);
            if(Array[i][0] > logk_limits[j+1]){
                hiIndex = i;
                
                break;
            } 
        }
        // printf("lo index: %d, hi index: %d \n", loIndex, hiIndex);    
        
        for(i=loIndex; i<hiIndex; i++){
            if((mu_lolimit<Array[i][1]) && (Array[i][1]<mu_hilimit)){      
                mean_modBin[j] += Array[i][0];
                   
                   Monopole[j] += Array[i][2];
                
                modes_perbin   += 1;
                
               // printf("mean modbin: %lf, monopole: %lf \n", mean_modBin[j], Monopole[j]);
            }
        }
        	    
        if(modes_perbin != 0)  mean_modBin[j]  /= modes_perbin;
        if(modes_perbin != 0)     Monopole[j]  /= modes_perbin;

        // Peacock and Nicholson 1991, pg 313. above eqn (20).
        // Result of summing over a shell in k space containing m modes, should be a Gaussian random variable with variance 2.m/N^2  
        
        if(fileOutput==1 && modes_perbin != 0)  {        
        	fprintf(output, "%e \t %e \t %d \n", mean_modBin[j], Monopole[j], modes_perbin);       
        }
        
        loIndex       = hiIndex;
    }
    
    if(fileOutput==1)  fclose(output);
    
    return 0;
}

int MultipoleCalc(int modBinNumb, double mean_modBin[], double Monopole[], double Quadrupole[], double** Array, int modeCount, char filepath[], double mu_lolimit, double mu_hilimit, int fileOutput){
    // P(k, mu_i) = P_0(k) + P_2(k)*(3.*mu_i*mu_i - 1.)/2
    // Least squares fit between theory prediction, P(k, mu_i) and measured.  (Linear regression).

    int          loIndex;
    int          hiIndex;
    int     modes_perbin;

    printf("\n\nPerforming multipole calculation. (to Quadrupole order)");
    
    // assign log k binning for pk and assign memory for binning. (logk_min, logk_max, # bins).
    // prep_pkbinning(-2., log10(modkMax), kbin_no);
    
    // Order by mod k to ensure binning is the mean between LowerBinIndex and UpperBinIndex.
    qsort(Array, modeCount, sizeof(Array[0]), FirstColumnCompare);
    
    loIndex      = 0;
    hiIndex      = 0;

    for(k=0; k<modBinNumb-1; k++){
        mean_modBin[k]    = 0.0;
        Monopole[k]       = 0.0;
        Quadrupole[k]     = 0.0;
    }
    
    for(i=0; i<modeCount; i++){
        if(Array[i][0] >= logk_limits[0]){
            loIndex = i; 
            break;
        }
    }
    
    if(fileOutput==1)  output = fopen(filepath, "w");
    
    for(j=0; j<modBinNumb-1; j++){
        //  Linear regression against P_i(k, \mu_i)  
        //  L_i  = 0.5*(3.*\mu_i**2 -1.)
        //  P_i  = P_i(k, \mu_i)
    
        double Li          = 0.0;
        double Pi          = 0.0;
    
        double Sum_Li      = 0.0;
        double Sum_Li2     = 0.0;
    
        double Sum_Pi      = 0.0;
        double Sum_PiLi    = 0.0;
        
        modes_perbin = 0;
        
        // Find the range of indices corresponding to the modes in a given k interval. 
        for(i=loIndex; i<modeCount; i++){
            if(Array[i][0] > logk_limits[j+1]){
                hiIndex = i;
                break;
            } 
        }
        
        for(i=loIndex; i<hiIndex; i++){
            if((mu_lolimit<Array[i][1]) && (Array[i][1]<mu_hilimit)){        
                mean_modBin[j] += Array[i][0];
                modes_perbin   += 1;
            
                         // L_i = 0.5*(3.*mu**2 -1.)
                     
                Li              = 0.5*(3.*pow(Array[i][1], 2.) - 1.);     
                Pi              = Array[i][2];
            
                if(j==2)  printf("\n%e \t %e \t %e \t %e", Array[i][0], Array[i][1], Array[i][2], Li);
            
                Sum_Li         += Li; 
                Sum_Li2        += Li*Li;
            
                Sum_Pi         += Pi;
                Sum_PiLi       += Pi*Li;
            }
         }
    
         mean_modBin[j]        /= modes_perbin;
        	    
         // For a matrix A, paramater vector, (P_0, P_2)^T, P and vector B
         // Required to invert AP  = B. 2x2 matrix inversion.
         	
         // A reads as (a b) for a = sum_Modes 1, b = sum_Modes L_i, c = sum_Modes L_i, d = sum_Modes L_i**2.
         //            (c d)        

         // and B = (b1, b2)^T for b1 = sum_Modes measured Pi, b2 = sum_Modes (measured Pi)*Li
     
         double detA;
         // det = ad - bc.
     
         detA                       = modes_perbin*Sum_Li2 - Sum_Li*Sum_Li;
     
         if(detA<pow(10., -10.))      printf("\n\nCannot decompose into monopole and quadrupole, detA = 0.0");
     
         if(j==2)  printf("\n\ndet A: %e", detA); 
     
         // (P_0, P_2)^T = (A^-1)B  = (1./detA)*(d*b1 - b*b2, -c*b1 + a*b2)^T.
     
                     Monopole[j]    = (1./detA)*( Sum_Li2*Sum_Pi - Sum_Li*Sum_PiLi);
                   Quadrupole[j]    = (1./detA)*(-Sum_Li*Sum_Pi  + modes_perbin*Sum_PiLi);
                   
         // for(jj=0; jj<= j; jj++)  printf("\n%d \t %e \t %e \t %e \t %d", j, mean_modBin[jj], Monopole[jj], Quadrupole[jj], modesperbin[jj]);
         
         // printf("\n\n");
         
         
         //for(k=LowerBinIndex; k<UpperBinIndex; k++){
           //if((mu_lolimit<Array[k][1]) && (Array[k][1]<mu_hilimit)){
             //Li                       = 0.5*(3.*pow(Array[k][1], 2.) - 1.);
         
             //kMonopole_expError[j]   += pow((*pt2Pk)(meanKBin[j])*pow(1. + beta*pow(Array[k][1], 2.), 2.)*pow(1. + 0.5*pow(meanKBin[j]*velDispersion*Array[k][1], 2.), -1.)*(Sum_Li2 - Li*Sum_Li), 2.);
           
             //kQuadrupole_expError[j] += pow((*pt2Pk)(meanKBin[j])*pow(1. + beta*pow(Array[k][1], 2.), 2.)*pow(1. + 0.5*pow(meanKBin[j]*velDispersion*Array[k][1], 2.), -1.)*(-1.*Sum_Li + + modesperbin[j]*Li), 2.);
           //}
         //}

           //kMonopole_expError[j]   *= pow(detA, -2.); 
         //kQuadrupole_expError[j]   *= pow(detA, -2.);
         
         loIndex   = hiIndex;
         
         if((fileOutput==1) && (detA>pow(10., -6.)))  fprintf(output, "%e \t %e \t %e \t %d \n", mean_modBin[j], Monopole[j], Quadrupole[j], modes_perbin);
     } 
    
     if(fileOutput==1) fclose(output);
       
     // free logk_limits, mean_modk, binnedPk, modes_perbin
     // free_pkRegression();
     
     return 0;
}


int DualBinning(int NumberModes, double** DualParamArray, double** BinnedDualParamArray){
    int m;

    int bin_no = 50;

    double firstBinLimits[bin_no];
    double secndBinLimits[bin_no];

    int modesPerBin[bin_no][bin_no];

    int firstColLowerBinIndex     = 0;
    int firstColUpperBinIndex     = 0;

    int secndColLowerBinIndex     = 0;
    int secndColUpperBinIndex     = 0;
    

    qsort(DualParamArray, NumberModes, sizeof(DualParamArray[0]), FirstColumnCompare);

    // Order by first then second column.
    printf("\nDual param array sorted.");
    
    for(j=0; j<bin_no; j++) firstBinLimits[j] = 0.01 + j*1./bin_no;
    for(j=0; j<bin_no; j++) secndBinLimits[j] = 0.01 + j*1./bin_no;


    // for(j=0; j<20; j++)  printf("\n%.3lf \t %.3lf \t %.3lf", DualParamArray[j][0], DualParamArray[j][1], DualParamArray[j][2]);
    
    for(j=0; j<bin_no; j++){
        for(k=0; k<bin_no; k++){
            BinnedDualParamArray[j][k] = 0.0;
            // mean_firstCol[j][k]        = 0.0;
            // mean_secndCol[j][k]        = 0.0;
            modesPerBin[j][k]          =   0;    
        }
    }
    
    for(j=0; j<NumberModes; j++){
        if(DualParamArray[j][0]   >= firstBinLimits[0]){
            firstColLowerBinIndex = j; 
        
            break;
        }
    }
    
    for(j=0; j<bin_no; j++){
        for(i=firstColLowerBinIndex; i<NumberModes; i++){
            if(DualParamArray[i][0] > firstBinLimits[j+1]){
                firstColUpperBinIndex = i;
                break;
            } 
        }
        
        secndColLowerBinIndex = firstColLowerBinIndex;
    
	    qsort(&DualParamArray[firstColLowerBinIndex], firstColUpperBinIndex - firstColLowerBinIndex, sizeof(DualParamArray[0]), SecondColumnCompare);
    
        // meets the lower bin limit
        for(i=secndColLowerBinIndex; i<NumberModes; i++){
            if(DualParamArray[i][1]   >= secndBinLimits[0]){
                secndColLowerBinIndex = i; 
        
                break;
            }
        }      
       
        // meets the upper bin limit. 
        for(k=0; k<bin_no; k++){
            for(m=secndColLowerBinIndex; m<firstColUpperBinIndex; m++){
                if(DualParamArray[m][1] > secndBinLimits[k+1]){
                    secndColUpperBinIndex = m;
                    break;
                } 
            }
        
            for(m=secndColLowerBinIndex; m<secndColUpperBinIndex; m++){
                BinnedDualParamArray[j][k]    += DualParamArray[m][2];
                // mean_firstCol[j][k]           += DualParamArray[m][0];
                // mean_secndCol[j][k]           += DualParamArray[m][1];
                modesPerBin[j][k]             += 1;    
            }
            
            if(modesPerBin[j][k] != 0)  BinnedDualParamArray[j][k] /= modesPerBin[j][k];
            // if(modesPerBin[j][k] != 0)  mean_firstCol[j][k]        /= modesPerBin[j][k];
            // if(modesPerBin[j][k] != 0)  mean_secndCol[j][k]        /= modesPerBin[j][k];
            
            secndColLowerBinIndex = secndColUpperBinIndex;
        }
        
        firstColLowerBinIndex = firstColUpperBinIndex;
    }
    
    return 0;
}
