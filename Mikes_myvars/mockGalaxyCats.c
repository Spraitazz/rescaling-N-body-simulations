int prep_DisplacementCalc(){
    xDisplacement       = malloc(n0*n1*n2*sizeof(double));
    yDisplacement       = malloc(n0*n1*n2*sizeof(double));
    zDisplacement       = malloc(n0*n1*n2*sizeof(double));
    
    outx                = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n0*n1*n2);  
    outy                = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n0*n1*n2);
    outz                = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n0*n1*n2);
    
    H_k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n0*n1*n2);

    iplan_x 			= fftw_plan_dft_3d(n0, n1, n2, outx, H_k, FFTW_BACKWARD, FFTW_ESTIMATE);
    iplan_y 			= fftw_plan_dft_3d(n0, n1, n2, outy, H_k, FFTW_BACKWARD, FFTW_ESTIMATE);
    iplan_z 			= fftw_plan_dft_3d(n0, n1, n2, outz, H_k, FFTW_BACKWARD, FFTW_ESTIMATE);

    return 0;
}


int DisplacementCalc(){
  int m0, m1, m2;

  int negkIndex;

  double Power, amplitude, phase, expectation;
  double k_x, k_y, k_z, kSq, kmodulus, mu, kIntervalx, kIntervaly, kIntervalz;

  int Index;
  
   kIntervalx = 2.0*pi/volume_limits[0];
   kIntervaly = 2.0*pi/volume_limits[1];
   kIntervalz = 2.0*pi/volume_limits[2];

  printf("\nGenerating displacements.\n");

  for(int k=0; k<n0; k++){  
      for(int j=0; j<n1; j++){
          for(int i=0; i<n2; i++){
	          m0 = k;
	          m1 = j;
	          m2 = i;

	          if(m2>n2/2)  m2                   -= n2;
	          if(m1>n1/2)  m1                   -= n1;
	          if(m0>n0/2)  m0                   -= n0;

          	  k_x                                = kIntervalx*m2;
	          k_y                                = kIntervaly*m1;
	          k_z                                = kIntervalz*m0;

	          Index                              = k*n1*n2 + j*n2 + i;

	          kSq                                = pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);

	          kmodulus                           = pow(kSq, 0.5);

	          //mu                                 = k_z/kmodulus;
	         // if(kmodulus < 0.000001)       mu   = 0.0;

          	  //expectation                        = (*pt2Pk)(kmodulus)/TotalVolume;
          	  expectation = splint_generic_VPk(&Pk_current, kmodulus);

	          // expectation                    *= pow(1. + beta*pow(mu, 2.), 2.);

	          // expectation                    /= 1. + 0.5*pow(kmodulus*mu*velDispersion, 2.);                                                            

	          // Power                          = -log(gsl_rng_uniform(gsl_ran_r))*expectation;

	          // amplitude                       = sqrt(Power);                                                                                            
	          amplitude                          = sqrt(expectation);
	          
				

	          phase                              = 2.*pi*gsl_rng_uniform(gsl_ran_r);
	          //if (Index < 1000) {
	          	//	printf("k: %lf, amplitude: %lf, phase: %lf \n", k, amplitude, phase);
	          //	}
	          

	          outx[Index][0]                     =  (k_x/kSq)*amplitude*sin(phase);
	          outx[Index][1]                     = -(k_x/kSq)*amplitude*cos(phase);
	          
	          //if (Index < 1000) {
	          	//printf("amplitude: %lf, phase: %lf \n", amplitude, phase);
	          	//printf("out x real: %lf, out x im: %lf \n", outx[Index][0], outx[Index][1]);
	        //  }
      
	          outy[Index][0]                     =  (k_y/kSq)*amplitude*sin(phase);
              outy[Index][1]                     = -(k_y/kSq)*amplitude*cos(phase);

	          outz[Index][0]                     =  (k_z/kSq)*amplitude*sin(phase);
              outz[Index][1]                     = -(k_z/kSq)*amplitude*cos(phase);
        }
      }
  }

  outx[0][0] = outx[0][1] = 0.0;
  outy[0][0] = outy[0][1] = 0.0;
  outz[0][0] = outz[0][1] = 0.0;
  
  // Hermitian condition. One hemi-sphere is independent, e.g. k_z >= 0.
  for(int k=n0-1; k>=n0/2; k--){
    for(int j=0; j<n1; j++){
        for(int i=0; i<n2; i++){
            negkIndex          = k*n1*n2 + j*n2 + i;
                       
            Index              = 0;
                
            // zero maps to zero on reflection through the origin.
            if(i!=0)  Index   += (n2 - i);
            if(j!=0)  Index   += (n1 - j)*n2;
                      Index   += (n0 - k)*n1*n2;

            outx[negkIndex][0]  =     outx[Index][0];
            outx[negkIndex][1]  = -1.*outx[Index][1];

            if(negkIndex == Index)    outx[Index][1] = 0.0; // purely real
            
            outy[negkIndex][0]  =     outy[Index][0];
            outy[negkIndex][1]  = -1.*outy[Index][1];

            if(negkIndex == Index)    outy[Index][1] = 0.0; // purely real
            
            outz[negkIndex][0]  =     outz[Index][0];
            outz[negkIndex][1]  = -1.*outz[Index][1];

            if(negkIndex == Index)    outz[Index][1] = 0.0; // purely real
        }
      }
    } 
    
    // in the k_z=0 plane one semi-circle is independent, k_y>0.         
    for(int j=n1-1; j>=n1/2; j--){
        for(int i=0; i<n2; i++){
            negkIndex          = j*n2 + i;
                       
            Index              = 0;
                
            // zero maps to zero on reflection through the origin.
            if(i!=0)  Index   += (n2 - i);
                      Index   += (n1 - j)*n2;

            outx[negkIndex][0]  =     outx[Index][0];
            outx[negkIndex][1]  = -1.*outx[Index][1];
            if(negkIndex == Index)    outx[Index][1] = 0.0;
            
            outy[negkIndex][0]  =     outy[Index][0];
            outy[negkIndex][1]  = -1.*outy[Index][1];
            if(negkIndex == Index)    outy[Index][1] = 0.0;
            
            outz[negkIndex][0]  =     outz[Index][0];
            outz[negkIndex][1]  = -1.*outz[Index][1];
            if(negkIndex == Index)    outz[Index][1] = 0.0;
        }
    }
    
    // on the line k_z=k_y=0, one half is independent, k_x>=0.
    for(int i=n2-1; i>=n2/2; i--){
        negkIndex          = i;
                       
        Index              = 0;
                
        // zero maps to zero on reflection through the origin.
        Index             += (n2 - i);

        outx[negkIndex][0]  =      outx[Index][0];
        outx[negkIndex][1]  =  -1.*outx[Index][1];
        if(negkIndex == Index)     outx[Index][1] = 0.0;
        
        outy[negkIndex][0]  =      outy[Index][0];
        outy[negkIndex][1]  =  -1.*outy[Index][1];
        if(negkIndex == Index)     outy[Index][1] = 0.0;
        
        outz[negkIndex][0]  =      outz[Index][0];
        outz[negkIndex][1]  =  -1.*outz[Index][1];
        if(negkIndex == Index)     outz[Index][1] = 0.0;
    }
    
    printf("starting ffts for displacement \n");

    fftw_execute(iplan_x);
    for(int j=0; j<n0*n2*n1; j++) {    
    	xDisplacement[j] = H_k[j][0];
    	if (H_k[j][1] > 1e-9) {
    		printf("imaginary displacement nonzero: %le \n", H_k[j][1]);
    	}
    }

    fftw_execute(iplan_y);
    for(int j=0; j<n0*n2*n1; j++)  yDisplacement[j] = H_k[j][0];
  
    fftw_execute(iplan_z);
    for(int j=0; j<n0*n2*n1; j++)  zDisplacement[j] = H_k[j][0];
    
   printf("done ffts \n");
        
    // fftw_destroy_plan(iplan_x);
    // fftw_destroy_plan(iplan_y);
    // fftw_destroy_plan(iplan_z);
    
    return 0;
}

/*
int periodic_bounds(double* x, int n){
  if(*x>AxisLimsArray[1][n]) *x -= (AxisLimsArray[1][n] - AxisLimsArray[0][n]);
  if(*x<AxisLimsArray[0][n]) *x += (AxisLimsArray[1][n] - AxisLimsArray[0][n]);			      

  return 0;			    
}
*/

/*
int poissonSample_lnNorm(){
    int GalNumber;
    
    double x, y, z, dx, dy, dz, nbar, lambda;
    
    nbar = pow(10., 6.)/TotalVolume;
    
    sprintf(filepath, "%s/Data/lnnorm_anisotropic/poissonSampled_lnNormal_fogOnly_500_particles_1000000_%d.dat", root_dir, loopCount);

    output = fopen(filepath, "w");    	 
    
    for(k=0; k<n0; k++){      
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){    
	            x  = AxisLimsArray[0][2] + i*CellSize;
	            y  = AxisLimsArray[0][1] + j*CellSize; 
	            z  = AxisLimsArray[0][0] + k*CellSize;

    	        Index = k*n1*n2 + j*n2 + i;
    
      	        // periodic_bounds(&x, 2);
	            // periodic_bounds(&y, 1);
	            // periodic_bounds(&z, 0);
    
                lambda    = CellVolume*nbar*(1. + densityArray[Index]);

                GalNumber = gsl_ran_poisson(gsl_ran_r, lambda);
    
                for(ii=0; ii<GalNumber; ii++){  
                    dx = CellSize*gsl_rng_uniform(gsl_ran_r);
	                dy = CellSize*gsl_rng_uniform(gsl_ran_r);
	                dz = CellSize*gsl_rng_uniform(gsl_ran_r);
	                
                    fprintf(output, "%e \t %e \t %e \n", x + dx, y + dy, z + dz);
                }
           }
        }
   }
   
   fclose(output);

   return 0;
}


int poissonSample_homogeneous(double maxGals){
    printf("\nPoisson sampling homogeneous.");

    double x, y, z;
    double GalNumber = 0.;

    sprintf(filepath, "%s/Data/stacpolly/poissonSampled_randoms_500_%.2e.dat", root_dir, maxGals); 

    output = fopen(filepath, "w");    	
    
    while(GalNumber<maxGals){
        x = (AxisLimsArray[1][2] - AxisLimsArray[0][2])*gsl_rng_uniform(gsl_ran_r) + AxisLimsArray[0][2];
        y = (AxisLimsArray[1][1] - AxisLimsArray[0][1])*gsl_rng_uniform(gsl_ran_r) + AxisLimsArray[0][1];
        z = (AxisLimsArray[1][0] - AxisLimsArray[0][0])*gsl_rng_uniform(gsl_ran_r) + AxisLimsArray[0][0];

        fprintf(output, "%e \t %e \t %e \n", x, y, z);

        GalNumber += 1.0;
   }
   
   fclose(output);

   return 0;
}


int load_homogeneous_rands(double maxGals, int load){
    sprintf(filepath, "%s/Data/stacpolly/poissonSampled_randoms_500_%.2e.dat", root_dir, maxGals); 

    inputfile   = fopen(filepath, "r");

    ch          = 0;
    rand_number = 0;

    do{
        ch = fgetc(inputfile);
        if(ch == '\n')
            rand_number += 1;
    } while (ch != EOF);

    printf("\n\n%d randoms.", rand_number);

    if(load == 1){
        rewind(inputfile);
    
        assign_randmemory();

        for(j=0; j<rand_number; j++)  fscanf(inputfile, "%le \t %le \t %le", &rand_x[j], &rand_y[j], &rand_z[j]);
    }
    
    fclose(inputfile);
    
    
    // randomiseCatalogue(rand_number, rand_x, rand_y, rand_z);
    
    // printf("\n\n");

    // for(j=0; j<rand_number; j++)  printf("\n%e \t %e \t %e", rand_x[j], rand_y[j], rand_z[j]);
    
    return 0;
}


int write_homogeneous_rands_window_gridded(int load, double sampling){
    // never load already smoothed randoms. 
    load_homogeneous_rands_window(load, sampling);
    
    int  loopCount      = 0;
    int  rand_count     = 0;
    int  nonempty_cells = 0;
    
    int* nonempty_indices;
    int* nonempty_cellx;
    int* nonempty_celly;
    int* nonempty_cellz;
    
    double x,y,z;
    
    for(j=0; j<n0*n1*n2; j++){
        if(Cell_SurveyLimitsMask[j] > 0.0){
            nonempty_cells += 1;
        }
    }
    
    printf("\nNumber of non-empty cells: %d", nonempty_cells);
    
    
    nonempty_indices = malloc(nonempty_cells*sizeof(int));
    
    nonempty_cellx   = malloc(nonempty_cells*sizeof(double));
    nonempty_celly   = malloc(nonempty_cells*sizeof(double));
    nonempty_cellz   = malloc(nonempty_cells*sizeof(double));
    
    for(k=0; k<n0; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
                Index = i + n2*j + n2*n1*k;
            
                if(Cell_SurveyLimitsMask[Index] > 0.0){
                    nonempty_indices[loopCount] = Index;
                    
                    nonempty_cellx[loopCount]   = (i+0.5)*CellSize;
                    nonempty_celly[loopCount]   = (j+0.5)*CellSize;
                    nonempty_cellz[loopCount]   = (k+0.5)*CellSize;
                    
                    loopCount += 1;
                }
            }
        }
    }
    

    sprintf(filepath, "/disk1/mjw/HOD_MockRun/Data/VIPERS_window2/randoms_W1_Nagoya_xyz_0.7_0.8_gridded.cat");

    output = fopen(filepath, "w");

    while(rand_count<rand_number){
        Index       = (int) gsl_rng_uniform_int(gsl_ran_r, nonempty_cells);        
    
        x           = nonempty_cellx[Index] + CellSize*(gsl_rng_uniform(gsl_ran_r) - 0.5);
        y           = nonempty_celly[Index] + CellSize*(gsl_rng_uniform(gsl_ran_r) - 0.5);
        z           = nonempty_cellz[Index] + CellSize*(gsl_rng_uniform(gsl_ran_r) - 0.5);
        
        rand_count += 1;
    
        fprintf(output, "%e \t %e \t %e \n", x, y, z);
    }

    fclose(output);
    
    return 0;
}


int deplete_homogeneous(){
    //  This will work terribly for density fields with \delta_Max >> 1.
    double nbar, maxDensity, alpha, acceptance;
    
    if(arrayMin(densityArray, n0*n1*n2) < -1.){
        printf("\n\nError, densities are not positive definite. exit.");
        
        return 0;
    }
    
    nbar       = 9.*pow(10., 4.)/TotalVolume;
    
    maxDensity = arrayMax(densityArray, n0*n1*n2);

    alpha      = ceil(1. + maxDensity);

    printf("\n\nalpha: %e", alpha);
    
    // poissonSample_homogeneous(alpha*TotalVolume*nbar);

    load_homogeneous_rands(alpha*TotalVolume*nbar, 1);


    sprintf(filepath, "%s/Data/stacpolly/poissonSampled_clustered_fog_500/poissonSampled_clustered_fog_500_%d.dat", root_dir, loopCount);

    output = fopen(filepath, "w");    	 

    for(j=0; j<rand_number; j++){ 
        xlabel                  = (int) floor((rand_x[j] - AxisLimsArray[0][0])/CellSize);    
        ylabel                  = (int) floor((rand_y[j] - AxisLimsArray[0][1])/CellSize);
        zlabel                  = (int) floor((rand_z[j] - AxisLimsArray[0][2])/CellSize);
    
        boxlabel                = (int)                 xlabel + n2*ylabel + n2*n1*zlabel;

        acceptance              = gsl_rng_uniform(gsl_ran_r);
        
        if(acceptance<(1. + densityArray[boxlabel])/alpha){
            fprintf(output, "%e \t %e \t %e \n", rand_x[j], rand_y[j], rand_z[j]);
        }
    }    
    
    return 0;
}


int load_clustered(int loadgals, double sampling){
    // sprintf(filepath, "%s/Data/stacpolly/poissonSampled_clustered_fog_500/poissonSampled_clustered_fog_500_%d.dat", root_dir, loopCount);
    sprintf(filepath, "%s/Data/stacpolly/NFW_halo_catalogue_%d.dat", root_dir, loopCount);
    
    inputfile   = fopen(filepath, "r");

    ch          = 0;
    Vipers_Num  = 0;

    do{
        ch = fgetc(inputfile);
        
        if(ch == '\n')
            Vipers_Num += 1;
    } while(ch != EOF);

    // last line of catalogue is not printing correctly.
    Vipers_Num -= 1;
    
    lowerSampling_randomisedCatalogue(sampling);

    printf("\n\n%d Vipers number", Vipers_Num);

    if(loadgals ==1){
        rewind(inputfile);

        xCoor          =  (double *)  realloc(xCoor,          Vipers_Num*sizeof(*xCoor));
        yCoor          =  (double *)  realloc(yCoor,          Vipers_Num*sizeof(*yCoor));
        zCoor          =  (double *)  realloc(zCoor,          Vipers_Num*sizeof(*zCoor));
        
        hod_disp       =  (double *)  realloc(hod_disp,       Vipers_Num*sizeof(*hod_disp));

        for(j=0; j<Vipers_Num; j++)   fscanf(inputfile, "%le \t %le \t %le \t %le", &xCoor[j], &yCoor[j], &zCoor[j], &hod_disp[j]);
    }

    fclose(inputfile);
    
    for(j=0; j<10; j++)  printf("\n%e \t %e \t %e \t %e", xCoor[j], yCoor[j], zCoor[j], hod_disp[j]);
    
    // randomiseCatalogue(Vipers_Num, xCoor, yCoor, zCoor, hod_disp);

    return 0;
}


int poissonSample_fog(){
    // homogeneous stays homogeneous under fingers of god. 
    sprintf(filepath, "%s/Data/lnnorm_anisotropic/poissonSampled_randoms_250_fog.dat", root_dir);

    output = fopen(filepath, "w");    	 

    int GalNumber;
    
    double x, y, z, nbar, lambda, dx, dy, dz, r_fog;
    
    nbar = 2.*pow(10., 4.)/TotalVolume;
    
    for(k=0; k<n0; k++){      
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){    
	            x         = AxisLimsArray[0][2] + i*CellSize;
	            y         = AxisLimsArray[0][1] + j*CellSize; 
	            z         = AxisLimsArray[0][0] + k*CellSize;
	            
    	        Index     = k*n1*n2 + j*n2 + i;
    
                lambda    = nbar*CellVolume*(1. + densityArray[Index]);
    
                GalNumber = gsl_ran_poisson(gsl_ran_r, lambda);
    
                for(ii=0; ii<GalNumber; ii++){  
                    dx    = CellSize*gsl_rng_uniform(gsl_ran_r);
	                dy    = CellSize*gsl_rng_uniform(gsl_ran_r);
	                dz    = CellSize*gsl_rng_uniform(gsl_ran_r);          
                
                    fprintf(output, "%e \t %e \t %e \n", x + dx, y + dy, z + dz);
                }
           }
       }
   }
   
   fclose(output);

   return 0;
}


int randomiseCatalogue(int objectNumber, double xvals[], double yvals[], double zvals[], double disp[]){
    shuffle_rows = realloc(shuffle_rows, objectNumber*sizeof(int));

    rows_randomOrder(shuffle_rows, objectNumber);
    
    
    shuffle(xvals,   shuffle_rows, objectNumber);
    
    shuffle(yvals,   shuffle_rows, objectNumber);
    
    shuffle(zvals,   shuffle_rows, objectNumber);
    
    // shuffle(disp,    shuffle_rows, objectNumber);
    
    // printf("\n\n");

    // for(j=0; j<Vipers_Num; j++)  printf("\n%e \t %e \t %e", xCoor[j], yCoor[j], zCoor[j]);
    
    // output = fopen(filepath, "w");
    
    // for(j=0; j<objectNumber; j++)  fprintf(output, "%e \t %e \t %e \t %e \n", xvals[j], yvals[j], zvals[j], disp[j]);
    
    // fclose(output);
    
    return 0;
}


int rows_randomOrder(int* shuffle_rows, int objectNumber){
    for(j=0; j<objectNumber; j++)  shuffle_rows[j] = j;
    
    gsl_ran_shuffle(gsl_ran_r, shuffle_rows, objectNumber, sizeof(int));

    return 0;
}


int shuffle(double array[], int shuffle_rows[], int N){
    double temp[N];
    
    for(k=0; k<N; k++)  temp[k]                = array[shuffle_rows[k]];
 
    for(k=0; k<N; k++)  array[k]               = temp[k];
    
    return 0;
}


int shuffletest(){
   double  array[100];
   int     shuffle_rows[100];
   
   for(j=0; j<100; j++) array[j] = 1.*j;
  
   shuffle(array, shuffle_rows, 100);
    
   return 0;
}


int toyTrees(){
    rand_number =  10000;
    Vipers_Num  =   5000; 
    
    assign_randmemory();

    xCoor          =  (double *)  realloc(xCoor,          Vipers_Num*sizeof(*xCoor));
    yCoor          =  (double *)  realloc(yCoor,          Vipers_Num*sizeof(*yCoor));
    zCoor          =  (double *)  realloc(zCoor,          Vipers_Num*sizeof(*zCoor));

    for(j=0; j<Vipers_Num; j++)  xCoor[j]   = (AxisLimsArray[1][0] - AxisLimsArray[0][0])*gsl_rng_uniform(gsl_ran_r); 
    for(j=0; j<Vipers_Num; j++)  yCoor[j]   = (AxisLimsArray[1][1] - AxisLimsArray[0][1])*gsl_rng_uniform(gsl_ran_r); 
    for(j=0; j<Vipers_Num; j++)  zCoor[j]   = (AxisLimsArray[1][2] - AxisLimsArray[0][2])*gsl_rng_uniform(gsl_ran_r); 

    for(j=0; j<rand_number; j++)  rand_x[j] = (AxisLimsArray[1][0] - AxisLimsArray[0][0])*gsl_rng_uniform(gsl_ran_r); 
    for(j=0; j<rand_number; j++)  rand_y[j] = (AxisLimsArray[1][1] - AxisLimsArray[0][1])*gsl_rng_uniform(gsl_ran_r); 
    for(j=0; j<rand_number; j++)  rand_z[j] = (AxisLimsArray[1][2] - AxisLimsArray[0][2])*gsl_rng_uniform(gsl_ran_r); 

    printf("\n\n%e \t %e \t %e", arrayMax(rand_x, rand_number), arrayMax(rand_y, rand_number), arrayMax(rand_z, rand_number));

    // printf("\n\nVipers gals: ");
    // for(j=0; j<Vipers_Num; j++)  printf("\n%e \t %e \t %e", xCoor[j], yCoor[j], zCoor[j]);
    
    // printf("\n\nRands: ");
    // for(j=0; j<rand_number; j++)  printf("\n%e \t %e \t %e", rand_x[j], rand_y[j], rand_z[j]);

    return 0;
}*/
