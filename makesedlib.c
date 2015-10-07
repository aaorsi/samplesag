/**
 * @file makesedlib.c
 * SEDs from a given stellar pop. model are read and stored in a convenient form.
 * For compatibility reasons, these SED are read from ASCII files.
 * @author Alvaro Orsi
 * @date 2011-07-06
 */
#include "proto.h"
#include "allvars.h"
#include "readparameterfile.h"


void read_sed(char *sed_file, float * w, float *h)
{
  
  FILE *fssp;
  int i,j,Inw2,nlc;
  char buf[SED_MAXCHAR];
  char * p, *q;
  float _age,_z,_w,_f; 
  double lsun = 3.826e33;
 
  if (!(fssp = fopen(sed_file,"r"))) {
    printf("Error (read_sed): cannot open file '%s'. Check!\n",sed_file);
    fflush(stdout);
    exit(EXIT_FAILURE);
  }
 
  if (strcmp(SSPModel,"CB07")==0)
  {
 
  fgets(buf,SED_MAXCHAR,fssp);
  p = strtok(buf," ");
  Ks = atoi(p);
  
//	printf("Ks = %ld\n",Ks);
//  AgeArray = (float *) malloc(Ks * sizeof(float));
  
  for (i = 0; i < Ks; i++) {
    p = strtok(NULL," ");
    AgeArray[i] = atof(p);
  }	
  
  fgets(buf,SED_MAXCHAR,fssp);
  fgets(buf,SED_MAXCHAR,fssp);
  fgets(buf,SED_MAXCHAR,fssp);
  fgets(buf,SED_MAXCHAR,fssp);
  fgets(buf,SED_MAXCHAR,fssp);
  
  fgets(buf,SED_MAXCHAR,fssp);
  
  p = strtok(buf," ");
  Inw = atoi(p);
  
  //	printf("Inw = %d\n",Inw);
  
  for (i = 0; i < Inw; i++) {		
    p = strtok(NULL," ");
    w[i] = atof(p);
  }
  
  for (i = 0; i<Ks; i++) {					
    fgets(buf,SED_MAXCHAR,fssp);
    p = strtok(buf," ");
    Inw2 = atoi(p);
    
    if (Inw != Inw2) {
      printf("Error (read_sed): something's wrong. Check\n");
      printf("Inw %d Inw2 %d\n",Inw,Inw2); fflush(stdout);
      exit(EXIT_FAILURE);
    }
    
    for (j = 0; j < Inw; j++) {
      p = strtok(NULL," ");	
      h[j+i*SED_MAXW] = atof(p);
    }
  }
  }
  
  if(strcmp(SSPModel,"Maraston05")==0)
  {

  	Ks = 67;
//    AgeArray = (float *) malloc(Ks * sizeof(float));
    char AgeFile[500];
  	FILE *AF; 
	sprintf(AgeFile,"%sAgegridSSP_Mar05",Sspdir);

	if (!(AF = fopen(AgeFile,"r"))) {
    printf("Error (read_sed): cannot open file '%s'. Check!\n",AgeFile);
    fflush(stdout);
    exit(EXIT_FAILURE);
    }
     
    for(i = 0;i<Ks;i++)
    {
		fscanf(AF,"%f",&_age);
		AgeArray[i] = _age * 1e9;
	}

 	fclose(AF);

	fscanf(fssp,"%f %f %f %f",&_age,&_z,&_w,&_f);	
	for (i = 0;i < Ks ; i++)
	{
		j = 0;
		while(_age*1e9 < AgeArray[i]+1 && _age*1e9 > AgeArray[i]-1)
		{
			w[j] = _w;
			h[j + i*SED_MAXW] = _f/lsun;
//			printf("%d h[j] = %f\n",j,h[j + i*SED_MAXW]);
			j++;
			fscanf(fssp,"%f %f %f %f",&_age,&_z,&_w,&_f);	
		}	
		if (i == 0)
			Inw = j-1;
	}

   }
		
  fclose(fssp);

  return;
}


void scale_sed(float *W0_, float *H0_, float *w, float *h,double Minw, \
	       double Maxw, int NBinSed)
{
  
  int i,j;
  double binsize,wlow,whigh,hlow,hhigh;
  int ScaleMethod;
  double binned_spec;
  //    double wly = 912.0;
  int imin;
  
  ScaleMethod = 2;
  
  /*  Test: Scale from wly upwards only 
      
  j = 0;
  while (W0_[j] < wly)
  {
  w[j] = W0_[j];
  h[j] = H0_[j];
  j++;
  }
  
  imin = j;
  printf("imin %d \n ",imin);
  //	exit(0);
  */
  
  binsize = (Maxw - Minw)/(NBinSed-1.0);
  //	binsize = (SED_MAXW - wly)/((NBinSed-imin)-1.0);
  
  j = 0;
  if (ScaleMethod == 1) {
    
    //		j = 0;
    for (i=0; i<NBinSed; i++) {
      w[i] = i * binsize + Minw;
      
      while (w[i] >= W0_[j] && j < Inw) {
	wlow = W0_[j];
	hlow = H0_[j];
	j++;
      }
      whigh = W0_[j];
      hhigh = H0_[j];				
      h[i] = (hhigh - hlow) / (whigh - wlow) * (w[i] - whigh) + hhigh;
      
      if (j >= Inw) {
	printf("Warning: scaled array exceeds original wavelength range\n");
	w[i] = W0_[Inw-1];
	h[i] = H0_[Inw-1];
      }
    }
  }
  else {
    //		j = 0;
    for (i = 0; i< NBinSed; i++) {
      binned_spec = 0;
      w[i] = i * binsize + Minw;
      
      while (W0_[j] < Minw - binsize/2.) {
	j++;
	if (j == Inw) {			
	  printf("scale_sed(): ERROR. Minw = %f from rebinned spectrum "
		 "outside SED range.\n", Minw); fflush(stdout);
	  exit(EXIT_FAILURE);
	}
      }
      
      while (W0_[j] <= w[i]+binsize/2. && W0_[j] > w[i] - binsize/2. && \
	     j < Inw) {		
	binned_spec += fabs(W0_[j] - W0_[j-1]) * (H0_[j] + H0_[j-1])/2.;
	j++;
      } 			
      h[i] = binned_spec/binsize;
    }
  }
  
  return;
}


/* To rescale the array of wavelengths and fluxes of each sed, this function
   assumes the wavelengths are in increasing order */
void init_sed()
{
  int i,j;
  char sed[256];
  
//  NBinSed 	= 200;
  NBinSed		= 200;
  Minw 		= 91;
  Maxw		= 25000;
  
  W0_ = (float *) malloc(SED_MAXW * sizeof(float));
  H0_ = (float *) malloc(NMetals * SED_MAXW * SED_MAXAGE * sizeof(float));
  
  W_scaled = (float *) malloc(NBinSed * sizeof(float));
  H_scaled = (float *) malloc(NMetals * NBinSed * SED_MAXAGE * sizeof(float));
 
  AgeArray = (float *)malloc(SED_MAXAGE * sizeof(float));
 
  for (i  = 0; i< NMetals; i++) {
    sprintf(sed,"%s%s",Sspdir,SSPFILE[i]);
    printf("Sed file: %s\n",sed);
    
    read_sed(sed,W0_,&H0_[SED_MAXAGE*SED_MAXW*i]);
    for (j= 0; j<Ks; j++)
      scale_sed(W0_,H0_+SED_MAXW*j + SED_MAXAGE*SED_MAXW*i,W_scaled,	\
		H_scaled+NBinSed*j + NBinSed*SED_MAXAGE*i,		\
		Minw,Maxw,NBinSed);
  }
  
  printf("SSP files read and scaled.\n");
  
  return;
}	

//int iz0,iz1,ia0,ia1;

//void compute_sed(float Z_SF, float ageSSP, double * flux_int)
void compute_sed(float Z_SF, float ageSSP, float * flux_int, int is)
{
  /* This function interpolates an SED using stored data as a function of 
     age and metallicity.
     Here flux_int corresponds to the interpolated SED  */
  
  int i,j;
  float z0,z1,age0,age1;
  double den,f00,f10,f01,f11,t1,t2,t3,t4;
  int iz0,iz1,ia0,ia1;

  ageSSP *= 1e9;
  
  if (Z_SF < ZArray[0]) {
    iz0 = 0;
    iz1 = 1;
    
    z0 = ZArray[iz0];	
    z1 = ZArray[iz1];
  }
  else {
    i = 0;	
    while(Z_SF >= ZArray[i] && i < NMetals-1 ) {
      z0 	= ZArray[i];
      iz0 = i;
      i++;
    }
    iz1 = i;
    z1 = ZArray[i];
  }
  
  if (ageSSP < AgeArray[0]) {
    ia0 = 0;
    ia1 = 1;
    age0 = AgeArray[ia0];	
    age1 = AgeArray[ia1];
  }
  else {	
    i = 0;	
    while(ageSSP >= AgeArray[i] && i < Ks-1 ) {
      ia0 = i;
      age0 = AgeArray[i];
      i++;
    }
    ia1 = i;
    age1 = AgeArray[i];
  }
  for (i = 0; i< NBinSed; i++) {
    
#ifdef IMF_VARIABLE
    f00 = H_scaled[i + NBinSed*ia0 + NBinSed*SED_MAXAGE*iz0 + NBinSed*SED_MAXAGE*SED_NZSSP*is];
    f10 = H_scaled[i + NBinSed*ia1 + NBinSed*SED_MAXAGE*iz0 + NBinSed*SED_MAXAGE*SED_NZSSP*is];
    f01 = H_scaled[i + NBinSed*ia0 + NBinSed*SED_MAXAGE*iz1 + NBinSed*SED_MAXAGE*SED_NZSSP*is];
    f11 = H_scaled[i + NBinSed*ia1 + NBinSed*SED_MAXAGE*iz1 + NBinSed*SED_MAXAGE*SED_NZSSP*is];
#else
    f00 = H_scaled[i + NBinSed*ia0 + NBinSed*SED_MAXAGE*iz0];
    f10 = H_scaled[i + NBinSed*ia1 + NBinSed*SED_MAXAGE*iz0];
    f01 = H_scaled[i + NBinSed*ia0 + NBinSed*SED_MAXAGE*iz1];
    f11 = H_scaled[i + NBinSed*ia1 + NBinSed*SED_MAXAGE*iz1];
#endif
    
    den = (z1 - z0) * (age1 - age0);
    
    t1 	= f00*(z1 - Z_SF)*(age1 - ageSSP)/den;
    t2  = f10*(Z_SF - z0)*(age1 - ageSSP)/den;
    t3  = f01*(z1 - Z_SF)*(ageSSP-age0)/den;
    t4  = f11*(Z_SF - z0)*(ageSSP-age0)/den;
    
    flux_int[i] = t1 + t2 + t3 + t4;
    
    if (flux_int[i] < 0) {
      //			printf("WARNING: SED interpolation gives flux < 0, flux_int[%d] = %f\n",i,flux_int[i]);
      //			printf("There must be something wrong.\n");
      //			printf("Setting flux_int[%d] = 0\n",i);	
      //			exit(0);
      flux_int[i] = 0.0;
    }
    
#ifdef DEBUG
    if (isnan(flux_int[i])) {
      printf("compute_sed() FATAL: flux_int[%d] = NaN. Check interpolation\n",i);
      fflush(stdout);
      exit(EXIT_FAILURE);
    }
#endif
  }		

#ifdef DEBUG	
  printf("Interpolation of SED done.\n");
  /*	for (i = 0;i<20;i++)
    {
    printf("flux_int[%d] = %g\n",i,flux_int[i]);
    if (isnan(flux_int[i]))
    printf("WARNING: NaN found!\n");
    }
    printf("\n");
  */
#endif	
  
  return;
}


/* 
   The following worKs to test the performance of the SED routines 
   
   int main()
{

	float z = 0.03;
	float age = 0.866;

	FILE *f0, *f1;
	char outfile0[] = {"/fdg/aaorsi/test/sed_int"};
	char outfile1[] = {"/fdg/aaorsi/test/seds0"};
	double *int_spec;
	int i;
	init_sed();
	
	int_spec = (double *) malloc(NBinSed * sizeof(double)); 
	compute_sed(z,age, int_spec);

	printf("SED computed. Now writing file %s\n",outfile0);

	if ( !(f0 = fopen(outfile0,"w")))
	{
		printf("Cannot open file %s. Check!\n",outfile0);
		exit(0);
	}

	fprintf(f0,"# lambda \t flux 00    flux 01   flux 10    flux 11   flux_int  \n");
	
	for (i = 0; i< NBinSed; i++)
	{	
		fprintf(f0,"%e  %e  %e  %e  %e  %e \n",w[i],h[i + ia0*NBinSed + iz0*NBinSed*MAXAGE],\
								h[i+ia0*NBinSed+iz1*NBinSed*MAXAGE],h[i + ia1*NBinSed + iz0*NBinSed*MAXAGE],\
								h[i+ia1*NBinSed+iz1*NBinSed*MAXAGE],int_spec[i]);

	}
	fclose(f0);

	printf("file %s written\n",outfile0);

	if ( !(f1 = fopen(outfile1,"w")))
	{
		printf("Cannot open file %s. Check!\n",outfile1);
		exit(0);
	}

	fprintf(f1,"# lambda \t flux 00    flux 01   flux 10    flux 11  \n");
	
	for (i = 0; i< Inw; i++)
	{
			fprintf(f1,"%g  %g  %g  %g  %g \n",W0[i],H0[i + ia0*SED_MAXW + iz0*SED_MAXW*MAXAGE],\
								H0[i+ia0*SED_MAXW+iz1*SED_MAXW*MAXAGE],H0[i+ia1*SED_MAXW+iz0*SED_MAXW*MAXAGE],\
								H0[i+ia1*SED_MAXW+iz1*SED_MAXW*MAXAGE]);

	}
	fclose(f1);

	printf("file %s written\n",outfile1);

	free(int_spec);
	return(0);
}	

*/
