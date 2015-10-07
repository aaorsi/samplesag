/**
 * @file filters.c
 * @brief This file contains all the subroutines related to the handling 
 * of filters 
 * @author Alvaro Orsi
 */ 
#include <stdlib.h>
#include <stdio.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <string.h>
#include <math.h>

#include "allvars.h"

void get_filterinfo(char *info_file)
{

  FILE * ff;
  char buf[SED_MAXC],info[SED_MAXC];
  char *p;
  int id,j,k; 
  
  
  if (!(ff = fopen(info_file,"r"))) {
    printf("DUMPHDF5_SED.C: Filter info file %s not found. Check!\n",info_file);
    exit(0);
  } 		

  j = 0;
  k = 0;
  // The first three lines are header+blank space, so they are irrelevant here.
  // If by any chance this changes, then this should be modified accordingly.
  
  fgets(buf,SED_MAXC,ff);
  fgets(buf,SED_MAXC,ff);
  fgets(buf,SED_MAXC,ff);
  
  while(!feof(ff)) {
    fgets(buf,SED_MAXC,ff);
    p = strtok(buf,"#");
    p = strtok(NULL,"#");
    id = k;
    
    //		printf(" id = %d Description: %s \n",id,p);
    
    for (j = 0; j < NMags; j++) {
      if ( id == idMag[j] ) {	
	Idinfo[j] = (char *) malloc((SED_MAXC+1) * sizeof(char));
	strcpy(Idinfo[j],p);
      }	
    }
    k++;
  }
  fclose(ff);
  
  printf("Magnitudes to be computed using the SED\n");
  for (j = 0;j<NMags;j++) {
    printf("%s",Idinfo[j]);
    fflush(stdout);
  }

  return;
}

void read_filterfile(char * ffile)
{
  // ffile is the filter file, fdata is an output containing data 
  
  FILE * ff, *fd;
  char buf[SED_MAXC],buf2[SED_MAXC],fname[SED_MAXC];	
  char *p;
  int debug = 1;
  int ifilter,isf;	
  
  //	Check for files

  if ( !(ff = fopen(ffile,"r"))) {
    printf("Cannot open file %s. Check!\n",ffile);
    exit(0);
  }

  ifilter = -1;	
  
  while (!feof(ff)) {	
    
    // Here it reads each line of the file and checKs whether it corresponds to a header or filter data.
    fgets(buf,SED_MAXC,ff);
    
    // Make sure we do not try to read the file when it has no more lines.
    if (feof(ff))
      break;
    strcpy(buf2,buf);
    p = strtok(buf," \t");
    
    // If header, then stores index and writes header to file fdata.
    if (strcmp(p,"#") == 0) {
      ifilter++;
      if (ifilter > 0) {
	Nsteps_filter[ifilter-1] = isf;
//	dprinti(Nsteps_filter[ifilter-1]);
      }
      strcpy(fname,buf2);
      isf = 0;
    }
    else {
      // Else, it stores wavelength and transmission in Filterlist array.
      
      //			Filterlist[isf + NFILTERL*(2*ifilter)] = atof(p);
      Filterlist[2*isf + 2*NFILTERL*ifilter] = atof(p);
      p = strtok(NULL," ");
      Filterlist[2*isf+1 + 2*NFILTERL*ifilter] = atof(p);
      //			Filterlist[isf+NFILTERL*(2*ifilter+1)] = atof(p);
      isf++;
    }
    if (feof(ff))
      break;
    
  }

  if (ifilter > 0) 
  {
  	Nsteps_filter[ifilter] = isf;
	dprinti(Nsteps_filter[ifilter]);
  }
  printf("Filter file read and stored \n");
  if (MAKE_FD_DATA ==1) {	
    printf("Data filter created\n");
    fclose(fd);
  }
  
  fclose(ff);

  return;
}


float trapz1(float * x, float * y, int n)
{
  int i;
  double res = 0;

  if (n < 1)
    return(0);
  
  for (i = 1; i<n; i++)
    res += fabs(x[i] - x[i-1]) * (y[i] + y[i-1])/2.;
  if (n == 1)
    res = x[0] * y[0] / 2;
  
  return(res);
}


/**
 * @param id Filter ID
 * @param spec Spectrum from which the magnitude is computed
 * @param frame Sets the magnitude in rest-frame (0) or observer-frame (1)
 * @param z Redshift
 *
 * This subroutine is based on function AB_MAG and FILTR_MASTER from the 
 * Charlot & Bruzual 2007 code.
 */	
float compute_mag_filter(int id, int nlambda, float * wlambda, 
			 float * spec, int frame, float z, int iff, int AB)
{
  FILE *fout;	
  char filename[] = {"test_mag/file_test"};
  int ns,i,j;
  float Zero,AB0, dl,rz,Mag;
  float * int_sed;
  float iw, ir;
  double whigh,wlow,hhigh,hlow;
  double c_in_angstroms = 2.99e18;
  float * rfa, *rfb;
  //	ns = Nsteps_filter[id];
  double f_mean;
  float corr;
  
  dl = 1.e-5;		// 10 pc in units of Mpc.
  //  This factor is SQRT(4*pi*(3.0856E24)^2/Lsun)
  
  AB0 = 5. * (8.0 + log10(1.7684 * dl));
  Zero = (AB == 1) ? AB0 : AB0 -MVega[id];	
  
  //	printf("MVega[%d] = %f\n",id,MVega[id]);
  
  /*	if (id == 16)
    {	
    printf("MVega[%d] = %f\n",id,MVega[id]);
    }
  */	
  int_sed = (float *) malloc(nlambda * sizeof(float));
  rfa		= (float *) malloc(nlambda * sizeof(float));
  rfb     = (float *) malloc(nlambda * sizeof(float));

  for (i = 0; i<nlambda; i++) {
    int_sed[i]  = 0;
    rfa[i]      = 0;
    rfb[i]		= 0;
  }
    
  if (iff == 1) {
    sprintf(filename,"%s.%d",filename,id);
    
    fout = fopen(filename,"w");
    
    fprintf(fout,"# Filter %d\n",id);
    for (i=0;i<nlambda;i++)
      fprintf(fout,"%g %g\n",Filterlist[2*i+2*NFILTERL*id],	\
	      Filterlist[2*i+1+2*NFILTERL*id]);
    fprintf(fout,"# Scaled Sed\n");
    for (i=0;i<NBinSed;i++)
      fprintf(fout,"%g %g\n",wlambda[i],spec[i]);
    fprintf(fout,"# Interpolated SED\n");
  }
  
  j = 0;
  for (i = 0 ;i < nlambda; i++) {
    iw = Filterlist[2*i+2*NFILTERL*id];
    ir = Filterlist[(2*i+1)+2*NFILTERL*id];
    
    if (frame == 1)
      iw = iw / (1.0 + z);

  	if (i == 0)
	{ 
		if (iw < wlambda[0])
		{
			Mag = 100.0;
			return(Mag);
		}
	}
 
    while(wlambda[j] <= iw  && j < SED_MAXW) {
      wlow = wlambda[j];
      hlow = spec[j];
      j++;
    } 
    whigh = wlambda[j];
    hhigh = spec[j];
    
    int_sed[i] = (hhigh - hlow) / (whigh - wlow) *		\
      (iw - whigh) + hhigh;	
    
    if (iff == 1)
      fprintf(fout,"%g %g\n",iw,int_sed[i]);
    // Here two aux. varIABles are defined
    rfa[i] = ir * int_sed[i] * iw;
    rfb[i] = ir / iw;
    
  }	
  if (iff == 1)
    fclose(fout);
  
  /* The mean flux is computed as follows:
     f_mean = Int( dLam  Rf_lam  F_lam  lam /c ) / Int ( dLam Rf_Lam /Lam) 
  */
  
  f_mean = trapz1(&Filterlist[NFILTERL*(2*id)],rfa,nlambda) /		\
    trapz1(&Filterlist[NFILTERL*(2*id)],rfb, nlambda)			\
    /c_in_angstroms;
  
  
  Mag = (f_mean == 0 ) ?  100. : Zero - 2.5 * log10(f_mean) - 48.60;	
  
  if (iff == 1) {
    printf("int rfa %g int rfb %g fmean %g \n",trapz1(&Filterlist[NFILTERL*(2*id)],rfa,nlambda), \
	   trapz1(&Filterlist[NFILTERL*(2*id)],rfb,nlambda),f_mean);  
    
    printf("AB0 %g Mag %g \n",AB0,Mag);
  }
    
  if (isnan(Mag) || isinf(Mag)) {
    printf("WARNING: Mag is Inf or NaN\n");
    printf("Mag = %f\n",Mag);
    exit(0);
  }		

  free(rfa);
  free(rfb);
  free(int_sed);
  
  return(Mag);
}


/**
 * @brief This subroutine computes the production rate of Lyman-continuum 
 * photons.
 * 
 * This corresponds to @f$\dot{N}_{Lyc} = \int_{\nu_0}^{\infty}) L_{\nu} / 
 * h\nu d\nu  = 1 / (hc) * \int_0^{912\AA{}} L_{\lambda} * 
 * \lambda d\lambda @f$. Notice the units of nlyc : [erg s@f$^{-1}@f$/1e30]
 */	
double compute_lyc(float * wlambda, float * spec)
{
  int i,j;
  double c_in_angstroms = 2.99792e18;
  double h_pl = 6.6262e-27;
  double lsun = 3.826e3;		//erg s-1 * 1e30
  double fact = lsun/(h_pl * c_in_angstroms);
  //	double fact = 1./h_pl * 1e-40;
  double wly  = 912.0;
  float *w, *f;
  double nlyc;	
  

  w = (float *) malloc(NBinSed * sizeof(float));
  f = (float *) malloc(NBinSed * sizeof(float));
  
  i = 0;
  nlyc = 0.0;
  
  while(wlambda[i] < wly) {
    w[i] = wlambda[i];
    f[i] = wlambda[i] * spec[i];
    i++;
  }
  
  w[i] = wlambda[i];
  f[i] = wlambda[i] * spec[i];
  
  nlyc = fact * trapz1(w,f,i);
  /*	
  // nlyc calculated as a simple trapezoidal rule
  nlyc += fabs(wlambda[i] - wlambda[i-1]) * (spec[i]*wlambda[i] + spec[i-1]*wlambda[i-1])/2.;
  
  i++;
  }
  
  nlyc *= fact;
  
  */
  
#ifdef DEBUG
  if (isnan(nlyc) || isinf(nlyc)) {
    printf("filters.c:compute_lyc WARNING: nlyc calculation gives Inf or NaN\n");
    exit(0);
  }
#endif
  
  free(w);
  free(f);
  
  return(log10(nlyc) + 30.0);
}

/**
 * @fn continuum_line
 * @brief This computes the continuum around a given emission line simply 
 * as the average continuum at either side of the line.
 *
 * 'lam_off' is the offset with respect to the line centre, and 'width' is
 * the width of the continuum used to get an average at each side 
 */
float continuum_line(float * wlambda, float * spec, float lambda0)
{
  float lam_off = 0;
  float width   = 100;
  long i, j;
  float wr0,wr1,wr2,wr3;
  float av1, av2,cont_line;
  float lsun = 3.826e-7; //Luminosity of the sun /1e40 
  /*
    wr0 = lambda0 - (lam_off + width/2.);
    wr1 = lambda0 - (lam_off - width/2.);
    wr2 = lambda0 + (lam_off -width/2.);
    wr3 = lambda0 + (lam_off + width/2.);
  */
  
  cont_line = 0;
  av1 = 0;
  av2 = 0;
  
  i = 0;
  j = 0;
  while (wlambda[i] < lambda0)
    i++;
  
  av1 = spec[i-1];
  av2 = spec[i];
  
  /*	while(wlambda[i] >= wr0 && wlambda[i] < wr1)
    {
    av1 += spec[i];
    i++;
    j++;
    }
    av1 = av1/( (float) j) ;	
    
    
    j = 0;
    
    while(wlambda[i] <	wr2)
    i++;
    
    while(wlambda[i] >= wr2 && wlambda[i] < wr3)
    {
    av2 += spec[i];
    i++;
    j++;
    }
    av2 = av2/( (float) j) ;	
  */
  
  cont_line = (0.5 * (av1 + av2))*lsun;
  return(cont_line);	
}


/**
 * @brief This function returns A(\lambda)/A(V), where \lambda = 1/x, 
 * based on the extinction curve of Cardelli et al. 1989. 
 */
float cardelli_ext(float x)
{
  float Al_AV;
  float a,b,y,Fa,Fb;
  float x2,x3;
  float x8s,x8c;
  float y2,y3,y4,y5,y6,y7;
  
  //  The following expressions correspond to Eqs. 2, 3, 4 and 5 from C89
  
  if(x >= 0.3 && x < 1.1) {	// Infrared
    a = 0.574 * pow(x,1.61);
    b = -0.527*pow(x,1.61);	
  }
  else if(x >=1.1 && x < 3.3) { // Optical
    y = x-1.82;
    y2 = y*y;
    y3 = y2*y;
    y4 = y3*y;
    y5 = y4*y;
    y6 = y5*y;
    y7 = y6*y;
    
    a = 1 + 0.17699*y - 0.50447*y2 - 0.02427*y3 +	       \
      0.72085*y4 + 0.01979*y5 - 0.77530*y6 +		       \
      0.32999*y7;
    b = 1.41338*y + 2.28305*y2 + 1.07233*y3 - 5.38434*y4 -	\
      0.62251*y5 + 5.30260*y6 - 2.09002*y7;
  }
  else if (x >= 3.3 && x < 8) { // Ultraviolet
    if (x < 5.9) {
      Fa = 0;	
      Fb = 0;
    }
    else {
      x2 = pow((x-5.9),2);	
      x3 = x2 * (x-5.9);
      Fa = -0.04473 * x2 - 0.009779*x3;
      Fb = 0.2130*x2 + 0.1207*x3;
    }
    
    a = 1.752 - 0.316*x - 0.104/(pow((x-4.67),2) + 0.341) + Fa;
    b = -3.090 + 1.825*x + 1.206/(pow((x-4.62),2)+0.263) + Fb;
  }
  else if (x >= 8 && x < 10) { //Far-UV
    x8s = (x-8)*(x-8);
    x8c = x8s*(x-8);
    a = -1.073 - 0.628*(x-8) + 0.137*x8s - 0.070*x8c;
    b = 13.670 + 4.257*(x-8) - 0.420*x8s + 0.374*x8c;
  }
  else {
    a = 0;
    b = 0;
  }
  
  Al_AV = a + b/SED_RV;
  
  /*	if (Al_AV == 0 && x < 10)
    {	
    printf("WARNING: A(lambda)/AV is zero.\n");
    }
  */
  
  return Al_AV;
}


/**
 * @brief This should be flexible enough to allow different implementations 
 * of dust extinction. The first one is equivalent to the one already 
 * implemented in SAG, and follows De Lucia et al. 2004, and Kauffmann et 
 * al. 1999	
 */
float dust_extinction(float lambda, float theta,int p,int magbin, 
		      float *spec ,float magB)
{
  
  int idB;
  float x;	
  float lumB;
  float tauB,tauV,tauL;
  float aL,ct;	
  
  //	idB = 16;	// This corresponds to the B-Band, but check in 
  // filters.info file just in case! 
  
  //	printf("B-band taken to be idMag = %d\n",idB);
  
  //		MagB = compute_mag_filter(idB,Nsteps_filter[idB],	\
    //				W_scaled,spec,0,ZZ[Snapshot],0,0);
    
    x = 1.0e4/lambda;	//x [\mu m-1]
    
    if (magB == 100.0 || x > 10.0) {
      aL = 0;
      return aL;
    }
    
    lumB = pow(10,-magB/2.5);
    tauB = TAUBSTAR * pow(lumB/LBSTAR,BETA);
    tauV = tauB/Ab_Av;
    
    tauL = cardelli_ext(1.0e4/lambda) * tauV;
    ct   = cos(theta);
    //	AL 	 = -2.5 * log(ct*(1 - exp(-tauL/ct))/tauL);
    aL 	 = (tauL == 0 ) ? 0 : ct*(1 - exp(-tauL/ct))/tauL;
    
    if (isnan(aL)) {
      printf("AL in dust_extinction() is NaN. Check!");
      exit(0);
    }
	
    return aL;
}

void get_Vega_AB(char *VegaFile)
{
  float *wVega, *sVega;
  FILE * ff;
  long i,j;
  float w,s;
  double tot,stot,scl;
  float vmag,bc,vegabol;
  //   Vega apparent V magnitude
  vmag	= 0.03;
  bc		= -0.25;
  vegabol	= vmag+bc;
  
  wVega = (float *) malloc(10*SED_MAXW * sizeof(float));
  sVega = (float *) malloc(10*SED_MAXW * sizeof(float));
  
  if (!(ff = fopen(VegaFile,"r"))) {
    printf("FILTERS.C: Vega spectrum file %s not found. Exiting...\n",VegaFile);
    exit(0);
  }
  
  i = 0;
  
  while (!feof(ff)) {	
    fscanf(ff,"%f %f\n",&w,&s);
    wVega[i] = w;
    sVega[i] = s;
    //		fscanf(ff,"%f %f\n",&wVega[i],&sVega[i]);
    i++;
  }
  
  fclose(ff);	
  tot 	= trapz1(wVega,sVega,i-1);
  stot 	= pow(10,-0.4*(vegabol-4.75));
  scl		= stot/tot;
  
  for (j = 0; j<i-1; j++)
    sVega[j] *= scl;
  
  for (i = 0;i<NMags; i++) {	
    MVega[idMag[i]] = compute_mag_filter(idMag[i],Nsteps_filter[idMag[i]], \
					  wVega,sVega,0,0,0,1);
    printf("%ld MVega[%d] = %f\n", i, idMag[i], MVega[idMag[i]]);
  }
  //	exit(0);
  free(wVega);
  free(sVega);	

  return;
}


/**
 * The d4000 parameter is defined as
 * @f$D_{4000} = \int_{4000}^{4100} f_\lambda d_\lambda / \int_{3850}^{3950} 
 * f_\lambda d_\lambda @f$
 */
float get_d4000(float *wlambda, float *spec)
{
  int i,j;
  float d4000,d1,d2;
  float l1 = 3850;	
  float l2 = 3950;
  float l3 = 4000;
  float l4 = 4100;
  
  
  i = 0;
  d4000 = 0;
  d1 = 0;
  d2 = 0;
  
  while (wlambda[i] < l1)
    i++;
  
  while (wlambda[i] >= l1 && wlambda[i] < l2) {
    d1 += fabs(wlambda[i] - wlambda[i-1])*(spec[i]*wlambda[i] + \
					   spec[i-1]*wlambda[i-1])/2.0;
    i++;
  }
  
  while (wlambda[i] < l3)
    i++;
  
  while (wlambda[i] >= l3 && wlambda[i] < l4) {
    d2 += fabs(wlambda[i] - wlambda[i-1])*(spec[i]*wlambda[i] + \
					   spec[i-1]*wlambda[i-1])/2.0;
    i++;
  }
  
  d4000 = d2/d1;
  
  return(d4000);
}
