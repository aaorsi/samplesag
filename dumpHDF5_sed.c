/**
 * @file dumpHDF5_sed.c
 */
#ifdef HDF5OUTPUT
#include "proto.h"
#include "allvars.h"
#include "readparameterfile.h"
#include <hdf5.h>
#include <hdf5_hl.h>


herr_t dumpHDF5_sed_popA(hid_t file_id, char *path, int magbin)
{
  int i, j, p, index,ised,mised,ii;
  double dbuf,ispec;
  char desc[256], label[256], subpath[256],comp[256];
  hid_t subgroup_id;
  float *buffer, *buffer2, *buffer3d, *bufferT,*buffer3, *buff4, *buff5, *buff6, *buff7, *buffW;
  int *ibuffer, *ibufferT;
  float *c_arr, *c_arr_ext;
  double *nlyc_arr;
  float currMag, dust_calc, ion_par, emlin, emlin_ext, qaux, magB, sfr;
  double line_lum, line_lum_ext;
  int idB = 16;
  hsize_t dims[2]={TotNumGalA,1};
  hsize_t dim3d[2]={TotNumGalA,3};
  hsize_t dimT[2]={TotNumGalA,OUTPUTS};
#ifdef OUTPUTSED
  hsize_t dimNBins[2]={TotNumGalA,NBinSed};
  hsize_t dimW[2]={NBinSed,1};
#endif
  herr_t status;
  
  buffer    = malloc(TotNumGalA*sizeof(float));
  nlyc_arr  = malloc(TotNumGalA*sizeof(double));
  buffer2   = malloc(TotNumGalA*sizeof(float));
  buffer3   = malloc(TotNumGalA*sizeof(float));
  c_arr	  = malloc(TotNumGalA*sizeof(float));
  c_arr_ext = malloc(TotNumGalA*sizeof(float));
  
  buff4	  = malloc(TotNumGalA*sizeof(float));
  buff5	  = malloc(TotNumGalA*sizeof(float));

#ifdef OUTPUTSED
  buff6	  = malloc(TotNumGalA*NBinSed*sizeof(float));
  buff7	  = malloc(TotNumGalA*NBinSed*sizeof(float));
  buffW   = malloc(NBinSed*sizeof(float));
#endif	

  // Magnitudes:
  // First, create subgroup called Magnitudes
  
  printf("Creating magnitudes and emission lines from SEDs\n");
  

  sprintf(subpath, "%s/Magnitudes", path);
  subgroup_id = H5Gcreate(file_id, subpath, H5P_DEFAULT, H5P_DEFAULT,
			  H5P_DEFAULT);
    
  mised = 2;

#ifdef OUTPUTSED

  sprintf(label,"%s/lambda",subpath);
  sprintf(desc, "wavelength array");
  for (ii = 0;ii < NBinSed; ii++)	
     buffW[ii] = W_scaled[ii];
  
  status = H5LTmake_dataset(file_id, label, 2, dimW, H5T_NATIVE_FLOAT, \
				buffW);
  dbuf = 0.0;
  status = set_float_variable_attr(file_id,label,desc,&dbuf,"Angstroms");

#endif

 
  for (ised = 0; ised<=mised; ised++) {
    if (ised == 0)
      sprintf(comp,"tot");
    if (ised == 1)
      sprintf(comp,"disk");
    if (ised == 2)
      sprintf(comp,"bulge");
    
    sprintf(label, "%s/D4000_%s",subpath,comp);
    sprintf(desc, "4000 A-break measurement.");
    
    for (p=1;p<= TotNumGalA; p++) {
      if (ised == 0) {
	for (ii = 0;ii < NBinSed; ii++)	{	
	  H_aux[ii] = GalA[p].Sed_disc[magbin][ii] \
	  + GalA[p].Sed_bulge[magbin][ii];
	}
      }
      if (ised == 1) {
	for (ii = 0; ii<NBinSed; ii++)
	  H_aux[ii] = GalA[p].Sed_disc[magbin][ii];
      }
      if (ised == 2) {
	for (ii = 0; ii<NBinSed; ii++)
	  H_aux[ii] = GalA[p].Sed_bulge[magbin][ii];
      }
      
      buffer[p-1] = get_d4000(W_scaled,H_aux);
    }
    
    status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_FLOAT, 
			      buffer);
    status = set_float_variable_attr(file_id,label,desc,&dbuf,"??");
    
    
    for (i=0; i<Nbands; i++) {
      if (Iframe[i] == 0) {
	sprintf(label, "%s/Mag_id%d_%s_r",subpath,Idband[i],comp);
	sprintf(desc, "%s, rest-frame",Idinfo[i]);
      }
      else {
	sprintf(label, "%s/Mag_id%d_%s_o",subpath,Idband[i],comp);
	sprintf(desc, "%s, observer-frame",Idinfo[i]);
      }
      
      for (p=1; p<=TotNumGalA; p++) {	
	if (ised == 0) {
	  for (ii = 0;ii < NBinSed; ii++) {	
	    H_aux[ii] = GalA[p].Sed_disc[magbin][ii];
	    H_aux[ii] += GalA[p].Sed_bulge[magbin][ii];
	  }
	}
	if (ised == 1) {
	  for (ii = 0; ii<NBinSed; ii++)
	    H_aux[ii] = GalA[p].Sed_disc[magbin][ii];
	}
	
	if (ised == 2) {
	  for (ii = 0; ii<NBinSed; ii++)
	    H_aux[ii] = GalA[p].Sed_bulge[magbin][ii];
	}
	
	currMag = compute_mag_filter(Idband[i],Nsteps_filter[Idband[i]],
				     W_scaled,H_aux,Iframe[i],ZZ[Snapshot],
				     0,IAB[i]);
	buffer[p-1] = currMag;
      }
      
      status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_FLOAT,
				buffer);
      dbuf = 0.0;
      status = set_float_variable_attr(file_id,label,desc,&dbuf,"??");
    }
    
    sprintf(label, "%s/D4000_%s_ext",subpath,comp);
    sprintf(desc, "4000 A-break + extinction  measurement.");
    
    for (p=1;p<= TotNumGalA; p++) {
      if (ised == 0) {
	for (ii = 0;ii < NBinSed; ii++) {	
	  H_aux[ii] = GalA[p].Sed_disc[magbin][ii];
	  H_aux[ii] += GalA[p].Sed_bulge[magbin][ii];
	}
      }
      if (ised == 1) {
	for (ii = 0; ii<NBinSed; ii++)
	  H_aux[ii] = GalA[p].Sed_disc[magbin][ii];
      }
      if (ised == 2) {
	for (ii = 0; ii<NBinSed; ii++)
	  H_aux[ii] = GalA[p].Sed_bulge[magbin][ii];
      }
      
      magB = compute_mag_filter(idB,Nsteps_filter[idB], W_scaled,H_aux,
				0,ZZ[Snapshot],0,0);
      
      for (j = 0; j < NBinSed; j++) {	
	dust_calc = dust_extinction(W_scaled[j],GalA[p].Inclination,p,magbin,H_aux,magB);
	Sed_ext[j] = H_aux[j] * dust_calc;
 #ifdef OUTPUTSED
if (ised ==0) {
	  index = indexrm(p-1,j,TotNumGalA,NBinSed) + NBinSed;
	  buff6[index] = H_aux[j];
	  buff7[index] = dust_calc;
}
#endif
  
   }
      
      buffer[p-1] = get_d4000(W_scaled,Sed_ext);
    }
    
    status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_FLOAT, 
			      buffer);
    status = set_float_variable_attr(file_id,label,desc,&dbuf,"??");

#ifdef OUTPUTSED
if (ised ==0) {
 	  sprintf(label, "%s/Sed",subpath);
	  sprintf(desc, "Spectral energy distribution");
	
      status = H5LTmake_dataset(file_id, label, 2, dimNBins, H5T_NATIVE_FLOAT, \
				buff6);
      dbuf = 0.0;
      status = set_float_variable_attr(file_id,label,desc,&dbuf,"??");

 	  sprintf(label, "%s/dust_ext",subpath);
	  sprintf(desc, "Dust Extinction");
	
      status = H5LTmake_dataset(file_id, label, 2, dimNBins, H5T_NATIVE_FLOAT, \
				buff7);
      dbuf = 0.0;
      status = set_float_variable_attr(file_id,label,desc,&dbuf,"??");
}
#endif

    
    for (i=0; i<Nbands;i++) {
      if (Iframe[i] == 0) {
	sprintf(label, "%s/Mag_ext_id%d_%s_r",subpath,Idband[i],comp);
	sprintf(desc, "%s + Dust extinction, rest-frame",Idinfo[i]);
      }
      else {
	sprintf(label, "%s/Mag_ext_id%d_%s_o",subpath,Idband[i],comp);
	sprintf(desc, "%s + Dust extinction, observer-frame",Idinfo[i]);
      }
            
      for (p=1; p<=TotNumGalA; p++) {
	if (ised == 0) {
	  for (ii = 0;ii < NBinSed; ii++) {	
	    H_aux[ii] = GalA[p].Sed_disc[magbin][ii];
	    H_aux[ii] += GalA[p].Sed_bulge[magbin][ii];
	  }
	}
	if (ised == 1) {
	  for (ii = 0; ii<NBinSed; ii++)
	    H_aux[ii] = GalA[p].Sed_disc[magbin][ii];
	}
	
	if (ised == 2) {
	  for (ii = 0; ii<NBinSed; ii++)
	    H_aux[ii] = GalA[p].Sed_bulge[magbin][ii];
	}
	magB = compute_mag_filter(idB,Nsteps_filter[idB],		\
				  W_scaled,H_aux,0,ZZ[Snapshot],0,0);
	
	for (j = 0; j < NBinSed; j++) {	
	  dust_calc = dust_extinction(W_scaled[j],GalA[p].Inclination,p,magbin,H_aux,magB);
	  Sed_ext[j] = H_aux[j] * dust_calc;
	}

	
	currMag = compute_mag_filter(Idband[i],Nsteps_filter[Idband[i]], \
				     W_scaled,Sed_ext,Iframe[i],ZZ[Snapshot],0,IAB[i]);
	
	buffer[p-1] = currMag;
      }
      status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_FLOAT, \
				buffer);
      dbuf = 0.0;
      status = set_float_variable_attr(file_id,label,desc,&dbuf,"??");
    }
  }
  
  status = H5Gclose(subgroup_id);

  // Emission Lines:
  // First, create subgroup called Emission Lines
  
  sprintf(subpath, "%s/EmissionLines", path);
  subgroup_id = H5Gcreate(file_id, subpath, H5P_DEFAULT, H5P_DEFAULT,
			  H5P_DEFAULT);
  
  for (ised = 0; ised<=mised; ised++) {
    if (ised == 0)
      sprintf(comp,"tot");
    if (ised == 1)
      sprintf(comp,"disk");
    if (ised == 2)
      sprintf(comp,"bulge");
        
    // 		First, Lyman continuum calculation
    for (p=1; p<=TotNumGalA; p++) {

      /*			nlyc_arr[p-1] = 0;
	
      if (GalA[p].ColdGas == 0)
      continue; 
      */
      if (ised == 0) {
	for (ii = 0;ii < NBinSed; ii++) {	
	  ispec = GalA[p].Sed_disc[magbin][ii];
	  H_aux[ii] = (ispec < 0.0) ? 0.0 : ispec;
	  //					H_aux[ii] = GalA[p].Sed[magbin][ii];
	  //					if (H_aux[ii] < 0.0) H_aux[ii] = 0.0;
	  H_aux[ii] += GalA[p].Sed_bulge[magbin][ii];
	}
	
	sfr = GalA[p].Sfr[Snapshot];
      }
      if (ised == 1) {
	for (ii = 0; ii<NBinSed; ii++) {	
	  ispec = GalA[p].Sed_disc[magbin][ii];
	  H_aux[ii] = (ispec < 0.0) ? 0.0 : ispec;
	  //					H_aux[ii] = GalA[p].Sed[magbin][ii] - GalA[p].Sed_bulge[magbin][ii];
	  //					if (H_aux[ii] < 0.0) H_aux[ii] = 0.0; 
	}
	
	sfr = GalA[p].SfrQuies[Snapshot];
      }
      if (ised == 2) {
	for (ii = 0; ii<NBinSed; ii++) {
	  ispec = GalA[p].Sed_bulge[magbin][ii];
	  H_aux[ii] = (ispec < 0.0) ? 0.0 : ispec;
	  // 	H_aux[ii] = (GalA[p].Sed_bulge[magbin][ii] < 0.0 ) ? 0.0 : GalA[p].Sed_bulge[magbin][ii];
	}
	sfr = GalA[p].SfrBulgeBurst[Snapshot];
      }
      
      nlyc_arr[p-1] = compute_lyc(W_scaled,H_aux);
      //			nlyc_arr[p-1] = Sfr * 10.2640 * 1.0e23 / 1.08;
    }
    
    sprintf(label, "%s/log_NLyc_%s",subpath,comp);
    sprintf(desc, "log10 of number rate of Lyman Continuum photons from %s component [s-1]",comp);
    status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_DOUBLE, 
			      nlyc_arr);
    status = set_float_variable_attr(file_id,label,desc,&dbuf,"??");
    
    //	    Emission Lines
    
    for (i = 0; i<NLines; i++) {
      for (p = 1;p<=TotNumGalA;p++) {
	
	//				nlyc_arr[p-1] = 0;
	
	//				if (GalA[p].ColdGas == 0)
	//				continue; 
	
	if (ised == 0) {
	  for (ii = 0;ii < NBinSed; ii++) {	
	    H_aux[ii] = GalA[p].Sed_disc[magbin][ii];
	    H_aux[ii] += GalA[p].Sed_bulge[magbin][ii];
	  }
	}
	if (ised == 1) {
	  for (ii = 0; ii<NBinSed; ii++)
	    H_aux[ii] = GalA[p].Sed_disc[magbin][ii];
	}
	if (ised == 2) {
	  for (ii = 0; ii<NBinSed; ii++)
	    H_aux[ii] = GalA[p].Sed_bulge[magbin][ii];
	}
	
	//				nlyc_arr[p-1] = compute_lyc(W_scaled,H_aux);
	
	magB = compute_mag_filter(idB,Nsteps_filter[idB],		\
				  W_scaled,H_aux,0,ZZ[Snapshot],0,0);
	
	Mean_q = 0;

	if (AgeEvolemlines == 0)
	{
		if (ised == 1){
			line_lum = (isnan(nlyc_arr[p-1]) || isinf(nlyc_arr[p-1])) ? -30 : \
			integ_line(GalA[p].Rscale,GalA[p].MetallicityColdGasSF,nlyc_arr[p-1],i,ised,\
			GalA[p].ColdGas,GalA[p].DiscMass);
		}
		if (ised ==2){
			line_lum = (isnan(nlyc_arr[p-1]) || isinf(nlyc_arr[p-1])) ? -30 : \
			integ_line(GalA[p].Rbulge,GalA[p].MetallicityColdGasSF,nlyc_arr[p-1],i,ised,\
			GalA[p].BulgeMass,GalA[p].BulgeMass);
		}

		if (ised ==1 || ised ==2)	
		line_lum = (line_lum < -20) ? 1e-30 : pow(10,line_lum);
	}

	dust_calc = dust_extinction(CentralWavelength[i],GalA[p].Inclination,p,magbin,H_aux,magB);
	if (ised ==1 || ised ==2)
	line_lum_ext = line_lum * dust_calc; 
	
	
	for (j = 0; j < NBinSed; j++) {	
	  dust_calc = dust_extinction(W_scaled[j],GalA[p].Inclination,p,magbin,H_aux,magB);
	  Sed_ext[j] = H_aux[j] * dust_calc;
	}
	
	c_arr[p-1] 		= continuum_line(W_scaled,H_aux,CentralWavelength[i]);
	c_arr_ext[p-1]  = continuum_line(W_scaled,Sed_ext,CentralWavelength[i]);
	if (ised ==1 || ised ==2)
	{		
		buffer[p-1] 	= line_lum;
		buffer2[p-1]	= line_lum_ext;
	
		if (i ==0)
		{	
		  buffer3[p-1]    = Mean_q;
		  if (LineMod == 2)
		  {
		    buff4[p-1] = Qlow;
			buff5[p-1] = Qhigh;
   		  }
		}
  	}

   }	

	   if (ised ==1 || ised ==2)
	  {		      
 	     sprintf(label, "%s/%s_%s",subpath,Emlines_name[i],comp);
   	     sprintf(desc, "Luminosity, units of 10^40 erg s-1 h-2.");
   	     status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_FLOAT, 
				buffer);
	     status = set_float_variable_attr(file_id,label,desc,&dbuf,"??");
      
      	 sprintf(label, "%s/%s_%s_ext",subpath,Emlines_name[i],comp);
         sprintf(desc, "Luminosity + dust extinction, units of 10^40 erg s-1 h-2.");
         status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_FLOAT, 
				buffer2);
         status = set_float_variable_attr(file_id,label,desc,&dbuf,"??");
	  }
      
      sprintf(label, "%s/Cont_%s_%s",subpath,Emlines_name[i],comp);
      sprintf(desc, "Continuum luminosity, units of 10^40 erg s-1 h-2 A-1.");
      status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_FLOAT, 
				c_arr);
      status = set_float_variable_attr(file_id,label,desc,&dbuf,"??");
      
      sprintf(label, "%s/Cont_%s_%s_ext",subpath,Emlines_name[i],comp);
      sprintf(desc, "Extincted continuum luminosity, units of 10^40 erg s-1 h-2 A-1.");
      status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_FLOAT, 
				c_arr_ext);
      status = set_float_variable_attr(file_id,label,desc,&dbuf,"??");
  
	  if (ised ==1 || ised ==2)
	  {    
      	if (i ==0) {
			sprintf(label, "%s/q_avg_%s",subpath,comp);
			sprintf(desc, "Average ionization parameter from %s component, units of cm/s",comp);
			status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_FLOAT, 
				  buffer3);
			status = set_float_variable_attr(file_id,label,desc,&dbuf,"??");

			if (LineMod == 2)
			{	
				sprintf(label, "%s/q_low_%s",subpath,comp);
  				sprintf(desc, "Lowest ionization parameter from %s component, units of cm/s",comp);
				status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_FLOAT, buff4);
				status = set_float_variable_attr(file_id,label,desc,&dbuf,"??");
	
				sprintf(label, "%s/q_high_%s",subpath,comp);
	  			sprintf(desc, "Highest ionization parameter from %s component, units of cm/s",comp);
				status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_FLOAT, buff5);
				status = set_float_variable_attr(file_id,label,desc,&dbuf,"??");
			}	
		}
	/*
	  sprintf(label, "%s/NLyc_%s",subpath,comp);
	  sprintf(desc, "Number rate of Lyman Continuum photons from %s component [10^-30 s-1]",comp);
	  status = H5LTmake_dataset(file_id, label, 2, dims, H5T_NATIVE_DOUBLE, 
	  nlyc_arr);
	  status = set_float_variable_attr(file_id,label,desc,&dbuf,"??");
	*/
	
      }
    }
  }
  
  status = H5Gclose(subgroup_id);
  printf("Emission Lines and Magnitudes computed for %d galaxies\n",
	 TotNumGalA); fflush(stdout);
  
  free(buffer);
  free(nlyc_arr);
  free(buffer2);
  free(buffer3);
  free(c_arr);
  free(c_arr_ext);
  free(buff4);	
  free(buff5);

#ifdef OUTPUTSED
  free(buff6);
  free(buff7);
  free(buffW);
#endif
  
  return status;
}
#endif
