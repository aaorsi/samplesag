// samplesag: Extracts galaxy properties from a SAG run, and computes magnitudes, luminosities and dust extinction
//
// Alvaro Orsi, aaorsi@cefca.es
//


#include <stdlib.h>
#include <stdio.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "allvars.h"
#include "proto.h"
#include "fitsio.h"
#include <gsl/gsl_rng.h>

// Read IDs   

int main(int argc, char *argv[])
{

  char RootDir[MAXLEN];
  char OutFile[MAXLEN];
  char LightConeFile[MAXLEN];
  char rho_zFile[MAXLEN];
  char FilterInfo[MAXLEN];
  char FilterData[MAXLEN];
  char VegaFile[MAXLEN];
  char buf[2*MAXLEN];
  char hdf5file[MAXLEN],RunName[MAXLEN];
  char OutFormat[256];
  int i,nTab,id,_i;
  int IGMAtt;
  float rTab[MAXZ],zTab[MAXZ];
  int *idFiltered;
  FILE *flc, *fout, *fbox;
  long j,k,ig;
  char strsnap[3];
  char strbox[3];
  int ibox;
  props *Gal;
  sagobj *SnapGal;  
  cond conditions[10];
  int currsnap,snap,idSnap;
  float E_z,ExpRate,rlow,rhigh,zlow,zhigh,mag;
  float *wlambda; 
  float junk;
  float redshift, *zapp, *zspace;
  float *Mags;
  float *taueff;  
  float *spec;  
  float *incl;
  const gsl_rng_type * T;


  gsl_rng *r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  gsl_rng_set(r, seed0);

  i = 1
  strcpy(RootDir,argv[i++]);
  strcpy(RunName,argv[i++]);
  nsnap = atoi(argv[i++]);
  int SnapList[nsnap];
  for (j = 0; j<nsnap ; j++)
    SnapList[i] = atoi(argv[i + j])
  
  i += j

  nvol = atoi(argv[i++])
  iv0  = atoi(argv[i++])
  ivf  = atoi(argv[i++])

  strcpy(OutputFormat,argv[i++]);
  strcpy(FilterInfo,argv[i++]);
  strcpy(FilterData,argv[i++]);
  strcpy(VegaFile,argv[i++]);

  NMags   = atoi(argv[i++]);

  int idMag[NMags];
  const char *NameMags[NMags];
  int ABArr[NMags];

  for (j = 0;j<NMags; j++)
    idMag[i] = atoi(argv[i+j]);

  i += j;
  for (j = 0; j<NMags; j++)
    sprintf(NameMags[j],"%s",argv[i + j]);

  i += j;
  for (j = 0; j<NMags; j++)
    ABArr[j] = atoi(argv[i + j]);

  i += j;

  strcpy(propcheck,argv[i++]);

//  Check consistency
  if (propcheck != "props")
  {
    printf("ERROR: main.c: When reading arguments, \nproperty list expected, got this instead: %s\n",propcheck);
    exit(1);
  }

  PropLength = argc - i;
 
  k = 0;

  for (j = 0; j<PropLength; j++)
  {
    // For now there are only two conditions implemented: gt and lt.
    if (strcmp(argv[i + j],"gt") == 0)
    {
      conditions[k].prop = PropArr[j-1];
      conditions[k].val  = atof(argv[i + j +1]);
      conditions[k].type = 1;
      i += 2;
      PropLength -= 1;
      k++;
    }
    else if (strcmp(argv[i + j],"lt") == 0)
    {
      conditions[k].prop = PropArr[j-1];
      conditions[k].val  = atof(argv[i + j +1]);
      conditions[k].type = -1;
      i += 2;
      PropLength -= 1;
      k++;
    }
    else
      sprintf(PropArr[j],"%s",argv[i + j]);
  } 

  nconditions = k;

  print_input();
  printf("SAMPLESAG\n\n");
  printf("INPUT PARAMETERS\n");

  for (i = 0;i<NMags;i++)
    printf("Band %d: %s\n",idMag[i],&NameMags[i]);
//    dprinti(idMag[i]);

  wlambda = (float *) malloc(SED_MAXW* sizeof(float));

  Filterlist = (float *) malloc(NFILTERMAX * NFILTERL * 2 * sizeof(float));
  Nsteps_filter = (int *) malloc(NFILTERMAX * sizeof(int));

  get_filterinfo(FilterInfo);
  read_filterfile(FilterData);
  get_Vega_AB(VegaFile);      

  currsnap = -1;  // This helps determining whether its necessary to open a new snapshot file

  Mags = (float *) malloc(NTotGals * NMags*2 * sizeof(float));
  incl = (float *) malloc(NTotGals * sizeof(float));
  taueff = (float *) malloc(NBinSed * sizeof(float));
  spec = (float *) malloc(NBinSed * sizeof(float));

  icount = 0
  if (currsnap != snap)
  {

    printf("Reading snapshot %d\n",snap);
    fflush(stdout);

    if (currsnap != -1)
    {
      for (i = 0;i<NGals; i++)
        free(SnapGal[i].Sed);
      free(SnapGal);
    }

    SnapGal = (sagobj *) malloc(MAXGALSNAP * sizeof(sagobj)); 
    currsnap = snap;
  
// Get Data from HDF5 snapshot file

    if (snap < 10)
      sprintf(strsnap,"00%d",snap);
    else
      sprintf(strsnap,"0%d",snap);

    ig = 0l;
    i = 0;  
    for (ibox = 1;ibox <= NBoxes; ibox++)
    {

      if (ibox < 10)
        sprintf(strbox,"00%d",ibox);
      else
        sprintf(strbox,"0%d",ibox);

        
      sprintf(hdf5file,"%s/%s/Subgalaxies/gal_itf_%s_SAG3_%s.hdf5",\
              RootDir,strbox,strsnap,RunName);
      
      if (!(fbox = fopen(hdf5file,"r")))
        continue;
      fclose(fbox); 
      
      read_saghdf5(hdf5file,SnapGal,wlambda,ig,icount);

      ig += NGals;
    }
      
  } 
    
  NGals = ig;

// A loop over all galaxies to compute their properties

  for (id = 0; id < NGals; id++)
  {

    if (IGMAtt)
    {
      igm_attenuation(taueff,wlambda,redshift,id);  
      spec = SnapGal[id].Sed; 
      for (_i = 0;_i <NBinSed; _i++)
        spec[_i] = SnapGal[id].Sed[_i] * exp(-1*taueff[_i]);
    }
    else
    {
      for (_i = 0;_i <NBinSed; _i++)
        spec[_i] = SnapGal[id].Sed[_i];
    }
   
//  Identify what properties are required. This is only done once
    if ( id == 0)
    { 
        Magdump   = get_magsprops(PropArr);   // magnitudes
        Linesdump = get_linesprops(PropArr);  // lines
        Copydump  = get_copiedprops(PropArr); // taken from standard SAG HDF5 dumps
    }

// Copy variables that dont require post-processing 

    if (Copydump == 1)
    {
      for (k = 0; k< NCopydump; k++)
        PropVal[id][CopyP[k]] = get_vars();
    }

    // Compute magnitudes
    if (Magdump == 1)
    { 
      for (k = 0; k< NMagdump; k++)
        Propval[id][MagP[k]] = compute_mag_filter(idMag[k],Nsteps_filter[idMag[k]],\
        wlambda,spec,1,redshift,0,1);
    }

    // Compute line luminosities
    if (Linesdump == 1)
    {
      nlyc_arr[id] = compute_lyc(wlambda,SnapGal[id].Sed); 
      for (k = 0; k < NLinedump; k++)
        Propval[id][LineP[k]] = integ_line(zcold[id],nlyc_arr[id],id, ised);
    }

 
  }

  if (strcmp(OutFormat,"HDF5") == 0 || \
    strcmp(OutFormat,"ALL") == 0)
  {
    write_hdf5(Gal,zapp,zspace,Mags,OutFile);
  }

  if (strcmp(OutFormat,"ASCII") == 0 || \
      strcmp(OutFormat,"ALL") == 0)
    write_ascii(Gal, zapp,zspace,Mags,OutFile);
  

  if (strcmp(OutFormat,"FITS") == 0 || \
    strcmp(OutFormat,"BOTH") == 0)
  {
    write_fits(Gal,zapp,zspace,Mags,OutFile);
  }


  free(Gal);
  free(SnapGal);

  free(wlambda);

  free(zapp);
  free(zspace);



}


