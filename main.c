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
	int currsnap,snap,idSnap;
	float E_z,ExpRate,rlow,rhigh,zlow,zhigh,mag;
	float *wlambda;	
	float junk;
	float redshift, *zapp, *zspace;
	float *Mags;
	float *taueff;	
	float *spec;	

	strcpy(RootDir,argv[1]);
	strcpy(RunName,argv[2]);
	strcpy(rho_zFile,argv[3]);
	strcpy(OutFile	,argv[4]);
	strcpy(FilterInfo,argv[5]);
	strcpy(FilterData,argv[6]);
	strcpy(VegaFile,argv[7]);
	strcpy(LightConeFile, argv[8]);
	NBinSed = atoi(argv[9]);
	NMags   = atoi(argv[10]);
	for (i = 0;i<NMags; i++)
		idMag[i] = atoi(argv[11+i]);
	NBoxes	= atoi(argv[NMags+11]);
	MinStellarMass = atof(argv[NMags+12]);
	IGMAtt = atoi(argv[NMags+13]);	
	strcpy(OutFormat,argv[NMags+14]);
	for (i = 0; i<NMags; i++)
		sprintf(NameMags[i],"%s",argv[NMags+15+i]);
//		strcpy(NameMags[i],argv[NMags+15+i]);

	printf("INPUT PARAMETERS\n");
	dprints(RootDir);
	dprints(rho_zFile);
	dprints(OutFile);
	dprints(LightConeFile);
	dprinti(NBinSed);
	dprinti(NBoxes);
	dprinti(NMags);
	for (i = 0;i<NMags;i++)
		printf("Band %d: %s\n",idMag[i],&NameMags[i]);
//		dprinti(idMag[i]);
	dprintf(MinStellarMass);

	wlambda = (float *) malloc(NBinSed * sizeof(float));


	NTotGals = 1000000;

	Filterlist = (float *) malloc(NFILTERMAX * NFILTERL * 2 * sizeof(float));
	Nsteps_filter = (int *) malloc(NFILTERMAX * sizeof(int));

	get_filterinfo(FilterInfo);
	read_filterfile(FilterData);
	get_Vega_AB(VegaFile);			

	currsnap = -1;	// This helps determining whether its necessary to open a new snapshot file

	Mags = (float *) malloc(NTotGals * NMags*2 * sizeof(float));
	taueff = (float *) malloc(NBinSed * sizeof(float));
	spec = (float *) malloc(NBinSed * sizeof(float));

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
			read_saghdf5(hdf5file,SnapGal,wlambda,ig);

			ig += NGals;
		}
			
	}	
		
	NGals = ig;

		//idSnap identifies the galaxy within the SnapGal structure (the unfiltered array)
	idSnap = Gal[j].id_snap;

	if (IGMAtt)
	{
		igm_attenuation(taueff,wlambda,redshift,idSnap);	
		spec = SnapGal[idSnap].Sed;	
		for (_i = 0;_i <NBinSed; _i++)
			spec[_i] = SnapGal[idSnap].Sed[_i] * exp(-1*taueff[_i]);
	}
	else
	{
		for (_i = 0;_i <NBinSed; _i++)
			spec[_i] = SnapGal[idSnap].Sed[_i];
	}

// 		Magnitudes in real-space, observer-frame

	for (k = 0; k<NMags; k++)
	{
		id = k + (NMags*2)*j;
		Mags[id] = compute_mag_filter(idMag[k],Nsteps_filter[idMag[k]],\
		wlambda,spec,1,redshift,0,1);

	}

	if (strcmp(OutFormat,"HDF5") == 0 || \
		strcmp(OutFormat,"ALL") == 0)
	{
		write_hdf5(Gal,zapp,zspace,Mags,OutFile);
	}

	if (strcmp(OutFormat,"ASCII") == 0 || \
		  strcmp(OutFormat,"ALL") == 0)
		write_ascee(Gal, zapp,zspace,Mags,OutFile);
	

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


