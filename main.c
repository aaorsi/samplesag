// Read lightcone file, read SAG file and extract SED of every galaxy, 
// and compute the observer-frame magnitudes from a list.

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

		
//	write_fits(Gal,zapp,zspace,Mags,OutFile);
//	exit(0);

/*#ifdef COMOVFILE

	printf("Creatinf look-up table %s\n",rho_zFile);
	if (!(flc = fopen(rho_zFile,"w")))
	{
		printF("Unable to create file %s\n",rho_zFile);
		exit(1);
	}

	nTab = 1000;
	int zmin = 0;
	int zmax = 20.0;
	float bin = (zmax - zmin)/(nTab -1.0);
	
	for (i = 0; i<nTab; i++)
	{
		
		zTab[i] = i*bin + zmin;
		rTab[i] = comoving_distance(zTab[i]);
		
		fprintf(flc,"%f %f",&zTab[i],&rTab[i]);
	}
	fclose(flc);
	printf("File %s created. Now run again, without compiling flag COMOVFILE\n",rho_zFile);
	exit(0); 	

#endif		
*/
	

	printf("Reading look-up table %s\n",rho_zFile);	
	if (!(flc = fopen(rho_zFile,"r")))
	{
		printf("Unable to read file %s\n",rho_zFile);
		exit(1);
	}

	j = 0;
	while(!feof(flc))
	{
		fscanf(flc,"%f %f",&zTab[j],&rTab[j]);
		j++;
	}
	
	nTab = j;
	
	printf("File %s read\n",rho_zFile);
	fclose(flc);
	
	printf("Reading light cone file \n%s\n",LightConeFile);	
	
	if (!(flc = fopen(LightConeFile,"r")))
	{
		printf("Unable to read file \n%s\n",LightConeFile);
		exit(1);
	}
	
	Gal = (props *) malloc(MAXGALS * sizeof(props));
	printf("Allocating %g GB for Lightcone galaxies\n",MAXGALS*sizeof(props)/1024.0/1024.0/1024.0);
//	Header of Lightcone file	
	fgets(buf,MAXLEN,flc);

	j = 0;
	while(!feof(flc))
	{
//		printf("%d\n",j);
		fscanf(flc,"%d %d %d %f %f %f %f %f %f %f %f %f %f",\
		&Gal[j].snapshot,&Gal[j].id_lc,&Gal[j].id_snap,&Gal[j].x,&Gal[j].y,&Gal[j].z,\
		&Gal[j].rho,&Gal[j].theta,&Gal[j].phi,&Gal[j].ra,&Gal[j].dec,&Gal[j].redshift,&Gal[j].v_r);

/*	
		fscanf(flc,"%d %ld %d %g %g %g %g %g %g %g %g %g %g",\
		&Gal[j].snapshot,&Gal[j].id_lc,&Gal[j].id_snap,&junk,&junk,&junk,\
		&Gal[j].rho,&junk,&junk,&junk,&junk,&Gal[j].redshift,&Gal[j].v_r);
*/	
		j++;
	}
	printf("%ld galaxies read succesfully\n",j);

	NTotGals = j-1;
	fclose(flc);

	NTotGals = 1000000;

	Filterlist = (float *) malloc(NFILTERMAX * NFILTERL * 2 * sizeof(float));
  	Nsteps_filter = (int *) malloc(NFILTERMAX * sizeof(int));

	get_filterinfo(FilterInfo);
	read_filterfile(FilterData);
	get_Vega_AB(VegaFile);			

/*	if (!(fout = fopen(OutFile,"w")))	
	{
		printf("Unable to create file %s\n",OutFile);
		exit(1);
	}

	fprintf(fout,"#Number of filters: %d\n",NMags);
	fprintf(fout,"# idLightcone idSnapFile zreal zapp rho rho_app Mags(real) Mags(app)\n");
*/

	 	

	currsnap = -1;	// This helps determining whether its necessary to open a new snapshot file

	zapp 	 = (float *) malloc(NTotGals * sizeof(float));
	zspace 	 = (float *) malloc(NTotGals * sizeof(float));

	Mags = (float *) malloc(NTotGals * NMags*2 * sizeof(float));

	taueff = (float *) malloc(NBinSed * sizeof(float));
	spec = (float *) malloc(NBinSed * sizeof(float));

/*	for(j = 0;j<NMags; j++);
	{		
		Mags[j]   = (float *) malloc(NTotGals * sizeof(float));
		Mags[j+6] = (float *)malloc(NTotGals * sizeof(float));
	}
*/

	for (j = 0; j<NTotGals; j++)
//	for (j = 0; j<1e6; j++)
	{
		snap = Gal[j].snapshot;	
		redshift = Gal[j].redshift;
		E_z = sqrt( pow((1 + redshift),3)*OmegaM + (1 - OmegaM));
	    ExpRate = 1./(1 + redshift);
		
		zspace[j] = Gal[j].rho + Gal[j].v_r/(ExpRate*H0*E_z);
		if (zspace[j] < 0.0) 
			zspace[j] = 0.0;

// 		Determine apparent redshift (in redshift space)				
//		by using look-up table

		k = 0;
		while(rTab[k] <= zspace[j])
		{
			rlow = rTab[k];
			zlow = zTab[k];
			k++;
			if (k == nTab)
			{
				printf("zspace of gal. %d is out of look-up table range. Check!\n",j);
				dprintf(zspace);
				dprintf(Gal[j].redshift);
				dprintf(Gal[j].rho);
				exit(1);
			}
		}
		rhigh = rTab[k];
		zhigh = zTab[k];

		zapp[j] = (zspace[j] - rlow) / (rhigh - rlow) * (zhigh - zlow) + zlow;

/*		if (zspace < Gal[j].rho)
		{
			dprintf(zspace);
			dprintf(Gal[j].rho);
			dprintf(redshift);
			dprintf(zapp);
		}
*/

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

				
				sprintf(hdf5file,"%s/%s/Subgalaxies/gal_itf_%s_SAG3_%s.hdf5",RootDir,strbox,strsnap,RunName);
			
				if (!(fbox = fopen(hdf5file,"r")))
				{
//					fclose(fbox);
					continue;
				}
				fclose(fbox);	
//				printf("Reading %s\n",hdf5file);

				read_saghdf5(hdf5file,SnapGal,wlambda,ig);

//		Filtering by stellar mass
/*				for (k = ig;k<ig+NGals;k++)
				{
					if(SnapGal[k].StellarMass >= MinStellarMass/1.0e10)
					{
						idFiltered[i] = k;
						i++;
					}	
				}
*/
				ig += NGals;
			}
			
			printf("Total number of filtered galaxies: %ld\n",ig);
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
//			spec = SnapGal[idSnap].Sed;	
			for (_i = 0;_i <NBinSed; _i++)
				spec[_i] = SnapGal[idSnap].Sed[_i];
		}

//		fprintf(fout,"%d\t%d\t%f\t%f\t%f\t%f\t",Gal[j].id_lc,idSnap,redshift,zapp,Gal[j].rho,zspace);				
// 		Magnitudes in real-space, observer-frame

		for (k = 0; k<NMags; k++)
		{
			id = k + (NMags*2)*j;
			Mags[id] = compute_mag_filter(idMag[k],Nsteps_filter[idMag[k]],\
			wlambda,spec,1,redshift,0,1);

/*			if (k <NMags-1)
				fprintf(fout,"%f ",mag);
			else
				fprintf(fout,"%f\t",mag);
*/
		}
// 			Magnitudes in redshift-space, observer-frame

		for (k = 0; k<NMags; k++)
		{
			id = k+NMags + (NMags*2)*j;
			Mags[id] = compute_mag_filter(idMag[k],Nsteps_filter[idMag[k]],\
			wlambda,spec,1,zapp[j],0,1);
/*
			if (k <NMags-1)
				fprintf(fout,"%f ",mag);
			else
				fprintf(fout,"%f\n",mag);
*/
		}
	}

	if (strcmp(OutFormat,"HDF5") == 0 || \
		strcmp(OutFormat,"BOTH") == 0)
	{
		write_hdf5(Gal,zapp,zspace,Mags,OutFile);
	}
	
	if (strcmp(OutFormat,"FITS") == 0 || \
		strcmp(OutFormat,"BOTH") == 0)
	{
		write_fits(Gal,zapp,zspace,Mags,OutFile);
	}
	

//	fclose(fout);
	free(Gal);
	free(SnapGal);

	free(wlambda);

	free(zapp);
	free(zspace);


//	free(FilterInfo);
//	free(FilterData);
	
//	free(Idinfo);

}


