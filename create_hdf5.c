/* Routine to create an hdf5 file containing both the input lightcone catalogue and the magnitudes computed
*/

#include <stdlib.h>
#include <stdio.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

int main(int argc, char *argv[])
{
	char outfile[MAXLEN],lightconefile[MAXLEN],magfile[MAXLEN];
	char buf[MAXLEN];
	hid_t file_id;

	strcpy(lightconefile,argv[1]);
	strcpy(magfile,argv[2]);
	strcpy(outfile,argv[3]);


	if (!(flc = fopen(LightConeFile,"r")))
	{
		printf("Unable to read file \n%s\n",lightconefile);
		exit(1);
	}

	Gal = (props *) malloc(MAXGALS * sizeof(props));
	printf("Allocating %g GB for Lightcone galaxies\n",MAXGALS*sizeof(props)/1024.0/1024.0/1024.0);
//	Header of Lightcone file	
	fgets(buf,MAXLEN,flc);

	j = 0;
	while(!feof(flc))
	{
		fscanf(flc,"%d %d %d %f %f %f %f %f %f %f %f %f %f",\
		&Gal[j].snapshot,&Gal[j].id_lc,&Gal[j].id_snap,&Gal[j].x,&Gal[j].y,&Gal[j].z,\
		&Gal[j].rho,&Gal[j].theta,&Gal[j].phi,&Gal[j].redshift,&Gal[j].v_r);
	
		j++;
	}
	NTotGals = j-1;
	fclose(flc);
	
	if (!(fout = fopen(magfile,"r")))	
	{
		printf("Unable to open file %s\n",magfile);
		exit(1);
	}

	fgets(buf,MAXLEN,fout);
	fgets(buf,MAXLEN,fout);
	j = 0;
	while(!feof(fout))
	{
		fscanf(fout,"%d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",\
		&Gal[j].idLC,&Gal[j].		
	
	


	file_id = H5Fcreate(outfile,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
	
	
