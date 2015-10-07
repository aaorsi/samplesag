#include <stdlib.h>
#include <stdio.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "allvars.h"
#include "proto.h"

void write_hdf5(props *Gal, float * zapp, float *zspace, float * Mags,char *OutFile_)
{
	long j,i,id;
	hid_t file_id;
	herr_t status;
	time_t currtime;
	char label[256];
    struct tm *loctime;    
	char OutFile[256];

	sprintf(OutFile,"%s.hdf5",OutFile_);

	file_id = H5Fcreate(OutFile,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
	if (file_id < 0)
	{
		printf("Unable to open file %s\n",OutFile);
		exit(1);
	}

	status = H5LTset_attribute_string(file_id,"/","Lightcone","LSST deep-drilling");

	currtime = time(NULL);
	loctime = localtime(&currtime);
	status = H5LTset_attribute_string(file_id, "/", "Created",
	asctime(loctime));

	status = H5LTset_attribute_long(file_id,"/","Total number of galaxies",&NTotGals,1);
	
	status = H5LTset_attribute_string(file_id,"/","Questions, comments, thanks to","Alejandra Munoz (amma.19@gmail.com)\n Alvaro Orsi (aaorsi@astro.puc.cl)");


	printf("Writing HDF5 file... \n");
	hsize_t dim1[2] = {NTotGals,1};
  	double dbuf = 0.0;
	float *buff1d;
	int *bufi1d;

	buff1d = (float *) malloc(NTotGals * sizeof(float));
	bufi1d = (int *) malloc(NTotGals * sizeof(int));

	for (j =0; j<NTotGals; j++)
		bufi1d[j] = Gal[j].id_snap;
	
	sprintf(label,"/sag_id");

	status = H5LTmake_dataset(file_id,label,2,dim1,H5T_NATIVE_INT,bufi1d);	
	status = set_float_variable_attr(file_id,label,"Galaxy id in semi-analytical model",&dbuf," ");

	for (j =0; j<NTotGals; j++)
		buff1d[j] = Gal[j].redshift;

	sprintf(label,"/redshift_realspace");

	status = H5LTmake_dataset(file_id,label,2,dim1,H5T_NATIVE_FLOAT,buff1d);
	status = set_float_variable_attr(file_id,label,"Galaxy redshift, real space",&dbuf," ");


	sprintf(label,"/redshift_zspace");

	status = H5LTmake_dataset(file_id,label,2,dim1,H5T_NATIVE_FLOAT,zapp);
	status = set_float_variable_attr(file_id,label,"Galaxy redshift including peculiar motion",&dbuf," ");

	for (j =0; j<NTotGals; j++)
		buff1d[j] = Gal[j].rho;

	sprintf(label,"/rho_real");
	status = H5LTmake_dataset(file_id,label,2,dim1,H5T_NATIVE_FLOAT,buff1d);
	status = set_float_variable_attr(file_id,label,"comoving distance, real space",&dbuf,"Mpc/h");
	
	sprintf(label,"/rho_zspace");
	status = H5LTmake_dataset(file_id,label,2,dim1,H5T_NATIVE_FLOAT,zspace);
	status = set_float_variable_attr(file_id,label,"comoving distance, redshift space",&dbuf,"Mpc/h");
/*
	for (j =0; j<NTotGals; j++)
		buff1d[j] = Gal[j].x;

	sprintf(label,"/x");
	status = H5LTmake_dataset(file_id,label,2,dim1,H5T_NATIVE_FLOAT,buff1d);
	status = set_float_variable_attr(file_id,label,"spatial coordinate",&dbuf,"Mpc/h");

	for (j =0; j<NTotGals; j++)
		buff1d[j] = Gal[j].y;

	sprintf(label,"/y");
	status = H5LTmake_dataset(file_id,label,2,dim1,H5T_NATIVE_FLOAT,buff1d);
	status = set_float_variable_attr(file_id,label,"spatial coordinate",&dbuf,"Mpc/h");

	for (j =0; j<NTotGals; j++)
		buff1d[j] = Gal[j].z;

	sprintf(label,"/z");
	status = H5LTmake_dataset(file_id,label,2,dim1,H5T_NATIVE_FLOAT,buff1d);
	status = set_float_variable_attr(file_id,label,"spatial coordinate",&dbuf,"Mpc/h");



	for (j =0; j<NTotGals; j++)
		buff1d[j] = Gal[j].theta;

	sprintf(label,"/theta");
	status = H5LTmake_dataset(file_id,label,2,dim1,H5T_NATIVE_FLOAT,buff1d);
	status = set_float_variable_attr(file_id,label,"polar angle",&dbuf,"rad");

	for (j =0; j<NTotGals; j++)
		buff1d[j] = Gal[j].phi;

	sprintf(label,"/phi");
	status = H5LTmake_dataset(file_id,label,2,dim1,H5T_NATIVE_FLOAT,buff1d);
	status = set_float_variable_attr(file_id,label,"azimuthal angle",&dbuf,"rad");

*/

	for (j =0; j<NTotGals; j++)
		buff1d[j] = Gal[j].ra;

	sprintf(label,"/ra");
	status = H5LTmake_dataset(file_id,label,2,dim1,H5T_NATIVE_FLOAT,buff1d);
	status = set_float_variable_attr(file_id,label,"right ascention",&dbuf,"");

	for (j =0; j<NTotGals; j++)
		buff1d[j] = Gal[j].dec;

	sprintf(label,"/dec");
	status = H5LTmake_dataset(file_id,label,2,dim1,H5T_NATIVE_FLOAT,buff1d);
	status = set_float_variable_attr(file_id,label,"declination",&dbuf,"");

	for (j =0; j<NTotGals; j++)
		buff1d[j] = Gal[j].v_r;

	sprintf(label,"/v_r");
	status = H5LTmake_dataset(file_id,label,2,dim1,H5T_NATIVE_FLOAT,buff1d);
	status = set_float_variable_attr(file_id,label,"radial velocity",&dbuf,"km/s");

	for (i = 0; i<NMags; i++)
	{
		sprintf(label,"/mag_%s_realspace",NameMags[i]);
	
		for (j =0; j<NTotGals; j++)
		{
			id = i + (NMags*2)*j;	
			buff1d[j] = Mags[id];
		}
		
		status = H5LTmake_dataset(file_id,label,2,dim1,H5T_NATIVE_FLOAT,buff1d);
		status = set_float_variable_attr(file_id,label,Idinfo[i],&dbuf,"magnitude, AB system, observer-frame, real space");
	
		sprintf(label,"/mag_%s_zspace",NameMags[i]);
	
		for (j =0; j<NTotGals; j++)
		{
			id = i+NMags + (NMags*2)*j;
			buff1d[j] = Mags[id];
		}
		
		status = H5LTmake_dataset(file_id,label,2,dim1,H5T_NATIVE_FLOAT,buff1d);
		status = set_float_variable_attr(file_id,label,Idinfo[i],&dbuf,"magnitude, AB system, observer-frame, redshift space");
	
	}
	

    status = H5Fclose(file_id);	
	
	free(buff1d);
	free(bufi1d);
}

int set_float_variable_attr(hid_t file_id, char *label, char *desc, 
                double *unit, char *comment)
{
  herr_t status;

  status = H5LTset_attribute_string(file_id, label, "Description", desc);
  if (status < 0) {
    fprintf(stderr, "Error (HDF5): cannot create attribute 'Description'"
        " - Exit\n"); fflush(stderr);
    exit(EXIT_FAILURE);
  }
  
  status = H5LTset_attribute_double(file_id, label, "Unit", unit, 1); 
  if (status < 0) {
    fprintf(stderr, "Error (HDF5): cannot create attribute 'Unit'"
        " - Exit\n"); fflush(stderr);
    exit(EXIT_FAILURE);
  }

  status = H5LTset_attribute_string(file_id, label, "Comment", comment);
  if (status < 0) {
    fprintf(stderr, "Error (HDF5): cannot create attribute 'Comment'"
        " - Exit\n"); fflush(stderr);
    exit(EXIT_FAILURE);
  }

  return EXIT_SUCCESS;
}


