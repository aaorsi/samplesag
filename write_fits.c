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

void write_fits(props *Gal, float *zapp, float *zspace, float *Mags,char *OutFile_)
{

    float *buff1d;
    int *bufi1d;
	char OutFile2[128];
	char label[256];
	
	sprintf(OutFile2,"%s.fit",OutFile_);

    buff1d = (float *) malloc((NTotGals+1) * sizeof(float));
    bufi1d = (int *) malloc((NTotGals+1) * sizeof(int));

    fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
    int status, hdutype;
    long firstrow, firstelem;
	long j,i,id;
		
    char extname[] = "Lightcone, LSST deep-drilling";           /* extension name */
	
    status = 0;

    /* open the FITS file containing a primary array and an ASCII table */

	if (fits_create_file(&fptr, OutFile2, &status)) printerror(status);

//    if ( fits_open_file(&fptr, OutFile2, READWRITE, &status) ) 
//      printerror( status );

//    if ( fits_movabs_hdu(fptr, 2, &hdutype, &status) ) /* move to 2nd HDU */
//         printerror( status );

	int tfields = 9 + 2*NMags;
	long nrows  = NTotGals;

	char *ttype[tfields];

	ttype[0]  = "sag_id";
	ttype[1] = "snapshot";
	ttype[2] = "redshift_realspace";
	ttype[3] = "redshift_zspace";
	ttype[4] = "rho_real";
	ttype[5] = "rho_zspace";
	ttype[6] = "ra";
	ttype[7] = "dec";
	ttype[8] = "v_r";
/*
	strcpy(ttype[0],"sag_id");
	strcpy(ttype[1],"snapshot");
	strcpy(ttype[2],"redshift_realspace");
	strcpy(ttype[3],"redshift_zspace");
	strcpy(ttype[4],"rho_real");
	strcpy(ttype[5],"rho_zspace");
	strcpy(ttype[6],"ra");
	strcpy(ttype[7],"dec");
	strcpy(ttype[8],"v_r");
*/

	char FullName[NMags*2][32];

	for (i = 0; i<NMags; i++)
	{
		sprintf(FullName[i],"Magnitude %s real-space",NameMags[i]);
		ttype[9+i] = FullName[i];
		sprintf(FullName[NMags+i],"Magnitude %s redshift-space",NameMags[i]);
		ttype[NMags+9+i] = FullName[NMags+i];
	}

	char *tform[tfields];
	tform[0] = "1J";
	for (i = 0; i< 8 + 2*NMags; i++)
		tform[1+i] = "1E";

	char *tunit[tfields];
	for (i = 0; i<9+2*NMags; i++)
		tunit[i] = "\0";

//	char *ttype[] = {"sag_id", "snapshot", "redshift_realspace", "redshift_zspace", "rho_real","ra", "dec", "v_r"}

    /* append a new empty binary table onto the FITS file */
    if ( fits_create_tbl( fptr, BINARY_TBL, nrows, tfields, ttype, tform,
                tunit, extname, &status) )
         printerror( status );

    firstrow  = 1;  /* first row in table to write   */
    firstelem = 1;  /* first element in row  (ignored in ASCII tables) */

	printf("Writing FITS file %s\n",OutFile2);
	fflush(stdout);

    for (j =0; j<NTotGals; j++)
	bufi1d[j+1] = Gal[j].id_snap;
	
    fits_write_col(fptr, TLONG, 1, firstrow, firstelem, nrows, bufi1d,
                   &status);
	
    for (j =0; j<NTotGals; j++)
	bufi1d[j+1] = Gal[j].snapshot;

    fits_write_col(fptr, TLONG, 2, firstrow, firstelem, nrows, bufi1d,
                   &status);

    for (j =0; j<NTotGals; j++)
	buff1d[j+1] = Gal[j].redshift;

    fits_write_col(fptr, TFLOAT, 3, firstrow, firstelem, nrows, buff1d,
                   &status);

    fits_write_col(fptr, TFLOAT, 4, firstrow, firstelem, nrows, zapp,
                   &status);

    for (j =0; j<NTotGals; j++)
	buff1d[j+1] = Gal[j].rho;

    fits_write_col(fptr, TFLOAT, 5, firstrow, firstelem, nrows, buff1d,
                   &status);

    fits_write_col(fptr, TFLOAT, 6, firstrow, firstelem, nrows, zspace,
                   &status);

    for (j =0; j<NTotGals; j++)
	buff1d[j+1] = Gal[j].ra;

    fits_write_col(fptr, TFLOAT, 7, firstrow, firstelem, nrows, buff1d,
                   &status);
	
    for (j =0; j<NTotGals; j++)
	buff1d[j+1] = Gal[j].dec;

    fits_write_col(fptr, TFLOAT, 8, firstrow, firstelem, nrows, buff1d,
                   &status);

    for (j =0; j<NTotGals; j++)
	buff1d[j+1] = Gal[j].v_r;

    fits_write_col(fptr, TFLOAT, 9, firstrow, firstelem, nrows, buff1d,
                   &status);

	for (i = 0; i<NMags; i++)
	{
		for (j =0; j<NTotGals; j++)
		{
			id = i + (NMags*2)*j;	
			buff1d[j+1] = Mags[id];
		}
		
    		fits_write_col(fptr, TFLOAT, 10+i, firstrow, firstelem, nrows, buff1d,
                   &status);
	
		for (j =0; j<NTotGals; j++)
		{
			id = i+NMags + (NMags*2)*j;
			buff1d[j+1] = Mags[id];
		}
		

    		fits_write_col(fptr, TFLOAT, NMags+10+i, firstrow, firstelem, nrows, buff1d,
                   &status);
	
	}

    if ( fits_close_file(fptr, &status) )       /* close the FITS file */
         printerror( status );
}

void printerror( int status)
{
    /*****************************************************/
    /* Print out cfitsio error messages and exit program */
    /*****************************************************/


    if (status)
    {
       fits_report_error(stderr, status); /* print error report */

       exit( status );    /* terminate the program, returning error status */
    }
    return;
}






