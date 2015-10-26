#include <stdlib.h>
#include <stdio.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <strings.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"


void read_saghdf5(char *fname, sagobj *SnapGal, float *wlambda, long id_arr, int icount)
{ 

  hid_t file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  hsize_t dims[2];
  int ng,i,j,id,k;
  float *DiscMass, *BulgeMass;
  float SMass;

  if (file_id < 0)
  {
    printf("ERROR (read_saghdf5) : Unable to open file %s\n",fname);
    exit(1);
  }

  if (H5LTget_dataset_info(file_id,"/SFR",dims,NULL,NULL) < 0)
  {
    printf("ERROR (read_saghdf5): Unable to read dimensions from SFR.\n");
    exit(1);
  }

  NGals = (size_t) (dims[0]*dims[1]);
  printf("INFO (read_saghdf5): Reading %d galaxies from file %s\n",NGals,fname);


  if (icount == 0)
  {
    if (H5LTread_dataset_float(file_id,"/SED/Magnitudes/lambda",wlambda) < 0)
    {
      printf("ERROR (read_saghdf5): unable to read lambda array. Exiting\n");
      exit(1);
    } 

    NBinSed = sizeof(wlambda)/sizeof(wlambda[0]);
  }

  Prop = (float *) malloc(NGals * sizeof(float));


// Read variables with conditions. Get ID list for galaxies satisfying conditions.

  for (ic = 0; ic< nconditions; ic++)
  {
    vars_in_subdirs(conditions[ic].prop,buf);

    if (H5LTread_dataset_float(file_id,buf,Prop) < 0)
    {
      printf("ERROR (read_saghdf5): Unable to read %s\n",Dump[k].Name);
      exit(1);
    }

    for (ig = 0; ig < NGals; ig++)
    {
      if (conditions[k].type == 1)
      { 
        if (Prop[ig] < conditions[k].val)
          gflag[ig]++;
      }
      else
      {
        if (Prop[ig] > conditions[k].val)
          gflag[ig]++;
      }
    }
  // In the end, any galaxy with gflag > 0 is not recorded.  
  }


  NGals = 0;
  for (k = 0; k< CopyDump; k++)
  {
    sprintf(location,"/%s%s",Dump[k].Dir,Dump[k].Name);
    if (H5LTread_dataset_float(file_id,location,Prop) < 0)
    {
      printf("ERROR (read_saghdf5): Unable to read %s\n",Dump[k].Name);
      exit(1);
    }

    for (i = 0;i<NGals;i++)
    {
      if (gflag[i] == 0)
        SnapGal[i+id_arr].Props[k] = Prop[i];
        NGals++;
    }
    
  }
  // Read SED only if magnitudes are requested
  if (Magdump > 0)
  {
    float * sed_arr = malloc(NGals * NBinSed * sizeof(float));

    if (H5LTread_dataset_float(file_id,"/SED/Magnitudes/Sed",sed_arr) < 0)
    {
      printf("ERROR (read_saghdf5): Unable to read sed.\n");
      exit(1);
    }

//    float *dust_ext = malloc(NGals*NBinSed * sizeof(float));
//    if (H5LTread_dataset_float(file_id,"/SED/Magnitudes/dust_ext",dust_ext) < 0)
//    {
//      printf("ERROR (read_saghdf5): Unable to read dust_ext.\n");
//      exit(1);
//    }

  }





//  No idea if this works... it must be tested.
/*  k = 0;
  for (i = 0;i < NGals; i++)
  {
    SMass = DiscMass[i] + BulgeMass[i];
    if( SMass >= Hubble_h * MinStellarMass/1.0e10)
    {
      SnapGal[id_arr+k].StellarMass = SMass;
      SnapGal[id_arr+k].Sed = (float *) malloc(NBinSed*sizeof(float));
      for (j = 0; j<NBinSed; j++)
      {
        id = j + i*NBinSed;
        SnapGal[id_arr+k].Sed[j] = dust_ext[id] * sed_arr[id];
      } 
      k++;
    }
  }
*/

  NGals = k;

  free(dust_ext);
  free(sed_arr);
  free(Prop);

  if (H5Fclose(file_id) < 0)
  {
    printf("ERROR (read_saghdf5): Error while closing file %s\n",fname);
    exit(1);
  } 

}
    
