#include <stdlib.h>
#include <stdio.h>
#include <hdf5.h>
#include <hdf5_hl.h>

#define MAXLEN 500

#include "allvars.h"
#include "proto.h"

void read_gadget(particle *part, char *fname)
{
  
/*  if(argc != 2)
    {
      printf("Usage: read_gadget <filename>\n");
      exit(1);
    }

//   Get name of input file 
  char fname[MAXLEN];
  strncpy(fname, argv[1], MAXLEN);

*/
  
  /* Open the file */
  hid_t file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if(file_id < 0)
    {
      printf("Unable to open file!\n");
      exit(1);
    }

  /* Get number of particles in this file */
  unsigned int npfile[6];
  if(H5LTget_attribute_uint(file_id, "Header", "NumPart_ThisFile", npfile) < 0)
    {
      printf("Failed to read number of particles!\n");
      exit(1);
    }
  
  /* Allocate storage */
  float *pos = malloc(3*sizeof(float)*npfile[1]);
  float *vel = malloc(3*sizeof(float)*npfile[1]);

  /* Read positions and velocities */
  if(H5LTread_dataset_float(file_id, "PartType1/Coordinates", pos) < 0)
    {
      printf("Unable to read positions!\n");
      exit(1);   
    }
  if(H5LTread_dataset_float(file_id, "PartType1/Velocities", vel) < 0)
    {
      printf("Unable to read velocities!\n");
      exit(1);   
    }

  /* Write out particle data */
  int i;
  for(i=0;i<npfile[1];i+=1)
  {
    part[i].x = pos[3*i+0];  
    part[i].y = pos[3*i+1];  
    part[i].z = pos[3*i+2];  

    part[i].vx = vel[3*i+0];
    part[i].vy = vel[3*i+1];
    part[i].vz = vel[3*i+2];
  }

 //   printf("%14.6f %14.6f %14.6f %14.6f %14.6f %14.6f\n",
//	   pos[3*i+0], pos[3*i+1], pos[3*i+2],
//	   vel[3*i+0], vel[3*i+1], vel[3*i+2]);
  
  free(pos);
  free(vel);

  /* Close file */
  H5Fclose(file_id);

}
