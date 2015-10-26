#include<string.h>
#include "allvars.h"
#include "proto.h"


// Adds path to variable in case it belongs to a subdir in the HDF5 file.
void vars_in_subdirs(char *Propj, char *path)
{
  // Halo dir
  if ((strcmp(Propj,"M200c") == 0) || (strcmp(Propj,"R200c") == 0) \
     || (strcmp(Propj,"Vmax") == 0) || (strcmp(Propj,"Vpeak") == 0) \
     || (strcmp(Propj,"cNFW") == 0))
     sprintf(path,"/Halo/");

  return(1);
}


int get_magsprops(char * PropArr)
{
  int dust_check,mgal;

// This function identifies which magnitudes will be computed. Returns 1 if there are magnitudes at all, 0 otherwise.
//  ** Magnitudes **    
  dust_check = 0;
  id = 0;
  for (j = 0; j< PropLength; j++)
  { 
    // Check if PropArr is a magnitude
    if (strstr(PropArr[j],"mag_") != 0)
    {
      // Identify the magnitude
      for (k = 0; k < NMags; k++) 
      {
        sprintf(buf,"mag_%s",NameMags[k]);
        if (strstr(PropArr[j],buf) != 0)
        {
          MagList[id].id = k;
          // identify the frame
          if (strstr(PropArr[j],"o") != 0)
          {
            sprintf(buf,"mag_%so",NameMags[k]);
            oframe = 1;
          } else
          {
            sprintf(buf,"mag_%sr",NameMags[k]);
            oframe = 0;
          }
          MagList[id].oframe = oframe;
          // Check if magnitude needs dust extinction
          if (strstr(PropArr[j],"x") != 0 && dust_check != 0)
          {
            MagList[id].dust = 1;
            // Inclination angle (random) and dust extinction for each galaxy:    
            incl[id] = gsl_rng_uniform (r);
            //  B-band magnitude is necessary to compute dust extinction for this model.
            magB = compute_mag_filter(idB,Nsteps_filter[idB], wlambda,SnapGal[id].Sed,
                                      0,redshift,0,0);
            dust_ext[id] = dust_extinction(wlambda,incl[id],id,magbin,SnapGal[id].Sed,magB);
            dust_check = 1;
          }
          id ++;
        }
      }
    }
  }

  return id;
}

int get_vars(char * PropArr)
{
/* Identify which variables are going to be copied from the hdf5 dump file */
  int id = 0;
  for (j = 0; j< PropLength< j++)
  {
    Propj = PropArr[j];
    // Avoid variables starting with mag_ or l_ since those are treated separately.
    if ((strstr(Propj,"mag_") == 0) && 
        (strstr(Propj,"l_") == 0))
    {
      strcpy(Dump[id].Name,Propj);
      vars_in_subdirs(Propj, Durmp[id].Dir);

/*      // Halo dir
      if ((strcmp(Propj,"M200c") == 0) || (strcmp(Propj,"R200c") == 0) \
      || (strcmp(Propj,"Vmax") == 0) || (strcmp(Propj,"Vpeak") == 0) \
      || (strcmp(Propj,"cNFW") == 0))
        strcpy(Dump[id].Dir,"/Halo/");
*/
    id++;
    }
  }
  return id:
}




int get_linesprops(char *PropArr)
{
// Identify lines to be computed.
  int j,k;
  char *buf;
  int id = 0;

  for (j = 0; j< PropLength; j++)
  {
    Propj = PropArr[j];
    // Check for emission lines in PropArr
    if (strstr(Propj,"l_") != 0)
    {
      for (k = 0; k< NLines; k++)
      {
        if strstr(Emlines_name[k],Propj) != 0 )
        {
          Lines[id].id = k;
          if (strstr(Propj,"_x") != 0)
            Lines[id].dust = 1
          
          id++;
        }
      }
    }

  }

  return id;
}

      









