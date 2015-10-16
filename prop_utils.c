

int get_magsprops(char * PropArr)
{
  int dust_check,mgal;

// This function identifies which magnitudes will be computed. Returns 1 if there are magnitudes at all, 0 otherwise.
//  ** Magnitudes **    
  dust_check = 0;
  mgal = 0;
  for (j = 0; j< PropLength; j++)
  { 
    // Check if PropArr is a magnitude
    if (strstr(PropArr[j],"mag_") != 0)
    {
      mgal = 1;
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
        }
      }
    }
  }
  return mgal;

}
