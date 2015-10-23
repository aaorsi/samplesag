#include <stdlib.h>
#include <stdio.h>
#include "allvars.h"

void read_args(int argc, char *argv)
{

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

}  
