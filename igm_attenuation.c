#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <math.h>
#include "allvars.h"
#include <hdf5.h>
#include <hdf5_hl.h>
#include "proto.h"

// This routine computs the IGM attenuation from Madau 95.

void igm_attenuation(float *tauArr, float *wlambda,  float z, int id)
{
	int i;
	float taueff;
	float z1 = 1 + z;
	float spec;
	float thex;
	for (i = 0; i<NBinSed; i++)
	{
		taueff = 0.0;
		spec = wlambda[i];
		if (1026.0 * z1 <= spec && spec <= 1216.0*z1)
			taueff = 0.0036 * pow(spec/1216.0,3.46);

		if (973.0 * z1 <= spec && spec <= 1026.0*z1)
			taueff = 0.0036 * pow(spec/1216.0,3.46) \
					  + 0.0017 * pow(spec/1026.0,3.46);

		if (950.0 * z1 <= spec && spec <= 973.0*z1)
			taueff = 0.0036 * pow(spec/1216.0,3.46) \
					  + 0.0017 * pow(spec/1026.0,3.46) \
					  + 0.0012 * pow(spec/973.0,3.46);
		
		if (912.0*z1 <= spec && spec <= 950.0)
			taueff = 0.0036 * pow(spec/1216.0,3.46) \
					  + 0.0017 * pow(spec/1026.0,3.46) \
					  + 0.0012 * pow(spec/973.0,3.46)  \
					  + 0.00093* pow(spec/950.0,3.46);
	 	
		if (912.0 <= spec && spec <= 912.0*z1)
		{
			thex = spec / 912.0;
			taueff = 0.0036 * pow(spec/1216.0,3.46) \
					  + 0.0017 * pow(spec/1026.0,3.46) \
					  + 0.0012 * pow(spec/973.0,3.46)  \
					  + 0.00093* pow(spec/950.0,3.46)  \
					  + 0.25*pow(thex,3)*(pow(z1,0.46) - pow(thex,0.46))  \
					  + 9.4*pow(thex,1.5)*(pow(z1,0.18) - pow(thex,0.18)) \
					  + 0.7*pow(thex,3)*(pow(z1,-1.32) - pow(thex,-1.32)) \
					  - 0.023*(pow(z1,1.68) - pow(thex,1.68));
		}
	

		tauArr[i] = taueff;

	}
	
	
}			

