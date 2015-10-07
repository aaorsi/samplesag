#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
//#include "defs.h"
//#include "functions.h"

double comoving_distance(double z)
{
  double dist;
  double integrand2(double a);
  double dqromb(double (*func)(double), double a, double b);
  
  double UnitLength_in_cm = 3.085678e21;
  double UnitMass_in_g = 1.989e43; 
  double UnitVelocity_in_cm_per_s = 1e5;
  double UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;
  double UnitTime_in_Megayears = UnitTime_in_s/SEC_PER_MEGAYEAR;
//  double Hubble = HUBBLE * UnitTime_in_s;
  //double dh = C_KM_S/H0_KM_S_MPC; // Para distancia [Mpc]
  double dh = C_KM_S/100.0; // Para distancia [Mpc/h]
  
  dist= dh*dqromb(integrand2, 1/(z+1), 1);

  return dist;
}

double integrand2(double a)
{
  double Omega = OmegaM;
  double OmegaLambda = 1 - OmegaM;;
  
  return 1/sqrt(Omega*a + (1-Omega-OmegaLambda)*a*a + OmegaLambda*a*a*a*a);
}



