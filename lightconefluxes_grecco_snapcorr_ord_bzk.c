#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "lib/nrsrc/nrutil.h"
#include "proto.h"

#include "allvars.h"

#define C_KM_S 2.9979e5
#define H0_KM_S_MPC 70.0
#define BOX 150.0 // Lado del cubo de simulacion [Mpc/h]
#define NZ 99 // Numero de snapshots
#define NZ0 98
#define SNAP1 14 // Primeras galaxias
#define SNAP0 17
#define LC 62.16
#define NCX1 3
#define NCX2 4
/*#define LC 124.4
#define NCX1 4
#define NCX2 3*/ // Para area 1 sq deg
#define NDCZ 9800

#define VAR_IMF 0

//USAR ESTE!!!

int iddcgal(double *xx, int n, double x)
{
  int j, jl, ju, jm;

  jl = 0;
  ju = n+1;

  while((ju-jl)>1) {
    jm = (ju+jl)/2.0;
    
    if((xx[n-1]>xx[0]) && (x>xx[jm-1])) jl = jm;
    else ju = jm;
  }

  j = jl;

  if(j>0 && j<n) {
    if(fabs(x-xx[j-1])>fabs(x-xx[j])) j = j+1;
  }

  if(j==0) j = j+1;

  return j;
}

double comoving_distance(double z)
{
  double dist;
  double integrand2(double a);
  
  double UnitLength_in_cm = 3.085678e21;
  double UnitMass_in_g = 1.989e43; 
  double UnitVelocity_in_cm_per_s = 1e5;
  double UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;
  double UnitTime_in_Megayears = UnitTime_in_s/SEC_PER_MEGAYEAR;
  double Hubble = HUBBLE * UnitTime_in_s;
  
  //double dh = C_KM_S/H0_KM_S_MPC; // Para distancia [Mpc]
  double dh = C_KM_S/100.0; // Para distancia [Mpc/h]
  
  dist= dh*dqromb(integrand2, 1/(z+1), 1);

  return dist;
}

double integrand2(double a)
{
  double Omega = 0.28;
  double OmegaLambda = 0.72;
  
  return 1/sqrt(Omega*a + (1-Omega-OmegaLambda)*a*a + OmegaLambda*a*a*a*a);
}

int main(void)
{
  FILE *os, *totals, *catalog, *catalog2, *p, *ngalsnap, *numgal, *magsap, *magsabs, *dcs, *mvir;
  int i, j, k, q, s, t, x;
  int totcube[NZ+1];
  int ncubex1[NZ+1];
  int ncubex2[NZ+1];
  int ncubex3[NZ+1];
  int nt[NZ+1]; // Numero de galaxias
  int nv[NZ+1]; // Numero de galaxias con StellarMass >= 5 M_Sun
  int ns[NZ+1]; // Numero de galaxias de muestra
  int nr[NZ+1];
  double eta;
  double theta;
  double epsilon;
  double etap;
  double epsilonp;
  double area;
  double ddccube;
  double ddccube0;
  double aext[NZ+1]; // Factor de escala
  double redshift[NZ+1]; // Redshift
  double dcbase[NZ+1]; // Distancia comovil [Mpc/h]
  double dabase[NZ+1]; // Distancia diametro angular [Mpc/h]
  double ddc[NZ+1];
  double letabase[NZ+1];
  double lthetabase[NZ+1];
  double lepsilonbase[NZ+1];
  double ztable[NDCZ];
  double dctable[NDCZ];
  double pos[3]; // Posicion original de galaxia [Mpc/h]
  double *f870, *mag_B, *mag_K, *mag_i, *mag_z, *mag_Kabs, *mvir_halofof, *mvir_subest;
  char dir[256];
  char dir2[256];
  char dir3[256];
  char fname1[256];
  char fname2[256];
  char fname3[256];
  char fname8[256];
  char fname11[256];
  char imf[256];
  
  if(VAR_IMF==0) sprintf(imf,"SalpIMF");
  else sprintf(imf,"VarIMF");
  
  sprintf(dir,"/data/ammunoz1/GRASIL/Input/GRECCO/Lightcone_%s/",imf);
  sprintf(dir2,"/data/ammunoz1/GRASIL/Pruebas/GRECCO/Lightcone_%s_Catalog/",imf);
  sprintf(dir3,"/data/ammunoz1/GRASIL/Pruebas/GRECCO/Lightcone_%s_Fluxes/",imf);
  //sprintf(dir3,"/data/ammunoz1/GRASIL/Pruebas/GRECCO/Lightcone_%s_Fluxes_IGM/",imf); // IGM
  
  sprintf(fname1,"%stotals_all.dat",dir);
  sprintf(fname2,"%scatalogfluxes_bzk.dat",dir2);
  //sprintf(fname2,"%scatalogfluxes_bzk_igm.dat",dir2); // IGM
  sprintf(fname3,"%sngalsnapfluxes_bzk.dat",dir2);
  //sprintf(fname3,"%sngalsnapfluxes_bzk_igm.dat",dir2); // IGM
  sprintf(fname8,"%scatalogfluxesrandom_bzk.dat",dir2);
  //sprintf(fname8,"%scatalogfluxesrandom_bzk_igm.dat",dir2); // IGM
  sprintf(fname11,"/home/ammunoz1/SUBMM2011/Comoving_Distance_GRECCO/dcz.dat");
  
  os = fopen("/data/PruebasSAM/d2/Clusters/GRECCO/output_list.dat","r"); // Factor de escala para cada snapshot
  totals = fopen(fname1,"r");
  catalog = fopen(fname2,"w");
  ngalsnap = fopen(fname3,"w");
  catalog2 = fopen(fname8,"w");
  dcs = fopen(fname11,"r");
 
  i = 0;
  while(fscanf(os,"%lf",&aext[i])!= EOF) {
    redshift[i] = 1/aext[i]-1;
    dcbase[i] = comoving_distance(redshift[i]);
    dabase[i] = dcbase[i]*aext[i];
    i++;
  }
  fclose(os);
  
  eta = NCX2*BOX/(dabase[SNAP1]*(1.0+redshift[SNAP1]));
  theta = LC/(dabase[SNAP1]*(1.0+redshift[SNAP1]));
  epsilon = NCX1*BOX/(dabase[SNAP1]*(1.0+redshift[SNAP1]));
  etap = -1.0*(eta-0.5*theta);
  epsilonp = -1.0*(epsilon-0.5*theta);
  
  printf("LC = %f Mpc/h, SNAP1 = %d, theta = %f rad = %f deg = %f arcmin, area = %f sq deg\n",LC,SNAP1,theta,(180/PI)*theta,(60*180/PI)*theta,(180/PI)*theta*(180/PI)*theta);
  
  i = SNAP1;
  while(fscanf(totals,"%d %d",&x,&totcube[i])!= EOF) {
    i++;    
  }
  fclose(totals);
  
  i = NDCZ-1;
  while(fscanf(dcs,"%d %lf %lf",&x,&ztable[i],&dctable[i])!= EOF) {
    i--;
  }
  fclose(dcs);
  
  q = 0;
  ddccube0 = 0.0;
  
  for(i=0;i<NZ;i++) nv[i] = 0;
  
  for(i=SNAP1;i<SNAP0;i++) {
    ddc[i] = dcbase[i]-dcbase[i+1];
    ncubex3[i] = (int)((ddc[i]+ddccube0)/BOX)+1;
    ddccube = ddc[i]+ddccube0-(ncubex3[i]-1)*BOX;
    ddccube0 = ddccube;
  }
  
  for(i=SNAP0;i<=NZ0;i++) {
    char snap[256];
    char snap_next[256];    
    char fname4[256];
    char fname5[256];
    char fname6[256];
    char fname7[256];
    char fname9[256];
    char fname10[256];
    char fname12[256];
    char caux[256];
    double x1;
    double x2;
    double x3; // Coordenada comovil x_3 de galaxia [Mpc/h]
    double dc; // Distancia comovil de galaxia [Mpc/h]
    double da; // Distancia diametro angular de galaxia [Mpc/h]
    double leta;
    double ltheta;
    double lepsilon;
    double x3max; // Maxima coordenada comovil x_3 para ultimo cubo del snapshot [Mpc/h]
    double aux, aux2;
    double f870_aux, mag_B_aux, mag_I_aux, mag_K_aux, mag_z_aux, mag_Kabs_aux;
    double mvir_halofof_aux, mvir_subest_aux;
    int id; // Id original de galaxia (SAM, StellarMass >= 5 M_Sun)
    int idz;
    int n, idorig;
    double alpha; // Ascencion recta de galaxia [arcsec]
    double delta; // Declinacion de galaxia [arcsec]
    double x1p, x2p, x3p, x1pp, x2pp, x3pp;
    double dcnew;
    double alphar;
    double deltar;
    double x1r;
    double x2r;
    double dcr;
    double zgal, zgalr;
    double x1min, x1max, x2min, x2max, x1minp, x1maxp, x2minp, x2maxp, x3minp, x3maxp, x1minpp, x1maxpp, x2minpp, x2maxpp, x3minpp, x3maxpp;
    
    snprintf(snap,3,"%d",i);
    snprintf(snap_next,3,"%d",i+1);
    sprintf(fname4,"%sGrasil_SAG3%s_Snap%s/allpos.dat",dir,imf,snap);
    sprintf(fname5,"%sGrasil_SAG3%s_Snap%s/NumGalaxies.dat",dir,imf,snap);
    sprintf(fname6,"%sGrasil_SAG3%s_Snap%s/allpos.dat",dir,imf,snap_next);
    sprintf(fname9,"%sGrasil_SAG3%s_Snap%s/magsap.dat",dir3,imf,snap);
    sprintf(fname10,"%sGrasil_SAG3%s_Snap%s/magsabs.dat",dir3,imf,snap);
    sprintf(fname12,"%sGrasil_SAG3%s_Snap%s/mvir.dat",dir,imf,snap);
    
    printf("Snapshot %d\n",i);

    if(i<NZ) ddc[i] = dcbase[i]-dcbase[i+1];
    if(i==NZ) ddc[i] = dcbase[i];
    
    ncubex1[i] = NCX1;
    ncubex2[i] = NCX2;
    ncubex3[i] = (int)((ddc[i]+ddccube0)/BOX)+1;
    
    ddccube = ddc[i]+ddccube0-(ncubex3[i]-1)*BOX;
    x3max = ddc[i]+ddccube0;
    
    letabase[i] = eta*dcbase[i];
    lthetabase[i] = theta*dcbase[i];
    lepsilonbase[i] = epsilon*dcbase[i];
    
    numgal = fopen(fname5,"r");
    magsap = fopen(fname9,"r");
    magsabs = fopen(fname10,"r");
    mvir = fopen(fname12,"r");
    
    fscanf(numgal,"%d %d %d %d",&nt[i],&nv[i],&ns[i],&nr[i]);
    //printf("%d %d %d %d\n",nt[i],nv[i],ns[i],nr[i]);
    fclose(numgal);
    
    f870 = (double*)malloc(nv[i]*sizeof(double));
    mag_B = (double*)malloc(nv[i]*sizeof(double));
    mag_K = (double*)malloc(nv[i]*sizeof(double));
    mag_i = (double*)malloc(nv[i]*sizeof(double));
    mag_z = (double*)malloc(nv[i]*sizeof(double));
    mag_Kabs = (double*)malloc(nv[i]*sizeof(double));
    mvir_halofof = (double*)malloc(nv[i]*sizeof(double));
    mvir_subest = (double*)malloc(nv[i]*sizeof(double));
    
    for(x=0;x<nv[i];x++) f870[x] = -1.0; // Para invalidar galaxias con flujo erroneo

    while(fscanf(magsap,"%s %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&caux,&x,&idorig,&aux,&aux,&mag_B_aux,&aux,&aux,&mag_I_aux,&aux,&aux,&mag_K_aux,&aux,&aux,&aux,&mag_z_aux,&aux,&aux,&aux,&aux,&aux,&aux,&f870_aux,&aux)!= EOF) {
    //while(fscanf(magsap,"%s %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&caux,&x,&idorig,&aux,&aux,&aux,&aux,&aux,&aux,&aux,&aux,&aux,&aux,&aux,&aux,&aux,&aux,&aux,&aux,&aux,&aux,&aux,&f870_aux,&aux)!= EOF) {
      f870[idorig] = f870_aux; // Asignar flujo a galaxias con id original
      mag_B[idorig] = mag_B_aux;
      mag_K[idorig] = mag_K_aux;
      mag_z[idorig] = mag_z_aux;
      mag_i[idorig] = (mag_I_aux-0.3780*mag_z[idorig]+0.3974)/(1-0.3780); // Lupton (2005): I = i-0.3780*(i-z)-0.3974 (http://www.sdss.org/dr5/algorithms/sdssUBVRITransform.html)
    }
    fclose(magsap);
    
    while(fscanf(magsabs,"%s %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&caux,&x,&idorig,&aux,&aux,&aux,&aux,&aux,&aux,&aux,&aux,&mag_Kabs_aux,&aux,&aux,&aux,&aux,&aux)!= EOF) {
      mag_Kabs[idorig] = mag_Kabs_aux;
    }
    fclose(magsabs);

    while(fscanf(mvir,"%d %d %d %d %lf %lf",&idorig,&x,&x,&x,&mvir_halofof_aux,&mvir_subest_aux)!= EOF) {
      mvir_halofof[idorig] = mvir_halofof_aux;
      mvir_subest[idorig] = mvir_subest_aux;
    }
    fclose(mvir);

    n = 0;

    for(s=1;s<=ncubex1[i];s++) {
    for(t=1;t<=ncubex2[i];t++) {
    for(j=1;j<=ncubex3[i];j++) {
      p = fopen(fname4,"r");

      for(k=0;k<totcube[i];k++) {
	fscanf(p,"%d %lf %lf %lf",&id,&pos[0],&pos[1],&pos[2]);
	
	x1 = pos[0]+(s-1)*BOX;
	x2 = pos[1]+(t-1)*BOX;
	x3 = pos[2]+(j-1)*BOX;
	dc = dcbase[i]+ddccube0-x3;
	da = dc*aext[i];
	leta = eta*dc;
	ltheta = theta*dc;
	lepsilon = epsilon*dc;
	
	// Obtener ascension recta y declinacion para galaxias dentro del beam y con flujo valido
	if(x3>=ddccube0 && x3<x3max && dc>0.0) {
	  if(x1>lepsilon-ltheta && x1<lepsilon) {
	    if(x2>leta-ltheta && x2<leta) {
	      if(f870[id]>-1.0) {
		x1p = x1;
		x2p = cos(etap)*x2+sin(etap)*dc;
		x3p = -1.0*sin(etap)*x2+cos(etap)*dc;
		x1pp = cos(epsilonp)*x1p+sin(epsilonp)*x3p;
		x2pp = x2p;
		x3pp = -1.0*sin(epsilonp)*x1p+cos(epsilonp)*x3p;
		
		alpha = 3600.0*(180/PI)*(x1pp/x3pp);
		delta = 3600.0*(180/PI)*(x2pp/x3pp);
		dcnew = x3pp;
		
		  n++;
		  q++;
		  
		  x1min = lepsilon-ltheta;
		  x1max = lepsilon;
		  x2min = leta-ltheta;
		  x2max = leta;
		  
		  x1minp = x1min;
		  x2minp = cos(etap)*x2min+sin(etap)*dc;
		  x3minp = -1.0*sin(etap)*x2min+cos(etap)*dc;
		  x1minpp = cos(epsilonp)*x1minp+sin(epsilonp)*x3minp;
		  x2minpp = x2minp;
		  x3minpp = -1.0*sin(epsilonp)*x1minp+cos(epsilonp)*x3minp;
		  
		  x1maxp = x1max;
		  x2maxp = cos(etap)*x2max+sin(etap)*dc;
		  x3maxp = -1.0*sin(etap)*x2max+cos(etap)*dc;
		  x1maxpp = cos(epsilonp)*x1maxp+sin(epsilonp)*x3maxp;
		  x2maxpp = x2maxp;
		  x3maxpp = -1.0*sin(epsilonp)*x1maxp+cos(epsilonp)*x3maxp;
		  
		  x1pp = (x1maxpp-x1minpp)*((float)rand()/(float)RAND_MAX)+x1minpp;
		  x2pp = (x2maxpp-x2minpp)*((float)rand()/(float)RAND_MAX)+x2minpp;
		  
		  alphar = 3600.0*(180/PI)*(x1pp/x3pp);
		  deltar = 3600.0*(180/PI)*(x2pp/x3pp);
		  dcr = x3pp;
		  
		  idz = iddcgal(dctable,NDCZ,dcnew)-1;
		  zgal = ztable[idz];
		  idz = iddcgal(dctable,NDCZ,dcr)-1;
		  zgalr = ztable[idz];
		  
		  fprintf(catalog,"%d %d %d %lf %lf %lf %.3e %.3lf %.3lf %.3lf %.3lf %.3lf %lf %.6e %.6e\n",i,q,id,alpha,delta,dcnew,f870[id],mag_B[id],mag_z[id],mag_K[id],mag_i[id],mag_Kabs[id],zgal,mvir_halofof[id],mvir_subest[id]);
		  fprintf(catalog2,"%d %d %d %lf %lf %lf %.3e %.3lf %.3lf %.3lf %.3lf %.3lf %lf %.6e %.6e\n",i,q,id,alphar,deltar,dcr,f870[id],mag_B[id],mag_z[id],mag_K[id],mag_i[id],mag_Kabs[id],zgalr,mvir_halofof[id],mvir_subest[id]);
		  //fprintf(catalog,"%d %d %d %lf %lf %lf %.3e %lf\n",i,q,id,alpha,delta,dcnew,f870[id],zgal);
		  //fprintf(catalog2,"%d %d %d %lf %lf %lf %.3e %lf\n",i,q,id,alphar,deltar,dcr,f870[id],zgalr);
	      }
	    }
	  }
	}
      }
      fclose(p);
    }
    }
    }
    
    ddccube0 = ddccube;
    
    if(i<NZ) {
      int id_next;
      double pos_next[3];

      ddc[i+1] = dcbase[i+1]-dcbase[i+2];
    
      ncubex1[i+1] = NCX1;
      ncubex2[i+1] = NCX2;
      ncubex3[i+1] = (int)((ddc[i+1]+ddccube0)/BOX)+1;

      ddccube = ddc[i+1]+ddccube0-(ncubex3[i+1]-1)*BOX;
      x3max = ddc[i+1]+ddccube0;

      letabase[i+1] = eta*dcbase[i+1];
      lthetabase[i+1] = theta*dcbase[i+1];
      lepsilonbase[i+1] = epsilon*dcbase[i+1];

    for(s=1;s<=ncubex1[i+1];s++) {
    for(t=1;t<=ncubex2[i+1];t++) {
    for(j=1;j<=ncubex3[i+1];j++) {
      p = fopen(fname6,"r");

      for(k=0;k<totcube[i+1];k++) {
        fscanf(p,"%d %lf %lf %lf",&id_next,&pos_next[0],&pos_next[1],&pos_next[2]);

        x1 = pos_next[0]+(s-1)*BOX;
        x2 = pos_next[1]+(t-1)*BOX;
        x3 = pos_next[2]+(j-1)*BOX;
        dc = dcbase[i+1]+ddccube0-x3;
        da = dc*aext[i+1];
        leta = eta*dc;
        ltheta = theta*dc;
        lepsilon = epsilon*dc;

        // Obtener ascension recta y declinacion para galaxias dentro del beam y con flujo valido
        if(x3>=ddccube0 && x3<x3max && dc>0.0) {
          if(x1>lepsilon-ltheta && x1<lepsilon) {
            if(x2>leta-ltheta && x2<leta) {
	      if(id_next<nv[i] && f870[id_next]>-1.0) {
		x1p = x1;
		x2p = cos(etap)*x2+sin(etap)*dc;
		x3p = -1.0*sin(etap)*x2+cos(etap)*dc;
		x1pp = cos(epsilonp)*x1p+sin(epsilonp)*x3p;
		x2pp = x2p;
		x3pp = -1.0*sin(epsilonp)*x1p+cos(epsilonp)*x3p;
		
		alpha = 3600.0*(180/PI)*(x1pp/x3pp);
		delta = 3600.0*(180/PI)*(x2pp/x3pp);
		dcnew = x3pp;
		
		if(dcnew-dcbase[i+1]>0.0) {
		
		  if(id_next<totcube[i]) {
		    n++;
		    q++;
		    
		    x1min=lepsilon-ltheta;
		    x1max=lepsilon;
		    x2min=leta-ltheta;
		    x2max=leta;
		    
		    x1minp = x1min;
		    x2minp = cos(etap)*x2min+sin(etap)*dc;
		    x3minp = -1.0*sin(etap)*x2min+cos(etap)*dc;
		    x1minpp = cos(epsilonp)*x1minp+sin(epsilonp)*x3minp;
		    x2minpp = x2minp;
		    x3minpp = -1.0*sin(epsilonp)*x1minp+cos(epsilonp)*x3minp;
		    
		    x1maxp = x1max;
		    x2maxp = cos(etap)*x2max+sin(etap)*dc;
		    x3maxp = -1.0*sin(etap)*x2max+cos(etap)*dc;
		    x1maxpp = cos(epsilonp)*x1maxp+sin(epsilonp)*x3maxp;
		    x2maxpp = x2maxp;
		    x3maxpp = -1.0*sin(epsilonp)*x1maxp+cos(epsilonp)*x3maxp;
		    
		    x1pp = (x1maxpp-x1minpp)*((float)rand()/(float)RAND_MAX)+x1minpp;
		    x2pp = (x2maxpp-x2minpp)*((float)rand()/(float)RAND_MAX)+x2minpp;
		    
		    alphar = 3600.0*(180/PI)*(x1pp/x3pp);
		    deltar = 3600.0*(180/PI)*(x2pp/x3pp);
		    dcr = x3pp;
		    
		    idz = iddcgal(dctable,NDCZ,dcnew)-1;
		    zgal = ztable[idz];
		    idz = iddcgal(dctable,NDCZ,dcr)-1;
		    zgalr = ztable[idz];
		    
		    fprintf(catalog,"%d %d %d %lf %lf %lf %.3e %.3lf %.3lf %.3lf %.3lf %.3lf %lf %.6e %.6e\n",i,q,id_next,alpha,delta,dcnew,f870[id_next],mag_B[id_next],mag_z[id_next],mag_K[id_next],mag_i[id_next],mag_Kabs[id_next],zgal,mvir_halofof[id_next],mvir_subest[id_next]);
		    fprintf(catalog2,"%d %d %d %lf %lf %lf %.3e %.3lf %.3lf %.3lf %.3lf %.3lf %lf %.6e %.6e\n",i,q,id_next,alphar,deltar,dcr,f870[id_next],mag_B[id_next],mag_z[id_next],mag_K[id_next],mag_i[id_next],mag_Kabs[id_next],zgalr,mvir_halofof[id_next],mvir_subest[id_next]);
		    //fprintf(catalog,"%d %d %d %lf %lf %lf %.3e %lf\n",i,q,id_next,alpha,delta,dcnew,f870[id_next],zgal);
		    //fprintf(catalog2,"%d %d %d %lf %lf %lf %.3e %lf\n",i,q,id_next,alphar,deltar,dcr,f870[id_next],zgalr);
		  }
                }
              }
            }
          }
        }
      }
      fclose(p);
    }
    }
    }

    }

    fprintf(ngalsnap,"%d %d\n",i,n);

    printf("snap %d: n=%d q=%d\n",i,n,q);
    
    free(f870);
    free(mag_B);
    free(mag_z);
    free(mag_K);
    free(mag_Kabs);
  }
  
  printf("Total galaxies = %d\n",q);
  
  fclose(catalog);
  fclose(catalog2);
  fclose(ngalsnap);
}


