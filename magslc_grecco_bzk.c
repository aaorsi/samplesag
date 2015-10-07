#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define NZ 98// Numero de snapshots
#define SNAP1 14 // Primeras galaxias
#define SNAP0 17
// Magnitud de Vega en AB
#define UXVEGA 0.810
#define BXVEGA -0.100
#define BVEGA -0.098
#define VVEGA 0.012
#define RVEGA 0.245
#define IVEGA 0.489
#define JVEGA 0.891
#define HVEGA 1.369
#define KVEGA 1.875
#define LVEGA 2.751
#define MVEGA 3.376
#define BWVEGA 0.0 // http://www.noao.edu/noao/noaodeep/
#define ZVEGA 0.54 // Taylor et al. (2009)
//#define ZVEGA 0.533 // Hewett et al. (2006)
#define KSVEGA 1.857 // Takagi et al. (2008)

#define VAR_IMF 0

int main(int argc, char *argv[])
{
  FILE *salida1, *salida2, *salida3, *salida4, *salida5, *salida6, *salida7, *salida8, *salida9, *salida10, *salida11;
  float UX, BX, B, V, R, I, J, H, K, L, M, Bw, z, Ks; // Magnitudes en AB
  float f75, f110, f170, f450, f850, f870, f21cm; // Flujos [mJy]
  float redshift;
  int nt[NZ]; // Numero de galaxias
  int nv[NZ]; // Numero de galaxias con StellarMass >= 5 M_Sun
  int ns[NZ]; // Numero de galaxias de muestra
  int nr[NZ];
  int i, j, k, n, r;
  
  k = 0;
  for(r=SNAP0-1;r<NZ;r++) {
    char snap[256];
    char dir[256];
    char dir2[256];
    char fname1[256];
    char fname2[256];
    char fname3[256];
    char fname4[256];
    char fname5[256];
    char fname6[256];
    char fname7[256];
    char fname10[256];
    char fname11[256];
    char imf[256];
    
    if(VAR_IMF==0) sprintf(imf,"SalpIMF");
    else sprintf(imf,"VarIMF");
        
    snprintf(snap,3,"%d",r+1);
    sprintf(dir,"/data/ammunoz1/GRASIL/Input/GRECCO/Lightcone_%s/Grasil_SAG3%s_Snap%s",imf,imf,snap);
    sprintf(dir2,"/data/ammunoz1/GRASIL/Pruebas/GRECCO/Lightcone_%s_Fluxes/Grasil_SAG3%s_Snap%s",imf,imf,snap);
    //sprintf(dir2,"/data/ammunoz1/GRASIL/Pruebas/GRECCO/Lightcone_%s_Fluxes_IGM/Grasil_SAG3%s_Snap%s",imf,imf,snap); // IGM

    sprintf(fname1,"%s/NumGalaxies.dat",dir);
    sprintf(fname2,"%s/magsabs.dat",dir2);
    sprintf(fname3,"%s/Grasil_SAG3%s_idgals.list",dir,imf);
    sprintf(fname4,"%s/magsunknown.dat",dir2);
    sprintf(fname5,"%s/fluxes_850_870.dat",dir2);
    //sprintf(fname5,"%s/fluxes_870.dat",dir2);
    sprintf(fname6,"%s/magsap.dat",dir2);
    sprintf(fname7,"%s/Grasil_SAG3%s_origs.list",dir,imf);
    //sprintf(fname10,"%s/magscolor.dat",dir2);
    //sprintf(fname11,"%s/mags_bw_z_ks.dat",dir2);
    //sprintf(fname11,"%s/mags_bzk.dat",dir2);
            
    salida1 = fopen(fname1,"r");
    salida2 = fopen(fname2,"w");
    salida3 = fopen(fname3,"r");
    salida4 = fopen(fname4,"w");
    salida5 = fopen(fname5,"w");
    salida6 = fopen(fname6,"w");
    salida7 = fopen(fname7,"r");
    //salida10 = fopen(fname10,"w");
    //salida11 = fopen(fname11,"w");
    
    fscanf(salida1,"%d %d %d %d",&nt[r],&nv[r],&ns[r],&nr[r]);
    
    j = 0;
    for(i=0;i<nr[r];i++) {
      char gal[256];
      char str[256];
      char str2[300];
      char fname8[256];
      char fname9[256];
      int idorig;
      
      fscanf(salida3,"%s",gal);
      fscanf(salida7,"%d",&idorig);
      sprintf(fname8,"%s/col%s",dir,gal);
      sprintf(fname9,"%s/colz%s",dir,gal);
      //sprintf(fname8,"%s/colspm%s",dir2,gal); // IGM
      //sprintf(fname9,"%s/colzspm%s",dir2,gal); // IGM
                  
      salida8 = fopen(fname8,"r");
      fgets(str,sizeof(str),salida8);
      fscanf(salida8,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f", &UX, &BX, &B, &V, &R, &I, &J, &H, &K, &L, &M, &Bw, &z, &Ks);
      
      if(UX<90.0) { // Para considerar solo flujos validos
        j++;
        k++;
	
	fprintf(salida2,"%s %d %d %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n", gal, i+1, idorig, UX+UXVEGA, BX+BXVEGA, B+BVEGA, V+VVEGA, R+RVEGA, I+IVEGA, J+JVEGA, H+HVEGA, K+KVEGA, L+LVEGA, M+MVEGA, Bw+BWVEGA, z+ZVEGA, Ks+KSVEGA); // Magnitudes absolutas en AB
	//fprintf(salida10,"%s %.3f %.3f\n", gal, V+VVEGA, (B+BVEGA)-(V+VVEGA));
	
	salida9 = fopen(fname9,"r");
	
	for(n=0;n<4;n++) fgets(str2,sizeof(str2),salida9);
	fscanf(salida9,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f", &redshift, &UX, &B, &V, &R, &I, &J, &H, &K, &L, &M, &Bw, &z, &Ks, &f75, &f110, &f170, &f450, &f850, &f870, &f21cm);
	fprintf(salida5,"%d %d %.3e %.3e\n", i+1, idorig, f850, f870); // Id original de galaxias (SAM, stellarMass >= 5 M_Sun) y flujo a 870 micrones [mJy]
	
	fprintf(salida6,"%s %d %d %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n", gal, i+1, idorig, redshift, UX+UXVEGA, B+BVEGA, V+VVEGA, R+RVEGA, I+IVEGA, J+JVEGA, H+HVEGA, K+KVEGA, L+LVEGA, M+MVEGA, Bw+BWVEGA, z+ZVEGA, Ks+KSVEGA, f75, f110, f170, f450, f850, f870, f21cm); // Magnitudes aparentes en AB y flujos [mJy]
	//fprintf(salida11,"%d %d %.3f %.3f %.3f\n", i+1, idorig, Bw+BWVEGA, z+ZVEGA, Ks+KSVEGA); // Magnitudes aparentes en AB y flujos [mJy]
	//fprintf(salida11,"%d %d %.3f %.3f %.3f\n", i+1, idorig, B+BVEGA, z+ZVEGA, K+KVEGA); // Magnitudes aparentes en AB y flujos [mJy]
	
	fclose(salida9);
      }
      else {
        fprintf(salida4,"%d %s\n",i+1,gal); // Galaxias con flujo no valido
      }
      
      fclose(salida8);
    }
    
    printf("%d Galaxies: total = %d valid = %d sample = %d ids = %d valid mags = %d\n",r+1,nt[r],nv[r],ns[r],nr[r],j);
    
    fclose(salida1);
    fclose(salida2);
    fclose(salida3);
    fclose(salida4);
    fclose(salida5);
    fclose(salida6);
    fclose(salida7);
    //fclose(salida10);
    //fclose(salida11);
  }
    
  printf("Valid magnitudes = %d\n",k);
  
  return(1);
}
