
#define dprinti(expr) printf(#expr " = %d\n", expr)
#define dprintf(expr) printf(#expr " = %f\n", expr)
#define dprints(expr) printf(#expr " = %s\n", expr)
#define MAXLEN 500				// Max. length of strings
#define MAXMAG 10				// Max. number of bands
#define MAXGALS 20000000		// Max. Number of gals. in Lightcone
#define MAXZ 2000				// Max. n. of rows in lookup table
#define MAXGALSNAP 5000000		// Max. N' of gals. in a snapshot
// Cosmological parameters
#define OmegaM 0.28
#define H0 100.0
#define Hubble_h 0.7
#define C_KM_S 2.9979e5

// SED and filter definitions
#define SED_RV       3.1
#define MAKE_FD_DATA 0
#define NFILTERMAX   300
#define NFILTERL     10000
#define SED_MAXC     500

#define SED_MAXW      2000
#define SED_MAXAGE    300
#define SED_MAXCHAR   150000
#define SED_NBANDSMAX 30

#define TAUBSTAR 0.8
#define BETA 0.5
/*#define Lbstar 3.7e8*/ /*In solar luminosity is 5.75e10 */
/* Si uso MB=19.6+5logh como en De Lucia et al (2004) es (27/12/06):*/
#define LBSTAR 1.4e8 /*In solar luminosity is 5.75e10 */
#define MU_DUST (1/3.)

/* Extinction curve of Cardelli et al. (1989) */ 
#define Au_Av 1.569
#define Ab_Av 1.337
#define Ar_Av 0.751
#define Ai_Av 0.479
#define Ak_Av 0.114




typedef struct {
	int snapshot;	// galaxy snapshot 
	long id_lc;		// id within lightcone
	int  id_snap;	// id within snapshot file

	float x;		// cartesian coordinates
	float y;
	float z;
	
	float rho;		//spherical coordinates
	float theta;
	float phi;

	float ra;		// angular coordinates
	float dec;
	
	float redshift;	// redshift in lightcone
	float v_r;		// radial velocity
} props;

typedef struct {
	int id;
	float *Sed;
	float StellarMass;
} sagobj;
	
int NBinSed;
int NMags;
long NGals;
int NBoxes;
long NTotGals;
int * Nsteps_filter;
int idMag[MAXMAG];
char NameMags[MAXMAG][16];
char *Idinfo[SED_NBANDSMAX];
float * Filterlist;
float MVega[SED_NBANDSMAX];
float MinStellarMass;
