
void get_filterinfo(char *);
void read_filterfile(char *);
float trapz1(float *,float *,int );
float compute_mag_filter(int ,int ,float *,float *,int ,float ,int ,int );
void get_Vega_AB(char *);
void read_saghdf5(char *, sagobj *, float *, long );
int set_float_variable_attr(hid_t , char *, char *, double *, char *);
void igm_attenuation(float *, float *, float , int );
void write_hdf5(props *, float *, float *,float *,char *);
void write_fits(props *, float *, float *,float *,char *);
void printerror(int );
